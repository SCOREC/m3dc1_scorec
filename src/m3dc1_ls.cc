/****************************************************************************** 

  (c) 2005-2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef M3DC1_MESHGEN
#ifdef M3DC1_TRILINOS
#include "m3dc1_ls.h"
#include "m3dc1_scorec.h"
#include "apf.h"
#include "apfNumbering.h"
#include "apfShape.h"
#include "apfMesh.h"
#include "PCU.h"
#include "m3dc1_mesh.h"
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

// ***********************************
// 		M3DC1_LINEAR SYSTEM
// ***********************************

m3dc1_ls* m3dc1_ls::_instance=NULL;
m3dc1_ls* m3dc1_ls::instance()
{
  if (_instance==NULL)
    _instance = new m3dc1_ls();
  return _instance;
}

m3dc1_ls::~m3dc1_ls()
{
  if (matrix_container!=NULL)
    matrix_container->clear();
  matrix_container=NULL;
  delete _instance;
}

void m3dc1_ls::add_matrix(int matrix_id, m3dc1_epetra* matrix)
{
  assert(matrix_container->find(matrix_id)==matrix_container->end());
  matrix_container->insert(std::map<int, m3dc1_epetra*>::value_type(matrix_id, matrix));
}

m3dc1_epetra* m3dc1_ls::get_matrix(int matrix_id)
{
  std::map<int, m3dc1_epetra*>::iterator mit = matrix_container->find(matrix_id);
  if (mit == matrix_container->end()) 
    return (m3dc1_epetra*)NULL;
  return mit->second;
}


//*******************************************************
apf::Numbering* get_owned_numbering()
//*******************************************************
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e;
  int start=0;
  apf::Numbering* n = createNumbering(m,"m3dc1_owned_numbering",apf::getConstant(0),1);
  
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (!is_ent_original(m,e)) continue;
    apf::number(n,e,0,0, start++);
  }
  m->end(it);
  std::cout<<"[p"<<PCU_Comm_Self()<<"] "<<__func__<<": #owned_nodes="<<start<<"\n";
  return n;
}



// ***********************************
// 		M3DC1_EPETRA
// ***********************************

Epetra_Map* createEpetraMap(global_ordinal_type num_dof, bool owned)
{
  global_ordinal_type num_node, global_id;
  if (owned) 
    num_node = static_cast<global_ordinal_type>(m3dc1_mesh::instance()->num_own_ent[0]);
  else // overlap
    num_node = static_cast<global_ordinal_type>(m3dc1_mesh::instance()->num_local_ent[0]);

  apf::DynamicArray<global_ordinal_type> dofIndices(num_node*num_dof);

  // loop over owned node
  apf::MeshEntity* e;
  apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
  int i=0;
  while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
  {
    if (owned && !is_ent_original(m3dc1_mesh::instance()->mesh,e)) continue;
    global_id = getNumber(get_global_numbering(), e, 0, 0); 
    for (global_ordinal_type j=0; j < num_dof; ++j)
      dofIndices[i*num_dof + j] = global_id*num_dof + j;
    ++i;
  }
  m3dc1_mesh::instance()->mesh->end(it);
  return new Epetra_Map(-1,dofIndices.getSize(),&dofIndices[0],0,Epetra_MpiComm(MPI_COMM_WORLD));
}

m3dc1_epetra::m3dc1_epetra(int i, int t, int s, FieldID f_id): id(i), matrix_type(t), scalar_type(s), field_id(f_id), num_solver_iter(0)
{
  _field = (*(m3dc1_mesh::instance()->field_container))[f_id]->get_field();
  global_ordinal_type num_dof = static_cast<global_ordinal_type>(apf::countComponents(_field));
  nge = static_cast<global_ordinal_type>(m3dc1_mesh::instance()->num_global_ent[0]*num_dof);
  _owned_map = createEpetraMap(num_dof,true);
  _overlap_map = createEpetraMap(num_dof,false);
  epetra_mat = new Epetra_CrsMatrix(Copy,*_overlap_map,(global_ordinal_type)nge,false/*ignoreNonLocalEntries*/);

/*
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES 
  if (!PCU_Comm_Self()) std::cout<<"nge = "<<nge<<", #my global_eqn = "<<_overlap_map->MyGlobalElements64()<<"(#elem = "<<m3dc1_mesh::instance()->num_local_ent[0]<<"), #global_elm = "<<_overlap_map->NumGlobalElements64()
<<", NumMyRows="<<epetra_mat->NumMyRows()<<"\n";
#else
  std::cout<<"[p"<<PCU_Comm_Self()<<"] nge = "<<nge<<" _overlap_map->Num(Global, My)Elements = "
	   <<_overlap_map->NumMyElements()
           <<", "<<_overlap_map->NumGlobalElements()
           <<", NumGlobal(NonZeros, Rows, Cols)="
	   <<epetra_mat->NumGlobalNonzeros()
           <<", "<<epetra_mat->NumGlobalRows()
	   <<", "<<epetra_mat->NumGlobalCols()
	   <<", NumMy(NonZeros, Rows, Cols)="
	   <<epetra_mat->NumMyNonzeros()
	   <<", "<<epetra_mat->NumMyRows()
	   <<", "<<epetra_mat->NumMyCols()<<"\n";
#endif
*/

#ifdef DEBUG_
  apf::MeshEntity* e;
  apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
  int local_id=0, global_id;
  while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
  {
    m3dc1_node_getglobalid (&local_id, &global_id);
    int start_dof_id=global_id*num_dof;
    for (int i=0; i<num_dof; ++i)
    { 
      if (epetra_mat->LRID(global_id*num_dof+i)!=local_id*num_dof+i)
      {
      std::cout<<"[p"<<PCU_Comm_Self()<<"] epetra_mat->LRID("<<global_id<<"*"<<num_dof<<"+"<<i<<")="
               <<epetra_mat->LRID(global_id*num_dof+i)
               <<", "<<local_id<<"*"<<num_dof<<"+"<<i<<"="<<local_id*num_dof+i
                <<"\n";
      }
      assert(epetra_mat->LRID(global_id*num_dof+i)==local_id*num_dof+i);
    }
    ++local_id;
  }
  m3dc1_mesh::instance()->mesh->end(it);
#endif

#ifdef DEBUG_
  if (!PCU_Comm_Self()) 
    std::cout<<"epetra matrix "<<id<<" created - field: "<<getName(_field)
            <<", num_glob_eq = "<<nge<<"\n";
#endif
}

m3dc1_epetra::~m3dc1_epetra()
{
  delete epetra_mat;
  delete _overlap_map;
  delete _owned_map;
}
#endif
#endif //#ifndef M3DC1_MESHGEN
