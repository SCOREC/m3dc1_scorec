/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_ghost.h"
#include "m3dc1_mesh.h"
#include "m3dc1_scorec.h"
#include "m3dc1_model.h"
#include "m3dc1_field.h"
#include "apf.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apfMDS.h"
#include "apfNumbering.h"
#include "apfShape.h"
#include "PCU.h"
#include <vector>
#include <set>
#include <string.h>
#include <iostream>

using namespace apf;

// m3dc1_ghost
// *********************************************************
m3dc1_ghost::m3dc1_ghost()
// *********************************************************
{
  mesh = NULL;
  nlayers = 0;
  //is_ghosted = false;
  field_container=NULL;
  reset();
  ordering_opt=M3DC1_NO_ORDER;
  local_entid_tag=own_partid_tag=num_global_adj_node_tag=num_own_adj_node_tag=NULL;
}

// *********************************************************
m3dc1_ghost::~m3dc1_ghost()
// *********************************************************
{
  if (!mesh)
    return;
  std::set<int> fields_keep; 
  clean(fields_keep);
  mesh->destroyNative();
  apf::destroyMesh(mesh);
}

void m3dc1_ghost:: clean(std::set<int>& fields_keep)
{
 // destroy field AND numbering
  if (field_container)
  {
    for (std::map<FieldID, m3dc1_field*>::iterator f_it=field_container->begin(); f_it!=field_container->end();)
    {
      if(fields_keep.find(f_it->first)!=fields_keep.end())
      {
        f_it++;
        continue;
      }
      // if(!PCU_Comm_Self()) std::cout<<" destroy field "<<getName(f_it->second->get_field())<<std::endl;
      FieldID id = f_it->first;
      std::map<FieldID, m3dc1_field*>::iterator it_next=++f_it;
      m3dc1_field_delete(&id);
      f_it=it_next;
    }
    //field_container->clear();
  }
  if(!fields_keep.size()) {delete field_container; field_container=0;}

  // destroy tag data
  for (int d=0; d<4; ++d)
  {
    removeTagFromDimension(mesh, local_entid_tag, d);
    removeTagFromDimension(mesh, own_partid_tag, d);
  }
  removeTagFromDimension(mesh, num_global_adj_node_tag, 0);
  removeTagFromDimension(mesh, num_own_adj_node_tag, 0);

  // destroy mesh
  //is_ghosted = false;
  global_ghost_state = uninit;
  mesh->destroyTag(local_entid_tag);
  mesh->destroyTag(own_partid_tag);
  mesh->destroyTag(num_global_adj_node_tag);
  mesh->destroyTag(num_own_adj_node_tag);
}

m3dc1_ghost* m3dc1_ghost::_instance=NULL;

m3dc1_ghost* m3dc1_ghost::instance()
{
  if (_instance==NULL)
    _instance = new m3dc1_ghost();
  return _instance;
}

void m3dc1_ghost::destroy()
{
  delete _instance;
  _instance = NULL;
}

// *********************************************************
void m3dc1_ghost::reset()
// *********************************************************
{
  for (int i=0; i<4; ++i)
  {
    num_local_ent[i] = 0;
    num_own_ent[i] = 0;
    num_global_ent[i] = 0;
  }
}

// **********************************************
bool is_ent_original_ghost(Mesh2* mesh, MeshEntity* e)
// **********************************************
{
  if (!mesh->hasMatching())
    return (PCU_Comm_Self() == get_ent_ownpartid_ghost(mesh, e));

  Matches ms;
  mesh->getMatches(e,ms);
  int self = PCU_Comm_Self();
  for (size_t i = 0; i < ms.getSize(); ++i)
  {
    if (ms[i].peer < self)
      return false;
    if ((ms[i].peer == self) &&
        (ms[i].entity < e))
      return false;
  }
  return true;
}


// **********************************************
int get_ent_ownpartid_ghost(Mesh2* mesh, MeshEntity* e)
// **********************************************
{
  int own_partid;

  if (mesh->hasTag(e, m3dc1_ghost::instance()->own_partid_tag))
    mesh->getIntTag(e, m3dc1_ghost::instance()->own_partid_tag, &own_partid);
  else
  {
    std::cout<<"[M3D-C1 WARNING] (p"<<PCU_Comm_Self()<<") not found own_partid_tag for e (dim"<<getDimension(mesh, e)
             <<", id "<<getMdsIndex(mesh, e)<<") \n";
    own_partid=mesh->getOwner(e);
  }
  return own_partid;
}

// **********************************************
MeshEntity* get_ent_owncopy_ghost(Mesh2* mesh, MeshEntity* e)
// **********************************************
{
  if (!(mesh->isShared(e))) // internal ent
    return e;

  int own_partid = get_ent_ownpartid_ghost(mesh, e);
  if (own_partid==PCU_Comm_Self()) return e;

  Copies remotes;
  mesh->getRemotes(e,remotes);
  return remotes[own_partid];
}

 struct entMsg_ghost
  {
    int pid;
    MeshEntity* ent;
    entMsg_ghost( int pid_p=0, MeshEntity* ent_p=NULL)
    {
      pid=pid_p;
      ent=ent_p;
    }
  };
  struct classcomp_ghost
  {
    bool operator() (const entMsg_ghost& lhs, const entMsg_ghost& rhs) const
    {
      if(lhs.ent==rhs.ent) return lhs.pid<rhs.pid;
      else return lhs.ent<rhs.ent;
    }
  };


// **********************************************
void m3dc1_ghost::set_node_adj_tag()
// **********************************************
{
  int value;
  int brgType = (num_local_ent[3])?3:2;

  apf::MeshEntity* e;
  apf::MeshIterator* it = mesh->begin(0);
  PCU_Comm_Begin();
  while ((e = mesh->iterate(it)))
  {
    int num_adj_node=0;
    Adjacent elements;
    getBridgeAdjacent(mesh, e, brgType, 0, elements);
    int num_adj = elements.getSize();

    for (int i=0; i<num_adj; i++)
    {
      if (is_ent_original_ghost(mesh, elements[i]))
        ++num_adj_node;
    }
    mesh->setIntTag(e, num_own_adj_node_tag, &num_adj_node);

    if (!mesh->isShared(e)) continue;
    // first pass msg size to owner
    int own_partid = get_ent_ownpartid_ghost(mesh, e);
    MeshEntity* own_copy = get_ent_owncopy_ghost(mesh, e);

    if (own_partid==PCU_Comm_Self()) continue;
    PCU_COMM_PACK(own_partid, own_copy);
    PCU_Comm_Pack(own_partid, &num_adj,sizeof(int));
  }
  mesh->end(it);

  PCU_Comm_Send();

  std::map<apf::MeshEntity*, std::map<int, int> > count_map;
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(e);
      PCU_Comm_Unpack(&value,sizeof(int));
      count_map[e][PCU_Comm_Sender()]=value;
    }
  }
  
  // pass entities to ownner
  std::map<apf::MeshEntity*, std::set<entMsg_ghost, classcomp_ghost> > count_map2;
  it = mesh->begin(0);
  PCU_Comm_Begin();
  while ((e = mesh->iterate(it)))
  {
    // pass entities to ownner

    std::vector<entMsg_ghost> msgs;
    Adjacent elements;
    getBridgeAdjacent(mesh, e, brgType, 0, elements);

    MeshEntity* ownerEnt=get_ent_owncopy_ghost(mesh, e);
    int own_partid = get_ent_ownpartid_ghost(mesh, e);
    for(int i=0; i<elements.getSize(); i++)
    {
      MeshEntity* ownerEnt2=get_ent_owncopy_ghost(mesh, elements[i]);
      int owner=get_ent_ownpartid_ghost(mesh, elements[i]);
      msgs.push_back(entMsg_ghost(owner, ownerEnt2));
      if(own_partid==PCU_Comm_Self()) 
      {
        count_map2[e].insert(*msgs.rbegin());
      }
    }

    if(own_partid!=PCU_Comm_Self())
    {
      PCU_COMM_PACK(own_partid, ownerEnt);
      PCU_Comm_Pack(own_partid, &msgs.at(0),sizeof(entMsg_ghost)*msgs.size());
    }
  }
  mesh->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(e);
      int sizeData = count_map[e][PCU_Comm_Sender()];
      std::vector<entMsg_ghost> data(sizeData);
      PCU_Comm_Unpack(&data.at(0),sizeof(entMsg_ghost)*sizeData);
      for (int i=0; i<data.size(); i++)
      {
        count_map2[e].insert(data.at(i));
      }
    }
  }

  for (std::map<apf::MeshEntity*, std::set<entMsg_ghost,classcomp_ghost> >::iterator mit=count_map2.begin(); mit!=count_map2.end(); ++mit)
  {
    e = mit->first;
    int num_global_adj =count_map2[e].size();
    mesh->setIntTag(mit->first, num_global_adj_node_tag, &num_global_adj);
  }
}


// *********************************************************
void m3dc1_ghost::initialize()
// *********************************************************
{
  if(!local_entid_tag) local_entid_tag = mesh->createIntTag("m3dc1_local_ent_id", 1);
  if(!own_partid_tag) own_partid_tag = mesh->createIntTag("m3dc1_own_part_id", 1);
  if(!num_global_adj_node_tag) num_global_adj_node_tag = mesh->createIntTag("m3dc1_num_global_adj_node", 1);
  if(!num_own_adj_node_tag) num_own_adj_node_tag = mesh->createIntTag("m3dc1_num_own_adj_node", 1);

  reset();
  MeshEntity* e;

  int counter = 0, own_partid, local_partid=PCU_Comm_Self();

  for (int d=0; d<4; ++d)
  {
    num_local_ent[d] = countEntitiesOfType(mesh, d);
    counter=0;
    MeshIterator* it = mesh->begin(d);
    while ((e = mesh->iterate(it)))
    {
      mesh->setIntTag(e, local_entid_tag, &counter);
      own_partid = mesh->getOwner(e);
      mesh->setIntTag(e, own_partid_tag, &own_partid);
      if (own_partid==local_partid)
        ++num_own_ent[d];
      ++counter;
    }
    mesh->end(it);
  }

  MPI_Allreduce(num_own_ent, num_global_ent, 4, MPI_INT, MPI_SUM, PCU_Get_Comm());
  set_node_adj_tag();

  // Copy the fields on the ghosted mesh into the field container
  if (!m3dc1_ghost::instance()->field_container)
    m3dc1_ghost::instance()->field_container = new std::map<FieldID, m3dc1_field*>;
  for (std::map<FieldID, m3dc1_field*>::iterator it =
	 m3dc1_mesh::instance()->field_container->begin();
       it !=  m3dc1_mesh::instance()->field_container->end();
       ++it) {
    int field_id = it->first;
    int num_values = it->second->get_num_value();
    int scalar_type = it->second->get_value_type();
    int num_dofs_per_value = it->second->get_dof_per_value();
    apf::Field* old_field = it->second->get_field();
    apf::Field* new_field = mesh->findField(apf::getName(old_field));
    
    m3dc1_ghost::instance()->field_container->insert(
        std::map<FieldID, m3dc1_field*>::value_type(field_id,
          new m3dc1_field(field_id,
                          new_field,
                          num_values,
                          scalar_type,
                          num_dofs_per_value)));
    apf::freeze(new_field);
  }

  // Set global tracking tag
  global_ghost_state = init;
}

