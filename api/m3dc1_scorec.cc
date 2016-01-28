/****************************************************************************** 

  (c) 2005-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "m3dc1_scorec.h"
#include "m3dc1_matrix.h"
#include "m3dc1_model.h"
#include "m3dc1_mesh.h"
#include "m3dc1_ghost.h"
#include "m3dc1_field.h"
#include <mpi.h>
#include <PCU.h>
#include <gmi_analytic.h>
#include <map>
#include "apfMDS.h"
#include "apfNumbering.h"
#include "Expression.h"
#include "m3dc1_slnTransfer.h"
#include "m3dc1_sizeField.h"
#include "ReducedQuinticImplicit.h"
#include "apfMesh2.h"
#include "apfOmega_h.h"

bool m3dc1_double_isequal(double A, double B)
{
  double maxDiff = 1e-5;
  double maxRelDiff = 1e-5;
// http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/ 
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    double diff = fabs(A - B);
    if (diff <= maxDiff)
        return true;
 
    A = fabs(A);
    B = fabs(B);
    double largest = (B > A) ? B : A;
 
    if (diff <= largest * maxRelDiff)
        return true;
    return false;
}

//*******************************************************
int m3dc1_scorec_init()
//*******************************************************
{ 
  PCU_Comm_Init();

  m3dc1_model :: instance()-> model = gmi_make_analytic();
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_scorec_finalize()
//*******************************************************
{ 
  // destroy global numbering and associated field
  // delete existing numbering
  apf::Numbering* old_n = m3dc1_mesh::instance()->mesh->findNumbering("m3dc1_global_node_id");
  if (old_n) destroyNumbering(old_n);
  int node_glb_order=NODE_GLB_ORDER; 
  m3dc1_field_delete (&node_glb_order);
  m3dc1_gfield_delete (&node_glb_order);

  //actually destroy badly designed singletons.
  //this was chosen over statically allocating the objects
  //and having the runtime deallocate them in order to avoid
  //possible issues linking to FORTRAN.
  //feel free to make them static objects and see if that works
  m3dc1_ghost::destroy();
  m3dc1_mesh::destroy();
  m3dc1_model::destroy();

  PCU_Comm_Free();
  return M3DC1_SUCCESS; 
}

/** plane functions */

//*******************************************************
int m3dc1_plane_setnum(int* num_plane)
//*******************************************************
{
  if (*num_plane<1 || PCU_Comm_Peers()%(*num_plane)) return M3DC1_FAILURE;
  m3dc1_model::instance()->set_num_plane(*num_plane);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_plane_getnum(int* num_plane)
//*******************************************************
{
  *num_plane = m3dc1_model::instance()->num_plane;
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_plane_getid(int * plane_id)
//*******************************************************
{
  *plane_id = m3dc1_model::instance()->local_planeid;
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_plane_setphirange(double* min_val, double* max_val)
//*******************************************************
{
  m3dc1_model::instance()->set_phi(*min_val, *max_val);
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_plane_setphi(int* planeid, double* phi)
//*******************************************************
{
  m3dc1_model::instance()->set_phi(*planeid, *phi);
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_plane_getphi(int* planeid, double* phi)
//*******************************************************
{
  m3dc1_model::instance()->get_phi(*planeid, phi);
  return M3DC1_SUCCESS; 
}

/** model functions */
//*******************************************************
int m3dc1_model_getplaneid(int * plane_id)
//*******************************************************
{
  *plane_id = m3dc1_model::instance()->local_planeid;
}

//*******************************************************
int m3dc1_model_getedge (int*  /* out */  left_edge, 
                         int*  /* out */  right_edge,
                         int*  /* out */  bottom_edge, 
                         int*  /* out */  top_edge)
//*******************************************************
{
  std::cout<<"m3dc1_model_getedge not implemented"; throw 1;
}

//*******************************************************
int m3dc1_model_getmincoord(double* x_min, double* y_min)
//*******************************************************
{
  *x_min = m3dc1_model::instance()->boundingBox[0];
  *y_min = m3dc1_model::instance()->boundingBox[1];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_getmaxcoord(double* x_max, double* y_max)
//*******************************************************
{
  double mincoord[3],maxcoord[3];
  *x_max = m3dc1_model::instance()->boundingBox[2];
  *y_max = m3dc1_model::instance()->boundingBox[3];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_load(const char* /* in */ model_file)
//*******************************************************
{  
  std::string str_model_name(model_file);
  m3dc1_model::instance()->load_analytic_model(model_file); 
  m3dc1_model::instance()->caculateBoundingBox();
  // save the num of geo ent on the oringal plane
  m3dc1_model::instance()->numEntOrig[0]=m3dc1_model::instance()->model->n[0];
  m3dc1_model::instance()->numEntOrig[1]=m3dc1_model::instance()->model->n[1];
  m3dc1_model::instance()->numEntOrig[2]=m3dc1_model::instance()->model->n[2];

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_print()
//*******************************************************
{
  if (PCU_Comm_Self() || m3dc1_model::instance()->local_planeid) 
    return M3DC1_SUCCESS;

  double min[3], max[3];
  gmi_iter* gf_it = gmi_begin(m3dc1_model::instance()->model, 2);
  gmi_ent* ge;
  int count = 0;
  while ((ge = gmi_next(m3dc1_model::instance()->model, gf_it))) 
  {
    gmi_set* gf_edges = gmi_adjacent(m3dc1_model::instance()->model, ge, 1);
    if (!PCU_Comm_Self())    std::cout<<"model face id "<<gmi_tag(m3dc1_model::instance()->model, ge)
             <<" - # adj model edges: "<<gf_edges->n<<"\n";
    if (gf_edges->n)
    {
    if (!PCU_Comm_Self())      std::cout<<"\tadj edge ID: ";
      for (int i=0; i<gf_edges->n; ++i)
    if (!PCU_Comm_Self())        std::cout<<gmi_tag(m3dc1_model::instance()->model,  gf_edges->e[i])<<" "; 
    if (!PCU_Comm_Self())      std::cout<<"\n";
    }
    gmi_free_set(gf_edges);
  }
  gmi_end(m3dc1_model::instance()->model, gf_it);

// verify geom edge
  gmi_iter* ge_it = gmi_begin(m3dc1_model::instance()->model, 1);
  while ((ge = gmi_next(m3dc1_model::instance()->model, ge_it))) 
  {
    M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,ge);
    if (!pn)   
    {
    if (!PCU_Comm_Self())       std::cout<<"["<<PCU_Comm_Self()<<"] model edge "<<gmi_tag(m3dc1_model::instance()->model, ge)
             <<" - failed with gmi_analytic_data retrieval\n";
    }
  }
  gmi_end(m3dc1_model::instance()->model, ge_it);

// verify gv_edges
  gmi_iter* gv_it = gmi_begin(m3dc1_model::instance()->model, 0);
  while ((ge = gmi_next(m3dc1_model::instance()->model, gv_it))) 
  {
    gmi_set* gv_edges = gmi_adjacent(m3dc1_model::instance()->model, ge, 1);
    if (!PCU_Comm_Self())
      std::cout<<"model vertex id "<<gmi_tag(m3dc1_model::instance()->model, ge)
             <<" - # adj model edges: "<<gv_edges->n<<"\n";
    if (gv_edges->n)
    {
    if (!PCU_Comm_Self())
      std::cout<<"\tadj edge ID: ";
      for (int i=0; i<gv_edges->n; ++i)
            if (!PCU_Comm_Self()) std::cout<<gmi_tag(m3dc1_model::instance()->model,  gv_edges->e[i])<<" "; 
          if (!PCU_Comm_Self()) std::cout<<"\n";
    }
    gmi_free_set(gv_edges);
  }
  gmi_end(m3dc1_model::instance()->model, gv_it);

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_setnumplane(int* num_plane)
//*******************************************************
{
  if (*num_plane<1 || PCU_Comm_Peers()%(*num_plane)) return M3DC1_FAILURE;
  m3dc1_model::instance()->set_num_plane(*num_plane);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_getnumplane(int* num_plane)
//*******************************************************
{
  *num_plane = m3dc1_model::instance()->num_plane;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_setpbc (int* /* in */ x_pbc, 
                        int* /* in */ y_pbc)
//*******************************************************
{
  m3dc1_model::instance()->xperiodic = *x_pbc;
  m3dc1_model::instance()->yperiodic = *y_pbc;
  return  M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_model_getpbc (int* /* out */ x_pbc, 
                        int* /* out */ y_pbc)
//*******************************************************
{
  *x_pbc = m3dc1_model::instance()->xperiodic;
  *y_pbc = m3dc1_model::instance()->yperiodic;
  return  M3DC1_SUCCESS;
}

/** mesh functions */
#include <parma.h>

void setWeight(apf::Mesh* m, apf::MeshTag* tag, int dim) {
  double w = 1.0;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
    m->setDoubleTag(e, tag, &w);
  m->end(it);
}

apf::MeshTag* setWeights(apf::Mesh* m) {
  apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
  setWeight(m, tag, 0);
  setWeight(m, tag, m->getDimension());
  return tag;
}

void clearTags(apf::Mesh* m, apf::MeshTag* t) {
  apf::removeTagFromDimension(m, t, 0);
  apf::removeTagFromDimension(m, t, m->getDimension());
}

//*******************************************************
int m3dc1_mesh_load(const char* mesh_file)
//*******************************************************
{ 
  if (m3dc1_model::instance()->local_planeid == 0) 
  {
    m3dc1_mesh::instance()->mesh = apf::loadMdsMesh(m3dc1_model::instance()->model, mesh_file);
    /* vertex load balancing */
    Parma_PrintPtnStats(m3dc1_mesh::instance()->mesh, "initial");
    // will not work not non-man geo 
    //m3dc1_mesh::instance()->mesh->verify();
    apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
    while(mesh->countFields())
    {
      apf::Field* f = mesh->getField(0);
      //if (!PCU_Comm_Self()) std::cout<<"[M3D-C1 INFO] (p"<<PCU_Comm_Self()<<") destroying field "<<getName(f)<<endl;
      destroyField(f);
    }
    while(mesh->countNumberings())
    {
      apf::Numbering* n = mesh->getNumbering(0);
      destroyNumbering(n);
    }
    apf::DynamicArray<apf::MeshTag*> tags;
    mesh->getTags(tags);
    for(int i=0; i<tags.getSize(); i++)
    {
      if (mesh->findTag("norm_curv")==tags[i]) continue;
      for(int idim=0; idim<4; idim++)
        apf::removeTagFromDimension(mesh, tags[i], idim);
      mesh->destroyTag(tags[i]);
    }
  }
  else {
    m3dc1_mesh::instance()->mesh = apf::makeEmptyMdsMesh(m3dc1_model::instance()->model, 2, false);
  }
  m3dc1_mesh::instance()->initialize();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_mesh_build3d (int* num_field, int* field_id,  
                        int* num_dofs_per_value)
//*******************************************************
{ 
  // switch COMM to GLOBAL COMM
  MPI_Comm groupComm = PCU_Get_Comm();
  PCU_Switch_Comm(m3dc1_model::instance()->oldComm);
  MPI_Comm_free(&groupComm);

  // initialize phi value and construct 3d
  int num_plane = m3dc1_model::instance()->num_plane;
  m3dc1_model::instance()->set_phi(0.0, 2.0*M3DC1_PI/num_plane*(num_plane-1));
  m3dc1_model::instance()->setupCommGroupsPlane();
  m3dc1_mesh::instance()->build3d(*num_field, field_id, num_dofs_per_value);
  // now construct 3d mesh
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_getnument (int* /* in*/ ent_dim, int* /* out */ num_ent)
//*******************************************************
{
  if (*ent_dim<0 || *ent_dim > 3)
    return M3DC1_FAILURE;
  *num_ent = m3dc1_mesh::instance()->num_local_ent[*ent_dim];
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_getnumownent (int* /* in*/ ent_dim, int* /* out */ num_ent)
//*******************************************************
{
  if (*ent_dim<0 || *ent_dim > 3)
    return M3DC1_FAILURE;
  *num_ent = m3dc1_mesh::instance()->num_own_ent[*ent_dim];
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_getnumglobalent (int* /* in*/ ent_dim, int* /* out */ num_ent)
//*******************************************************
{
  if (*ent_dim<0 || *ent_dim > 3)
    return M3DC1_FAILURE;
  *num_ent = m3dc1_mesh::instance()->num_global_ent[*ent_dim];
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_setordering (int* option)
//*******************************************************
{
  if (*option < M3DC1_NO_ORDER || *option > M3DC1_ADJ_SOLVER_ORDER)
    return M3DC1_FAILURE;
  m3dc1_mesh::instance()->ordering_opt = *option;
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_mesh_getordering (int* option)
//*******************************************************
{
  *option = m3dc1_mesh::instance()->ordering_opt;
  return M3DC1_SUCCESS;
}

/* mesh entity functions */
//*******************************************************
int m3dc1_node_getglobalid (int* /* in */ ent_id, int* /* out */ global_ent_id)
//*******************************************************
{
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, *ent_id);
  if (!e)
    return M3DC1_FAILURE;

  apf::Numbering* n = get_global_numbering();
  *global_ent_id = getNumber(n, e, 0, 0);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getgeomclass (int* /* in */ ent_dim, int* /* in */ ent_id, 
            int* /* out */ geom_class_dim, int* /* out */ geom_class_id)
//*******************************************************
{ 
  apf::MeshEntity* ent = getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  assert(ent);
  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->mesh->toModel(ent));
  *geom_class_dim = gmi_dim(m3dc1_model::instance()->model,gent);
  *geom_class_id = gmi_tag(m3dc1_model::instance()->model,gent);
  // if 3D mesh, need to return the classification on the original plane
  if( m3dc1_mesh::instance()->mesh->getDimension() ==3 )
  {
    int numEntOrig[3];
    int numPlane = m3dc1_model::instance()->num_plane;
    memcpy( numEntOrig, m3dc1_model::instance()->numEntOrig, sizeof(numEntOrig));
    *geom_class_id-=1;
    switch (*geom_class_dim)
    {
      case 3: *geom_class_id%=numEntOrig[2]; break;
      case 2: *geom_class_id%=(numEntOrig[1]+numEntOrig[2]); break;
      case 1:  *geom_class_id%=(numEntOrig[0]+numEntOrig[1]); break;
      case 0: *geom_class_id%=(numEntOrig[0]);
    }
    *geom_class_id+=1;
  }
}

//*******************************************************
int m3dc1_ent_getadj (int* /* in */ ent_dim, int* /* in */ ent_id, 
                      int* /* in */ adj_dim, int* /* out */ adj_ent, 
                      int* /* in */ adj_ent_allocated_size, int* /* out */ num_adj_ent)
//*******************************************************
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  if (!e || *adj_dim==*ent_dim)
    return M3DC1_FAILURE;

  if (*adj_dim>*ent_dim) // upward
  {
    apf::Adjacent adjacent;
    m3dc1_mesh::instance()->mesh->getAdjacent(e,*adj_dim,adjacent);
    *num_adj_ent = adjacent.getSize();
    if (*adj_ent_allocated_size<*num_adj_ent)
      return M3DC1_FAILURE;
    for (int i=0; i<*num_adj_ent; ++i)
      adj_ent[i] = getMdsIndex(m3dc1_mesh::instance()->mesh, adjacent[i]);
  }
  else if (*adj_dim<*ent_dim) 
  {
    apf::Downward downward;
    *num_adj_ent = m3dc1_mesh::instance()->mesh->getDownward(e, *adj_dim, downward);
    if (*adj_ent_allocated_size<*num_adj_ent)
      return M3DC1_FAILURE;
    for (int i=0; i<*num_adj_ent; ++i)
      adj_ent[i] = getMdsIndex(m3dc1_mesh::instance()->mesh, downward[i]);
    //adjust the order to work with m3dc1
    if (m3dc1_mesh::instance()->mesh->getDimension()==3 && *ent_dim==3 &&*adj_dim==0 &&adj_ent[0]>adj_ent[3])
    {
      int buff[3];
      memcpy(buff, adj_ent, 3*sizeof(int));
      memcpy(adj_ent, adj_ent+3, 3*sizeof(int));
      memcpy(adj_ent+3, buff, 3*sizeof(int));
    }
  }
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_ent_getnumadj (int* /* in */ ent_dim, int* /* in */ ent_id, 
                         int* /* in */ adj_dim, int* /* out */ num_adj_ent)
//*******************************************************
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  if (!e || *adj_dim==*ent_dim)
    return M3DC1_FAILURE;

  if (*adj_dim>*ent_dim) // upward
  {
    apf::Adjacent adjacent;
    m3dc1_mesh::instance()->mesh->getAdjacent(e,*adj_dim,adjacent);
    *num_adj_ent = adjacent.getSize();
  }
  else if (*adj_dim<*ent_dim) 
  {
    apf::Downward downward;
    *num_adj_ent = m3dc1_mesh::instance()->mesh->getDownward(e, *adj_dim, downward);
  }
  return M3DC1_SUCCESS; 
}

//*******************************************************
int m3dc1_ent_getownpartid (int* /* in */ ent_dim, int* /* in */ ent_id, 
                            int* /* out */ owning_partid)
//*******************************************************
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  if (!e)
    return M3DC1_FAILURE;
  *owning_partid = get_ent_ownpartid(m3dc1_mesh::instance()->mesh, e);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_ismine (int* /* in */ ent_dim, int* /* in */ ent_id, 
                            int* /* out */ ismine)
//*******************************************************
{
  *ent_id -= 1; //index change from Fortran to C
 
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e = getMdsEntity(m, *ent_dim, *ent_id);
  if (!e)
    return M3DC1_FAILURE;

  if (is_ent_original(m,e)) 
     *ismine = 1;   // 
  else
     *ismine = 0; 
  return M3DC1_SUCCESS;
}

// node-specific functions
//*******************************************************
int m3dc1_node_getcoord (int* /* in */ node_id, double* /* out */ coord)
//*******************************************************
{
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, *node_id);
  if (!e)
    return M3DC1_FAILURE;
  apf::Vector3 xyz;
  m3dc1_mesh::instance()->mesh->getPoint(e, 0, xyz);
  for (int i=0; i<3; ++i)
    coord[i] = xyz[i]; 
  return M3DC1_SUCCESS;
}

//*******************************************************
void get_gv_edges(gmi_ent* gvertex, std::vector<gmi_ent*>& gedges)
//*******************************************************
{
  gedges.clear();
  gmi_set* gv_edges = gmi_adjacent(m3dc1_model::instance()->model, gvertex, 1);
  for (int i=0; i<gv_edges->n; ++i)
    gedges.push_back(gv_edges->e[i]);
  gmi_free_set(gv_edges);
  assert(gedges.size()>=1);
}

//*******************************************************
int m3dc1_node_getnormvec (int* /* in */ node_id, double* /* out */ xyzt)
//*******************************************************
{
  apf::MeshEntity* vt = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, *node_id);
  assert(vt);
  xyzt[2]=0.0;
  //cout<<"nodnormalvec_ "<<*iNode<<" "<<vt<<endl;
  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->mesh->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  if(gType !=  1 && gType !=  0)
  {
    xyzt[0] = xyzt[1] = 0.0;
    return M3DC1_SUCCESS;
  }

  apf::MeshTag* norm_curv_tag = m3dc1_mesh::instance()->mesh->findTag("norm_curv");
  if (norm_curv_tag && m3dc1_mesh::instance()->mesh->hasTag(vt, norm_curv_tag))
  {
    double norm_curv[3];
    m3dc1_mesh::instance()->mesh->getDoubleTag(vt, norm_curv_tag, &norm_curv[0]);
    xyzt[0] = norm_curv[0]; 
    xyzt[1] = norm_curv[1];
    return M3DC1_SUCCESS;
  }
  else
  { // if norm/curv is not attached, evaluate
    apf::Vector3 param(0,0,0);
    m3dc1_mesh::instance()->mesh->getParam(vt,param);
    // geo node avage on the connected edges
    if (gType == 0) // node is on the 
    {
      apf::Vector3 vcd_t;
      double vcd[3];
      m3dc1_mesh::instance()->mesh->getPoint(vt, 0, vcd_t);
      for(int i=0; i<3; i++) 
        vcd[i]=vcd_t[i];
      std::vector<gmi_ent*> gEdges;
      get_gv_edges(gent, gEdges);
      int numEdgePlane=0;
      double normalvec[3]={0.,0.,0.};
      xyzt[0]=xyzt[1]=xyzt[2]=0;
      if (gEdges.size()<2) 
        std::cout<<"["<<PCU_Comm_Self()<<"] "<<__func__<<" ERROR: #adjEdge of gVertex="<<gEdges.size()<<" (it should be minimum 2) \n";
      assert(gEdges.size()>=2);
      for(int i=0;i<gEdges.size();i++)
      {
        gmi_ent* pe = gEdges.at(i);
        double cd[3]={0,0,0};
        double paraRange[2];
        gmi_range(m3dc1_model::instance()->model, pe, 0, paraRange);
        M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,pe);
        if(!pn) continue;
        numEdgePlane++;
        M3DC1::evalCoord(paraRange[0],cd, pn);
        if(checkSamePoint2D(vcd,cd))
        {
          M3DC1::evalNormalVector(pn[0],pn[1], paraRange[0], normalvec);
        }
        else
        {
          evalNormalVector(pn[0],pn[1], paraRange[1], normalvec);
        }
        xyzt[0]+=normalvec[0];
        xyzt[1]+=normalvec[1];
        xyzt[2]+=normalvec[2];
      }
      if (numEdgePlane!=2) 
        std::cout<<"["<<PCU_Comm_Self()<<"] "<<__func__<<" ERROR: numEdgePlane="<<numEdgePlane<<" (it should be 2) \n";
      assert(numEdgePlane==2);
      double arclen=sqrt(xyzt[0]*xyzt[0]+xyzt[1]*xyzt[1]+xyzt[2]*xyzt[2]);
      assert(arclen>0);
      xyzt[0]=xyzt[0]/arclen;
      xyzt[1]=xyzt[1]/arclen;
      xyzt[2]=xyzt[2]/arclen;
    }
    else
    {
      apf::Vector3 param(0,0,0);
      m3dc1_mesh::instance()->mesh->getParam(vt,param);
      M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,gent);
      evalNormalVector(pn[0],pn[1], param[0], xyzt);
    }
  }
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_node_getcurv (int* /* in */ node_id, double* /* out */ curv)
//*******************************************************
{
  apf::MeshEntity* vt = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, *node_id);
  assert(vt);

  apf::MeshTag* norm_curv_tag = m3dc1_mesh::instance()->mesh->findTag("norm_curv");
  if (norm_curv_tag && m3dc1_mesh::instance()->mesh->hasTag(vt, norm_curv_tag))
  {
    double norm_curv[3];
    m3dc1_mesh::instance()->mesh->getDoubleTag(vt, norm_curv_tag, &norm_curv[0]);
    *curv = norm_curv[2]; 
    return M3DC1_SUCCESS;
  }

  *curv=0.0;
  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->mesh->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  if(gType==0)
  {
    apf::Vector3 vcd_t;
    double vcd[3];
    m3dc1_mesh::instance()->mesh->getPoint(vt, 0, vcd_t);
    for(int i=0; i<3; i++)
      vcd[i]=vcd_t[i];
    std::vector<gmi_ent*> gEdges;
    get_gv_edges(gent, gEdges);
    int numEdgesPlane=0;
    double curv_tmp;
    for(int i=0;i<gEdges.size();i++)
    {
      gmi_ent* pe = gEdges.at(i);
      double cd[3]={0,0,0};
      double paraRange[2];
      gmi_range(m3dc1_model::instance()->model, pe, 0, paraRange);
      M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,pe);
      if(!pn) continue;
      numEdgesPlane++;
      evalCoord(paraRange[0], cd, pn);
      if(checkSamePoint2D(vcd,cd))
      {
        evalCurvature(pn[0],pn[1], paraRange[0], &curv_tmp); 
      }
      else
      {
        evalCurvature(pn[0],pn[1], paraRange[1], &curv_tmp);
      }
      *curv+=curv_tmp;
    }
    assert(numEdgesPlane==2);
    *curv/=numEdgesPlane;
  }
  else if(gType==1)
  {
      apf::Vector3 param(0,0,0);
      m3dc1_mesh::instance()->mesh->getParam(vt,param);
      M3DC1::Expression** pn=(M3DC1::Expression**) gmi_analytic_data(m3dc1_model::instance()->model,gent);
      evalCurvature(pn[0],pn[1], param[0], curv);
  }
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_node_isongeombdry (int* /* in */ node_id, int* /* out */ on_geom_bdry)
//*******************************************************
{
  apf::MeshEntity* vt = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, *node_id);
  assert(vt);
  gmi_ent* gent= (gmi_ent*)(m3dc1_mesh::instance()->mesh->toModel(vt));
  int gType = gmi_dim(m3dc1_model::instance()->model,gent);
  *on_geom_bdry=(gType==0||gType==1);
  return M3DC1_SUCCESS;
}

// region-specific function
//*******************************************************
int m3dc1_region_getoriginalface( int * /* in */ elm, int * /* out */ fac)
//*******************************************************
{
  apf::MeshEntity* ent = getMdsEntity(m3dc1_mesh::instance()->mesh, 3, *elm);
  apf::Downward downward;
  int num_adj_ent = m3dc1_mesh::instance()->mesh->getDownward(ent, 2, downward);
  assert(num_adj_ent==5);
  int triFace[2];
  int counter=0;
  for(int i=0; i<num_adj_ent; i++)
  {
    int num_adj_ent;
    apf::Downward downward2;
    int num_edge= m3dc1_mesh::instance()->mesh->getDownward(downward[i], 1, downward2);
    if(num_edge==3) triFace[counter++]= getMdsIndex(m3dc1_mesh::instance()->mesh,downward[i]);    
  }
  assert(counter==2);
  *fac = std::min(triFace[0],triFace[1]);
  return M3DC1_SUCCESS;
}

/** field manangement */
int fieldIdMax=0;
//*******************************************************
int m3dc1_field_getnewid ( FieldID* /*out*/field_id )
//*******************************************************
{
  *field_id = fieldIdMax+1;
  return M3DC1_SUCCESS;
}

// *scalar_type is either M3DC1_REAL or M3DC1_COMPLEX
int m3dc1_field_create (FieldID* /*in*/ field_id, const char* /* in */ field_name, int* /*in*/ num_values, 
int* /*in*/ scalar_type, int* /*in*/ num_dofs_per_value)
{
  if (!m3dc1_mesh::instance()->field_container)
    m3dc1_mesh::instance()->field_container=new std::map<FieldID, m3dc1_field*>;

  // shape evaluation will be performed outside the APF
  // only need to tell APF all dofs are attached to mesh vertex
  int components = (*num_values)*(*scalar_type+1)*(*num_dofs_per_value);
  apf::Field* f = createPackedField(m3dc1_mesh::instance()->mesh, field_name, components);
  m3dc1_mesh::instance()->field_container->insert(std::map<FieldID, m3dc1_field*>::value_type(*field_id, new m3dc1_field(*field_id, f, *num_values, *scalar_type, *num_dofs_per_value)));
  apf::freeze(f); // switch dof data from tag to array

#ifdef DEBUG
  if (!PCU_Comm_Self()) 
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": ID "<<*field_id<<", #values "<<*num_values<<", #dofs "<<countComponents(f)<<", name "<<field_name<<"\n";
#endif

  if (*field_id>fieldIdMax) fieldIdMax=*field_id;
  double val[2]={0,0};
  m3dc1_field_assign(field_id, val, scalar_type);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_delete (FieldID* /*in*/ field_id)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container)
    return M3DC1_FAILURE;
  if (!m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;

  apf::Field* f = (*m3dc1_mesh::instance()->field_container)[*field_id]->get_field();
#ifdef DEBUG
  if (!PCU_Comm_Self()) 
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": ID "<<*field_id<<", name "<<getName(f)<<"\n";
#endif

  destroyField(f);

  // remove f from field container
  delete (*m3dc1_mesh::instance()->field_container)[*field_id];
  m3dc1_mesh::instance()->field_container->erase(*field_id);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getinfo(FieldID* /*in*/ field_id, 
                        char* /* out*/ field_name, int* num_values, 
                        int* scalar_type, int* total_num_dof)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf ->get_field();
  strcpy(field_name, getName(f));
  *num_values = mf -> get_num_value();
  *scalar_type = mf ->get_value_type();
  *total_num_dof = countComponents(f);
  if (*scalar_type) *total_num_dof/=2;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_exist(FieldID* field_id, int * exist)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    *exist = 0;
  else
    *exist = 1;
  return M3DC1_SUCCESS;
}

//*******************************************************
void m3dc1_field_synchronize(apf::Field* f)
//*******************************************************
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e;       

  int num_dof, n = countComponents(f);
  double* sender_data = new double[n];
  double* dof_data = new double[n]; 

  PCU_Comm_Begin();

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (!is_ent_original(m,e) || !m->isShared(e))
      continue;
    getComponents(f, e, 0, dof_data);

    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&(dof_data[0]),n*sizeof(double));
    }
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* r;
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&(sender_data[0]),n*sizeof(double));
      setComponents(f, r, 0, sender_data);
    }
  delete [] dof_data;
  delete [] sender_data;
}

//*******************************************************
int m3dc1_field_sync (FieldID* /* in */ field_id)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field_synchronize((*m3dc1_mesh::instance()->field_container)[*field_id]->get_field());
#ifdef DEBUG
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

// send non-owned copies' dof to owner copy and add them up
//*******************************************************
void m3dc1_field_accumulate(apf::Field* f)
//*******************************************************
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  apf::MeshEntity* e;       

  int num_dof, own_partid, n = countComponents(f);
  double* dof_data = new double[n];
  double* sender_data = new double[n];
  apf::MeshEntity* own_e;
  apf::MeshEntity* r;
  std::map<apf::MeshEntity*, std::vector<double> > save_map;

  PCU_Comm_Begin();

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    own_partid=get_ent_ownpartid(m, e);
    if (own_partid==PCU_Comm_Self()) continue;

    own_e = get_ent_owncopy(m, e);

    getComponents(f, e, 0, &(dof_data[0]));
      
    PCU_COMM_PACK(own_partid, own_e);
    PCU_Comm_Pack(own_partid,&(dof_data[0]),n*sizeof(double));
  }
  m->end(it);

  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&(sender_data[0]),n*sizeof(double));
      for (int i = 0; i < n; ++i)
save_map[r].push_back(sender_data[i]);      
    }

  for (std::map<apf::MeshEntity*, std::vector<double> >::iterator mit=save_map.begin(); mit!=save_map.end(); ++mit)
  {
    e = mit->first;
    getComponents(f, e, 0, dof_data);
    int num_data = mit->second.size()/n;
    for (int i=0; i<num_data;++i)
    {
      for (int j=0; j<n; ++j)
dof_data[j] += mit->second[i*n+j];
    }
    setComponents(f, e, 0, dof_data);
  } 
  delete [] dof_data;
  delete [] sender_data;
}

//*******************************************************
int m3dc1_field_sum (FieldID* /* in */ field_id)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field_accumulate((*m3dc1_mesh::instance()->field_container)[*field_id]->get_field());
  m3dc1_field_synchronize((*m3dc1_mesh::instance()->field_container)[*field_id]->get_field());
#ifdef DEBUG
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_sumsq (FieldID* /* in */ field_id, double* /* out */ sum)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  
  apf::Field* f = (*(m3dc1_mesh::instance()->field_container))[*field_id]->get_field();

#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  *sum=0.;
  int num_dof = countComponents(f);

  double* dof_data= new double[num_dof];
  apf::MeshEntity* e;
  apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
  while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
  {
    if (!is_ent_original(m3dc1_mesh::instance()->mesh,e)) continue;
    getComponents(f, e, 0, dof_data);
    for(int i=0; i<num_dof; ++i)
      *sum+=dof_data[i]*dof_data[i];
  }
  m3dc1_mesh::instance()->mesh->end(it);
  delete [] dof_data;
  return M3DC1_SUCCESS;
}

/** field dof functions */
//*******************************************************
int m3dc1_field_getlocaldofid (FieldID* field_id, 
         int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  *start_dof_id=0;
  *end_dof_id_plus_one=num_dof*m3dc1_mesh::instance()->num_local_ent[0];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getowndofid (FieldID* field_id, 
         int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();

  int num_own_ent = m3dc1_mesh::instance()->num_own_ent[0];
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  
  int start_id = num_own_ent;
  PCU_Exscan_Ints(&start_id,1);

  *start_dof_id=start_id*num_dof;
  *end_dof_id_plus_one=*start_dof_id+num_own_ent*num_dof;
  return M3DC1_SUCCESS;
}
 
//******************************************************* 
int m3dc1_field_getglobaldofid ( FieldID* field_id, 
         int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  assert(mf->get_num_value()*mf->get_dof_per_value()==num_dof);  

  *start_dof_id=0;
  *end_dof_id_plus_one=*start_dof_id+num_dof*m3dc1_mesh::instance()->num_global_ent[0];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getnumlocaldof (FieldID* field_id, int* /* out */ num_local_dof)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  *num_local_dof = (m3dc1_mesh::instance()->num_local_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getnumowndof (FieldID* field_id, int* /* out */ num_own_dof)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  *num_own_dof = (m3dc1_mesh::instance()->num_own_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getnumglobaldof (FieldID* field_id, int* /* out */ num_global_dof)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  *num_global_dof = (m3dc1_mesh::instance()->num_global_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_getdataptr (FieldID* field_id, double** pts)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  apf::Field* f = (*(m3dc1_mesh::instance()->field_container))[*field_id]->get_field();
  if (!isFrozen(f)) freeze(f);
  *pts=getArrayData(f);
  return M3DC1_SUCCESS;
}

// add field2 to field1
//*******************************************************
int m3dc1_field_add(FieldID* /*inout*/ field_id1, FieldID* /*in*/ field_id2)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id1, &isnan);
  assert(isnan==0);
  m3dc1_field_isnan(field_id2, &isnan);
  assert(isnan==0);
#endif
  m3dc1_field * mf1 = (*(m3dc1_mesh::instance()->field_container))[*field_id1];
  m3dc1_field * mf2 = (*(m3dc1_mesh::instance()->field_container))[*field_id2];
  int dofPerEnt1 = mf1->get_num_value()*mf1->get_dof_per_value(); 
  int dofPerEnt2 = mf2->get_num_value()*mf2->get_dof_per_value();
  assert(mf1->get_value_type()==mf2->get_value_type());
  std::vector<double> dofs1(dofPerEnt1*(1+mf1->get_value_type())), dofs2(dofPerEnt2*(1+mf2->get_value_type()));
  int dofMin = std::min(dofPerEnt1,dofPerEnt2);
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  int dofPerEntDummy[2];
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id1, dofPerEntDummy, &dofs1[0]);
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id2, dofPerEntDummy+1, &dofs2[0]);
    for(int i=0; i<dofMin*(1+mf1->get_value_type()); i++)
      dofs1.at(i)+=dofs2.at(i);
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id1, &dofPerEnt1, &dofs1[0]);
  } 
#ifdef DEBUG
  m3dc1_field_isnan(field_id1, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_mult(FieldID* /*inout*/ field_id, double* fac, int * scalar_type)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  double dofs[FIXSIZEBUFF], dofsNew[FIXSIZEBUFF];
  m3dc1_field* mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  int value_type = mf->get_value_type();
  assert(dofPerEnt<=sizeof(dofs)/sizeof(double)*(1+value_type));
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, dofs);
    if(*scalar_type==0)
    {
      for(int i=0; i<dofPerEnt*(1+value_type); i++)
dofsNew[i]=*fac*dofs[i];
    }
    else
    {
      for(int i=0; i<dofPerEnt; i++)
      {
dofsNew[2*i]=fac[0]*dofs[2*i]-fac[1]*dofs[2*i+1];
dofsNew[2*i+1]=fac[0]*dofs[2*i+1]+fac[1]*dofs[2*i];
      }
    }
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofsNew[0]);
  }
#ifdef DEBUG
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_assign(FieldID* /*inout*/ field_id, double* fac, int * scalar_type)
//*******************************************************
{
  m3dc1_field* mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  if (dofPerEnt==0) return M3DC1_FAILURE;

  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  std::vector<double> dofs(dofPerEnt*(1+mf->get_value_type()), fac[0]);
  if(*scalar_type)
    for(int i=0; i<dofPerEnt; i++)
      dofs.at(2*i+1)=fac[1];
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofs[0]);
  }
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_copy(FieldID* /* out */ field_id1, FieldID* /* in */ field_id2)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id1, &isnan);
  assert(isnan==0);
#endif
  m3dc1_field * mf1 = (*(m3dc1_mesh::instance()->field_container))[*field_id1];
  m3dc1_field * mf2 = (*(m3dc1_mesh::instance()->field_container))[*field_id2];
  apf::Field* f1 =  mf1->get_field();
  apf::Field* f2 =  mf2->get_field();
  int dofPerEnt1 = mf1->get_num_value()*mf1->get_dof_per_value();
  int dofPerEnt2 = mf2->get_num_value()*mf2->get_dof_per_value();
  assert(mf1->get_value_type()==mf2->get_value_type());
  std::vector<double> dofs1(dofPerEnt1*(1+mf1->get_value_type())), dofs2(dofPerEnt2*(1+mf2->get_value_type()));
  int dofMin = std::min(dofPerEnt1,dofPerEnt2);
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  int dofPerEntDummy[2];
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id1, dofPerEntDummy, &dofs1[0]);
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id2, dofPerEntDummy+1, &dofs2[0]);
    for(int i=0; i<dofMin*(1+mf1->get_value_type()); i++)
      dofs1.at(i)=dofs2.at(i);
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id1, &dofPerEnt1, &dofs1[0]);
  }
#ifdef DEBUG
  m3dc1_field_isnan(field_id2, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_retrieve (FieldID* /* in */ field_id, double * /*out*/ data, int * /* in */size)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  //m3dc1_field_print(field_id);
  int num_local_dof, num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
  m3dc1_field_getnumlocaldof (field_id, &num_local_dof);

  double* pts=NULL;
  m3dc1_field_getdataptr (field_id, &pts);
  memcpy(data, pts, *size*(1+value_type)*sizeof(double));
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_set (FieldID* /* in */ field_id, double * /*in*/ data, int * /* in */size)
//*******************************************************
{
  int num_local_dof, num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
  m3dc1_field_getnumlocaldof (field_id, &num_local_dof);

  double* pts=NULL;
  m3dc1_field_getdataptr (field_id, &pts);
  memcpy(pts, data, *size*(1+value_type)*sizeof(double));
#ifdef DEBUG
  int isnan;
  m3dc1_field_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_insert(FieldID* /* in */ field_id, int /* in */ * local_dof, 
         int * /* in */ size, double* /* in */ values, int * type, int * op)
//*******************************************************
{
  int num_local_dof, num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
#ifdef DEBUG
  m3dc1_field_getnumlocaldof (field_id, &num_local_dof);
  assert(*local_dof<num_local_dof);
  if(!value_type) assert(!(*type)); // can not insert complex value to real vector
  for(int i=0; i<*size*(1+(*type)); i++)
  {
#ifdef REPLACENANWITHZERO
    if(values[i]!=values[i]) values[i]=0;
#else
    assert(values[i]==values[i]);
#endif
  }
#endif
  std::vector<double> values_convert(*size*(1+value_type),0);
  if(!(*type)&&value_type) // real into complex
  {
    for(int i=0; i<*size; i++)
    {
      values_convert.at(2*i)=values[i];
    }
  }
  else
  {
    for(int i=0; i<*size*(1+value_type); i++)
    {
      values_convert.at(i)=values[i];
    }
  }
  double * dataptr;
  int ibegin=*local_dof*(1+value_type);
  m3dc1_field_getdataptr(field_id, &dataptr);
  if(*op==0) // set value
  {
   for(int i=0; i<*size*(1+value_type); i++)
     dataptr[ibegin+i]=values_convert.at(i);
  }
  else
  {
    for(int i=0; i<*size*(1+value_type); i++)
      dataptr[ibegin+i]+=values_convert[i];
  }
  return M3DC1_SUCCESS;
}
#define FIELDVALUELIMIT 1e100
bool value_is_nan(double val)
{
  return val!=val ||fabs(val) >FIELDVALUELIMIT;
}

//*******************************************************
int m3dc1_field_isnan(FieldID* /* in */ field_id, int * isnan)
//*******************************************************
{
  *isnan=0;
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  int dofPerEnt;
  double dofs[FIXSIZEBUFF];
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, dofs);
    for(int i=0; i<dofPerEnt; i++)
      if(value_is_nan(dofs[i])) 
        *isnan=1;
  }
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_write(FieldID* field_id, const char* file_name)
//*******************************************************
{ 
  m3dc1_field* mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  if (dofPerEnt==0) return M3DC1_FAILURE;
  FILE * fp =fopen(file_name, "w");
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  std::vector<double> dofs(dofPerEnt*(1+mf->get_value_type()),0);
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofs[0]);
    for(int i=0; i<dofPerEnt; i++)
    {
      for(int j=0; j<1+mf->get_value_type(); j++)
        fprintf(fp, "%lf ", dofs[i*(1+mf->get_value_type())+j]);
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
}


//*******************************************************
int m3dc1_field_compare(FieldID* field_id_1, FieldID* field_id_2)
//*******************************************************
{
  apf::Field* f_1 = (*(m3dc1_mesh::instance()->field_container))[*field_id_1]->get_field();
  double* field_data_1 =getArrayData(f_1);

  apf::Field* f_2 = (*(m3dc1_mesh::instance()->field_container))[*field_id_2]->get_field();
  double* field_data_2 =getArrayData(f_2);

  int num_dof_1=countComponents(f_1);
  int num_dof_2=countComponents(f_2);
  if (num_dof_1!=num_dof_2) 
  {
    if (!PCU_Comm_Self()) 
      cout<<"[M3DC1 INFO] "<<__func__<<": #dof mismatch "<<getName(f_1)
          <<"- "<<num_dof_1<<", "<<getName(f_2)<<"- "<<num_dof_2<<"\n";
    return M3DC1_FAILURE;
  }
  int ierr = M3DC1_SUCCESS;
  for (int i=0; i<num_dof_1*m3dc1_mesh::instance()->num_local_ent[0]; ++i)
  {  if (!m3dc1_double_isequal(field_data_1[i], field_data_2[i])) 
    {
     cout<<"[M3DC1 ERROR] "<<__func__<<": "<<getName(f_1)<<"["<<i<<"]="<<field_data_1[i]
          <<", "<<getName(f_2)<<"["<<i<<"]="<<field_data_2[i]<<"\n";
      ierr=M3DC1_FAILURE;
      break;
    }
  }
  int global_ierr;
  MPI_Allreduce(&ierr, &global_ierr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
  if (global_ierr==M3DC1_FAILURE)
  {
    if (!PCU_Comm_Self())
      cout<<"[M3DC1 INFO] "<<__func__<<": dof value mismatch of fields "<<getName(f_1)
          <<" and "<<getName(f_2)<<"\n";
    
    return M3DC1_FAILURE;
  }
  else
    return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_field_print(FieldID* field_id)
//*******************************************************
{ 
  apf::Field* f = (*(m3dc1_mesh::instance()->field_container))[*field_id]->get_field();
  double* field_data =getArrayData(f);
  
  apf::MeshEntity* e;
  int global_ent_id, num_dof=countComponents(f);
  double* dof_data = new double[num_dof];

  switch (num_dof)
  {
    case 1: {
              int i=0;
              apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
              while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<"]\n";
               }
               m3dc1_mesh::instance()->mesh->end(it);
               break;
            }
    case 2: {
              int i=0;
              apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
              while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<"]\n";
               }
               m3dc1_mesh::instance()->mesh->end(it);
               break;
            }
    case 3: {
              int i=0;
              apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
              while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<"]\n";
               }
               m3dc1_mesh::instance()->mesh->end(it);
               break;
            }
    case 4: {
              int i=0;
              apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
              while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<"]\n";
               }
               m3dc1_mesh::instance()->mesh->end(it);
               break;
            }
    case 6: {
              int i=0;
              apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
              while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<", "<<dof_data[4]
               <<", "<<dof_data[5]
               <<"]\n";
               }
               m3dc1_mesh::instance()->mesh->end(it);
               break;
          }
    case 8: {
              int i=0;
              apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
              while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<", "<<dof_data[4]
               <<", "<<dof_data[5]
               <<", "<<dof_data[6]
               <<", "<<dof_data[7]
               <<"]\n";
               }
               m3dc1_mesh::instance()->mesh->end(it);
               break;
            }
    case 12: {
              int i=0;
              apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
              while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
  //            for (int i=0; i<m3dc1_mesh::instance()->num_local_ent[0]; ++i)
              { 
                assert(e ==getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i++));
               // if (m3dc1_mesh::instance()->mesh->getOwner(e)!=PCU_Comm_Self()) continue;
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<", "<<dof_data[4]
               <<", "<<dof_data[5]
               <<", "<<dof_data[6]
               <<", "<<dof_data[7]
               <<", "<<dof_data[8]
               <<", "<<dof_data[9]
               <<", "<<dof_data[10]
               <<", "<<dof_data[11]
               <<"]\n";
               }
               m3dc1_mesh::instance()->mesh->end(it);
               break;
             }
    case 18: {
              int i=0;
              apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
              while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<", "<<dof_data[4]
               <<", "<<dof_data[5]
               <<", "<<dof_data[6]
               <<", "<<dof_data[7]
               <<", "<<dof_data[8]
               <<", "<<dof_data[9]
               <<", "<<dof_data[10]
               <<", "<<dof_data[11]
               <<", "<<dof_data[12]
               <<", "<<dof_data[13]
               <<", "<<dof_data[14]
               <<", "<<dof_data[15]
               <<", "<<dof_data[16]
               <<", "<<dof_data[17]
               <<"]\n";
               }
               m3dc1_mesh::instance()->mesh->end(it);
               break;
             }
    case 24: {
              int i=0;
              apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
              while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<", "<<dof_data[4]
               <<", "<<dof_data[5]
               <<", "<<dof_data[6]
               <<", "<<dof_data[7]
               <<", "<<dof_data[8]
               <<", "<<dof_data[9]
               <<", "<<dof_data[10]
               <<", "<<dof_data[11]
               <<", "<<dof_data[12]
               <<", "<<dof_data[13]
               <<", "<<dof_data[14]
               <<", "<<dof_data[15]
               <<", "<<dof_data[16]
               <<", "<<dof_data[17]
               <<", "<<dof_data[18]
               <<", "<<dof_data[19]
               <<", "<<dof_data[20]
               <<", "<<dof_data[21]
               <<", "<<dof_data[22]
               <<", "<<dof_data[23]
               <<"]\n";
               }
               m3dc1_mesh::instance()->mesh->end(it);
               break;
            }
      default: if (!PCU_Comm_Self()) std::cout<<__func__<<" failed for field "<<getName(f)<<": does support "<<num_dof<<" dofs\n";
              return M3DC1_FAILURE;
  } // switch
  return M3DC1_SUCCESS;
}
//*******************************************************
int m3dc1_ent_getlocaldofid(int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
                       int* /* out */ start_dof_id, int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (*ent_dim!=0)
    return M3DC1_FAILURE;

  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  if (!e)
    return M3DC1_FAILURE;

  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int num_val =  mf->get_num_value();
  int dof_per_value = mf->get_dof_per_value();
  int dof_per_node = dof_per_value * num_val;
  *start_dof_id = *ent_id*dof_per_node;
  *end_dof_id_plus_one = *start_dof_id +dof_per_node;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getglobaldofid (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
         int* /* out */ start_global_dof_id, int* /* out */ end_global_dof_id_plus_one)
//*******************************************************
{
  if (*ent_dim!=0)
    return M3DC1_FAILURE;

  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);

  if (!e)
    return M3DC1_FAILURE;

  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;

  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int num_val =  mf->get_num_value();
  int dof_per_value = mf->get_dof_per_value();
  int dof_per_node = dof_per_value * num_val;

  int global_id = getNumber(get_global_numbering(), e, 0, 0);
  *start_global_dof_id = global_id*dof_per_node;
  *end_global_dof_id_plus_one =*start_global_dof_id + dof_per_node;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getnumdof (int* /* in */ ent_dim, int* /* in */ ent_id, 
         FieldID* field_id, int* /* out */ num_dof)
//*******************************************************
{
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  if (*ent_dim!=0)
    return M3DC1_FAILURE;

  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  *num_dof =  mf->get_num_value() * mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_setdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
                          int* /* out */ num_dof, double* dof_data)
//*******************************************************
{
  assert(*ent_dim==0);
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  if (!e)
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();

#ifdef DEBUG
  assert(countComponents(f)==*num_dof*(1+mf->get_value_type()));
  for(int i=0; i<*num_dof*(1+mf->get_value_type()); i++)
    assert(!value_is_nan(dof_data[i]));
#endif
  setComponents(f, e, 0, dof_data);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_ent_getdofdata (int* /* in */ ent_dim, int* /* in */ ent_id, FieldID* field_id, 
                          int* /* out */ num_dof, double* dof_data)
//*******************************************************
{
  assert(*ent_dim==0);
  apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, *ent_dim, *ent_id);
  if (!e)
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();

  getComponents(f, e, 0, dof_data);
  *num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
#ifdef DEBUG
  for(int i=0; i<*num_dof*(1+mf->get_value_type()); i++)
    assert(!value_is_nan(dof_data[i]));
  int start_dof_id,end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid(ent_dim, ent_id,field_id, &start_dof_id, &end_dof_id_plus_one);
  double* data;
  m3dc1_field_getdataptr(field_id, &data);
  int start=start_dof_id*(1+mf->get_value_type());
  for( int i=0; i< *num_dof; i++)
    assert(data[start++]==dof_data[i]);
#endif
  return M3DC1_SUCCESS;
}

#ifndef M3DC1_MESHGEN
/** matrix and solver functions */
std::map<int, int> matHit;
int getMatHit(int id) { return matHit[id];};
void addMatHit(int id) { matHit[id]++; }

//*******************************************************
int m3dc1_matrix_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID *field_id)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);

  if (mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" already created\n";
    return M3DC1_FAILURE; 
  }
  // check field exists
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: field with id "<<*field_id<<" doesn't exist\n";
    return M3DC1_FAILURE; 
  }

  if (*matrix_type==M3DC1_MULTIPLY) // matrix for multiplication
  {
    matrix_mult* new_mat = new matrix_mult(*matrix_id, *scalar_type, *field_id);
    m3dc1_solver::instance()->add_matrix(*matrix_id, (m3dc1_matrix*)new_mat);
  }
  else 
  {
    matrix_solve* new_mat= new matrix_solve(*matrix_id, *scalar_type, *field_id);
    m3dc1_solver::instance()->add_matrix(*matrix_id, (m3dc1_matrix*)new_mat);
  }

  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_freeze(int* matrix_id) 
//*******************************************************
{
  double t1 = MPI_Wtime();
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  mat->assemble();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_delete(int* matrix_id)
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  typedef std::map<int, m3dc1_matrix*> matrix_container_map;
  m3dc1_solver::instance()->matrix_container->erase(matrix_container_map::key_type(*matrix_id));
  delete mat;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_insert(int* matrix_id, int* row, 
         int* col, int* scalar_type, double* val)
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  if (mat->get_status()==M3DC1_FIXED)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" is fixed\n";
    return M3DC1_FAILURE;
  }
  if (*scalar_type==M3DC1_COMPLEX)
    mat->set_value(*row, *col, INSERT_VALUES, val[0], val[1]);
  else
    mat->set_value(*row, *col, INSERT_VALUES, *val, 0);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_add (int* matrix_id, int* row, int* col, 
                      int* scalar_type, double* val) //globalinsertval_
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  if (mat->get_status()==M3DC1_FIXED)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" is fixed\n";
    return M3DC1_FAILURE;
  }
  if (*scalar_type==M3DC1_COMPLEX)
    mat->set_value(*row, *col, ADD_VALUES, val[0], val[1]);
  else
    mat->set_value(*row, *col, ADD_VALUES, *val, 0);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_matrix_setbc(int* matrix_id, int* row)
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  if (mat->get_type()!=M3DC1_SOLVE)
  { 
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for multiplication (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  int field = mat->get_fieldOrdering();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
#ifdef DEBUG
  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
#endif
  int row_g = start_global_dof_id+*row%total_num_dof;
  //std::cout<<PCU_Comm_Self()<<" m3dc1_matrix_setbc "<<*matrix_id<<" row "<<*row<<" row_g "<<row_g<<std::endl;
  (dynamic_cast<matrix_solve*>(mat))->set_bc(row_g);
}

//*******************************************************
int m3dc1_matrix_setlaplacebc(int * matrix_id, int *row,
         int * numVals, int *columns, double * values)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  if (mat->get_type()!=M3DC1_SOLVE)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for multiplication (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  std::vector <int> columns_g(*numVals);
  int field = mat->get_fieldOrdering();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
#ifdef DEBUG
  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
  for (int i=0; i<*numVals; i++)
    assert(columns[i]>=start_dof_id&&columns[i]<end_dof_id_plus_one);
#endif
  int row_g = start_global_dof_id+*row%total_num_dof;
  for(int i=0; i<*numVals; i++)
  {
    columns_g.at(i) = start_global_dof_id+columns[i]%total_num_dof;
  }
  (dynamic_cast<matrix_solve*>(mat))->set_row(row_g, *numVals, &columns_g[0], values);
}

int m3dc1_matrix_solve(int* matrix_id, FieldID* rhs_sol) //solveSysEqu_
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  if (mat->get_type()!=M3DC1_SOLVE)
  { 
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for multiplication (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  if (!PCU_Comm_Self())
     std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<", field "<<*rhs_sol<<"\n";

  (dynamic_cast<matrix_solve*>(mat))->solve(*rhs_sol);

  addMatHit(*matrix_id);
}

//*******************************************************
int m3dc1_matrix_multiply(int* matrix_id, FieldID* inputvecid, 
         FieldID* outputvecid) 
//*******************************************************
{  
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  if (mat->get_type()!=M3DC1_MULTIPLY)
  { 
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported with matrix for solving (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
#ifdef DEBUG_
  int isnan;
  if(!PCU_Comm_Self()) std::cout<<" mult matrix "<<*matrix_id<<" with input vec "<<*inputvecid<<" output vec "<<*outputvecid<<std::endl;
  m3dc1_field_isnan(inputvecid, &isnan);
  assert(isnan==0);
#endif

#ifdef PRINTSOLVEMATRIX
  char filename[256];
  int time=getMatHit(*matrix_id);
  sprintf(filename, "mat%din%d.m",*matrix_id, time);
  m3dc1_field_write(inputvecid, filename);
#endif
  (dynamic_cast<matrix_mult*>(mat))->multiply(*inputvecid, *outputvecid);
#ifdef PRINTSOLVEMATRIX
  sprintf(filename, "mat%dout%d.m",*matrix_id, time);
  m3dc1_field_write(outputvecid, filename);
#endif
#ifdef DEBUG_
  m3dc1_field_isnan(outputvecid, &isnan);
  assert(isnan==0);
#endif
  addMatHit(*matrix_id);
}

//*******************************************************
int m3dc1_matrix_flush(int* matrix_id)
//*******************************************************
{
  double t1 = MPI_Wtime();
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  mat->flushAssembly();
}

//*******************************************************
int m3dc1_matrix_getiternum(int* matrix_id, int * iter_num)
//*******************************************************
{ 
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  *iter_num = dynamic_cast<matrix_solve*> (mat)->iterNum;
}

//*******************************************************
int m3dc1_matrix_insertblock(int* matrix_id, int * ielm, 
          int* rowIdx, int * columnIdx, double * values)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  int field = mat->get_fieldOrdering();
  // need to change later, should get the value from field calls ...
  int dofPerVar = 6;
  char field_name[256];
  int num_values, value_type, total_num_dof; 
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  dofPerVar=total_num_dof/num_values;
  int nodes[6];
  int ent_dim=0;
  int ielm_dim = 2;
  int nodes_per_element=sizeof(nodes)/sizeof(int), nodes_per_element_get;

  if (m3dc1_mesh::instance()->mesh->getDimension()==3) ielm_dim =3;
  m3dc1_ent_getadj (&ielm_dim, ielm, &ent_dim, nodes, &nodes_per_element, &nodes_per_element_get);
  nodes_per_element=nodes_per_element_get;
  int start_global_dof_id,end_global_dof_id_plus_one;
  int start_global_dof,end_global_dof_id;
  // need to change later, should get the value from field calls ...
  int scalar_type = mat->get_scalar_type();
  assert(scalar_type==value_type);
  int numDofs = total_num_dof;
  int numVar = numDofs/dofPerVar;
  assert(*rowIdx<numVar && *columnIdx<numVar);
  int rows[1024], columns[1024];
  assert(sizeof(rows)/sizeof(int)>=dofPerVar*nodes_per_element);
  if(mat->get_type()==0)
  {
    int localFlag=0;
    matrix_mult* mmat = dynamic_cast<matrix_mult*> (mat);
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      if(mmat->is_mat_local()) m3dc1_ent_getlocaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      else m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      for(int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnIdx)*dofPerVar+i;
      }
    }
    mmat->add_values(dofPerVar*nodes_per_element, rows,dofPerVar*nodes_per_element, columns, values);
  }
  else
  {
    matrix_solve* smat = dynamic_cast<matrix_solve*> (mat);
    int nodeOwner[6];
    int columns_bloc[6], rows_bloc[6];
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      m3dc1_ent_getownpartid (&ent_dim, nodes+inode, nodeOwner+inode);
      m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      rows_bloc[inode]=nodes[inode]*numVar+*rowIdx;
      columns_bloc[inode]=nodes[inode]*numVar+*columnIdx;
      for(int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnIdx)*dofPerVar+i;
      }
    }
    int numValuesNode = dofPerVar*dofPerVar*nodes_per_element*(1+scalar_type);
    int offset=0;
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      if(nodeOwner[inode]!=PCU_Comm_Self()&&!m3dc1_solver::instance()->assembleOption)
        smat->add_blockvalues(1, rows_bloc+inode, nodes_per_element, columns_bloc, values+offset);
      else 
        smat->add_values(dofPerVar, rows+dofPerVar*inode, dofPerVar*nodes_per_element, columns, values+offset);
      offset+=numValuesNode;
    }
  }
}

//*******************************************************
int m3dc1_matrix_write(int* matrix_id, const char* file_name)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  char new_name[256];
  sprintf(new_name, "%s%d",file_name, *matrix_id);
  mat->write(new_name);
}


//*******************************************************
int m3dc1_matrix_print(int* matrix_id)
//*******************************************************
{
  m3dc1_matrix* mat = m3dc1_solver::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  int row, col, csize, sum_csize=0, index=0;

  vector<int> rows;
  vector<int> n_cols;
  vector<int> cols;
  vector<double> vals;

  mat->get_values(rows, n_cols, cols, vals);
  for (int i=0; i<rows.size(); ++i)
    sum_csize += n_cols[i];
  assert(vals.size()==sum_csize);

  if (!PCU_Comm_Self()) 
    std::cout<<"[M3DC1 INFO] "<<__func__<<": printing matrix "<<*matrix_id<<"\n";

  for (int i=0; i<rows.size(); ++i)
  {
    row = rows[i];
    csize = n_cols[i];
    for (int j=0; j<csize; ++j)
    {
      std::cout<<"["<<PCU_Comm_Self()<<"]\t"<<row<<"\t"<<cols[index]<<"\t"<<vals[index]<<"\n";
      ++index;
    }
  }
  assert(index == vals.size());
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_matrix_setassembleoption(int * op)
//*******************************************************
{
  m3dc1_solver::instance()->assembleOption= *op;
}

//*******************************************************
int m3dc1_field_sum_plane (FieldID* /* in */ field_id)
//*******************************************************
{
  MPI_Comm icomm= m3dc1_model::instance()->getMPICommPlane();
  int num_vtx,num_dof=0, vertex_type=0;
  m3dc1_field_getnumlocaldof(field_id, &num_dof);
  char field_name[256];
  int num_values, value_type, total_num_dof;
  m3dc1_field_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);

  m3dc1_mesh_getnument (&vertex_type, &num_vtx);

  int data_size=total_num_dof*num_vtx*(1+value_type);
  assert(total_num_dof*num_vtx==num_dof);
  double * thevec = NULL;
  m3dc1_field_getdataptr(field_id, &thevec);
  double * sendbuf = new double [data_size];
  m3dc1_field_retrieve (field_id, sendbuf, &num_dof);
  MPI_Allreduce (sendbuf,thevec,data_size,MPI_DOUBLE,MPI_SUM,icomm) ;
  m3dc1_field_sync(field_id);
  delete []sendbuf;
}

int adapt_time=0;
void group_complex_dof (apf::Field* field, int option);
int adapt_by_field (int * fieldId, double* psi0, double * psil)
{
  FILE *fp = fopen("sizefieldParam", "r");
  if(!fp)
  {
    std::cout<<" file sizefieldParam not found "<<std::endl;
    throw 1;
  }
  double param[13];
  set<int> field_keep;
  field_keep.insert(*fieldId);
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  // read the size field parameters
  for(int i=0; i<13; i++)
    fscanf(fp, "%lf ", &param[i]);
  fclose(fp);
  apf::Field* psiField = (*(m3dc1_mesh::instance()->field_container))[*fieldId]->get_field();

  int node_glb_order=NODE_GLB_ORDER;
  m3dc1_field_delete (&node_glb_order);

  // delete all the matrix
  std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }

  int valueType = (*(m3dc1_mesh::instance()->field_container))[*fieldId]->get_value_type();
  SizeFieldPsi sf (psiField, *psi0, *psil, param, valueType);
  double mmax[2], mmin[2];
  m3dc1_model_getmaxcoord(mmax,mmax+1);
  m3dc1_model_getmincoord(mmin,mmin+1);

  ReducedQuinticImplicit shape;
  vector<apf::Field*> fields;
  std::map<FieldID, m3dc1_field*> :: iterator it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    assert(valueType==complexType);
    if(complexType) group_complex_dof(field, 1);
    if(isFrozen(field)) unfreeze(field);
    if(!PCU_Comm_Self()) std::cout<<"Solution transfer: add field "<<apf::getName(field)<<std::endl;
    fields.push_back(field);
    it++;
  }
  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    if(!PCU_Comm_Self()) std::cout<<"destroying numbering "<<getName(n)<<endl;
    apf::destroyNumbering(n);
  }
  ReducedQuinticTransfer slnTrans(mesh,fields, &shape);
  ma::Input* in = ma::configure(mesh,&sf,&slnTrans);
  in->maximumIterations = 9;
  in->shouldRunPostZoltan = true;
  ma::adapt(in);
  reorderMdsMesh(mesh);

  m3dc1_mesh::instance()->initialize();
  it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if(complexType) group_complex_dof(field, 0);
    if(!isFrozen(field)) freeze(field);
#ifdef DEBUG
    int isnan;
    int fieldId= it->first; 
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    m3dc1_field_synchronize(field);
    get_global_numbering();
#ifdef DEBUG
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    it++;
  }
}

double absSize[2]={0,1}, relSize[2]={0.3, 1.5};
int set_mesh_size_bound (double* abs_size, double * rel_size)
{
  for(int i=0; i<2; i++)
  {
    absSize[i] =  abs_size[i];
    relSize[i] = rel_size[i];
  }
}

void smooth_size_field (apf::Field* sizeField)
{
  int entDim=0, numVert=0;
  m3dc1_mesh_getnument (&entDim, &numVert);
  for(int i=0; i<numVert; i++)
  {
    vector<apf::MeshEntity*> nodes;
    apf::MeshEntity* e = apf:: getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
    assert(e);
    nodes.push_back(e);
    double sizeOrg=0;
    getComponents(sizeField, e, 0, &sizeOrg);
    apf::Adjacent adjacent;
    m3dc1_mesh::instance()->mesh->getAdjacent(e,1,adjacent);
    for(int i=0; i<adjacent.getSize(); i++)
    {
      apf::Downward downward;
      m3dc1_mesh::instance()->mesh->getDownward(adjacent[i], 0, downward);
      nodes.push_back(downward[0]==e?downward[1]:downward[0]);
    }
    double size=0;
    for(int i=0; i<nodes.size(); i++)
    {
      double buff;
      getComponents(sizeField, nodes[i], 0, &buff);
      size+=1./buff;
    }
    size/=nodes.size();
    size=1./size;
    setComponents(sizeField, e, 0, &size);
  }
}

void group_complex_dof (apf::Field* field, int option)
{
  int num_dof_double = countComponents(field);
  assert(num_dof_double/6%2==0);
  int num_dof = num_dof_double/2;
  vector<double> dofs(num_dof_double);
  vector<double> newdofs(num_dof_double);
  int entDim=0, numVert=0;
  m3dc1_mesh_getnument (&entDim, &numVert);
  for(int i=0; i<numVert; i++)
  {
    apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
    getComponents(field, e, 0, &(dofs[0]));
    for(int j=0; j<num_dof/6; j++)
    {
      if(option)
      {
        for(int k=0; k<6; k++)
        {
          newdofs.at(2*j*6+k)=dofs.at(2*j*6+2*k);
          newdofs.at(2*j*6+6+k)=dofs.at(2*j*6+2*k+1);
        }
      }
      else
      {
        for(int k=0; k<6; k++)
        {
          newdofs.at(2*j*6+2*k)=dofs.at(2*j*6+k);
          newdofs.at(2*j*6+2*k+1)=dofs.at(2*j*6+6+k);
        }
      }
    }
    setComponents(field, e, 0, &(newdofs[0]));
  }
}

double p=4;
int set_adapt_p (double * pp) {p=*pp;}
int adapt_by_error_field (double * errorData, double * errorAimed, int * max_adapt_node, int * option)
{
  apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
  apf::Field* sizeField = createPackedField(m3dc1_mesh::instance()->mesh, "size_field", 1);
  SizeFieldError sf (m3dc1_mesh::instance()->mesh, sizeField, *errorAimed);

  int entDim=0, numVert=0;
  m3dc1_mesh_getnument (&entDim, &numVert);
  // first sum error ^ (2d/(2p+d))
  double d=2;
  double errorSum=0;
  if(*option)
  {
    for(int i=0; i<numVert; i++)
    {
      if(is_ent_original(m3dc1_mesh::instance()->mesh,getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i)))
        errorSum+=pow(errorData[i],d/(p+d/2.0));
    }
    double errorSumBuff=errorSum;
    MPI_Allreduce(&errorSumBuff, &errorSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    errorSum = *errorAimed*(*errorAimed)/errorSum;
    errorSum = pow(errorSum,1./(2.*p));
  }
  else errorSum=pow(*errorAimed,1./(p+d/2.));
  double size_estimate=0;
  for(int i=0; i<numVert; i++)
  {
    apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
    assert(e);
    if(!is_ent_original(m3dc1_mesh::instance()->mesh,e)) continue;
    double size = sf.getSize(e);
    double targetSize = errorSum*pow(errorData[i],-1./(p+d/2.));
    size_estimate+=max(1.,1./targetSize/targetSize);
  }
  double size_estimate_buff=size_estimate;
  MPI_Allreduce(&size_estimate_buff, &size_estimate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int numNodeGlobl=0, dim=0;
  m3dc1_mesh_getnumglobalent(&dim, &numNodeGlobl);
  cout<<" numVert "<<numNodeGlobl<<" size_estimate "<<size_estimate;
  if(size_estimate>*max_adapt_node) errorSum*=sqrt(size_estimate/(*max_adapt_node));
  for(int i=0; i<numVert; i++)
  {
    apf::MeshEntity* e =getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
    assert(e);
    double size = sf.getSize(e);
    assert(errorData[i]==errorData[i]);
    double targetSize = errorSum*pow(errorData[i],-1./(p+d/2.));
    if(targetSize>relSize[1]) targetSize=relSize[1]; // not too much coarsening
    if(targetSize<relSize[0]) targetSize=relSize[0]; // not too much refining
    targetSize*=size;
    if(targetSize>absSize[1]) targetSize=absSize[1];
    if(targetSize<absSize[0]) targetSize=absSize[0];
    setComponents(sizeField, e, 0, &targetSize);
  }
  // only implemented for one process
  //if(PCU_Comm_Peers()==1)
  {
    smooth_size_field(sizeField);
    smooth_size_field(sizeField);
  }
  m3dc1_field_synchronize(sizeField);
  int node_glb_order=NODE_GLB_ORDER;
  m3dc1_field_delete (&node_glb_order);

  // delete all the matrix
  std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }

  ReducedQuinticImplicit shape;
  vector<apf::Field*> fields;
  fields.push_back(sizeField);
  std::map<FieldID, m3dc1_field*> :: iterator it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if(complexType) group_complex_dof(field, 1);
    if(isFrozen(field)) unfreeze(field);
    fields.push_back(field);
    it++;
  }
  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    apf::destroyNumbering(n);
  }

  ReducedQuinticTransfer slnTrans(mesh,fields, &shape);
  ma::Input* in = ma::configure(mesh,&sf,&slnTrans);
  in->maximumIterations = 5;
  in->shouldRunPostZoltan = true;
  ma::adapt(in);
  reorderMdsMesh(mesh);
  
  m3dc1_mesh::instance()->initialize();
  it=m3dc1_mesh::instance()->field_container->begin();
  while(it!=m3dc1_mesh::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if(complexType) group_complex_dof(field, 0);
    if(!isFrozen(field)) freeze(field);
#ifdef DEBUG
    int isnan;
    int fieldId= it->first; 
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    m3dc1_field_synchronize(field);
    get_global_numbering();
#ifdef DEBUG
    m3dc1_field_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    it++;
  }
 destroyField(sizeField);
}

#endif // #ifndef M3DC1_MESHGEN

static double square(double x) {return x * x;}

int m3dc1_field_printcompnorm(FieldID* /* in */ field_id, const char* info)
{
  /* assumes DOFs are real */
  double* array=NULL;
  m3dc1_field_getdataptr (field_id, &array);
  int dof_in_array;
  m3dc1_field_getnumlocaldof(field_id, &dof_in_array);
  m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];

  apf::Field* f = mf ->get_field();
  int dof_per_node = countComponents(f);
  int dof_per_value = C1TRIDOFNODE;
  int values_per_node = dof_per_node / dof_per_value;
  vector<double> norms(values_per_node);
  int values_in_array = dof_in_array / dof_per_value;
  int nnodes = dof_in_array / dof_per_node;
  double* p = array;
  for (int i = 0; i < nnodes; ++i)
    for (int j = 0; j < values_per_node; ++j)
      for (int k = 0; k < dof_per_value; ++k)
        norms.at(j) += square(*p++);
  assert(nnodes * values_per_node * dof_per_value == dof_in_array);
  vector<double> buff = norms;
  MPI_Allreduce(&buff[0], &norms[0], values_per_node,
      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int psize;
  MPI_Comm_size(MPI_COMM_WORLD,&psize);
  if(PCU_Comm_Self() == psize-1)
  {
    std::cout<< "norm of vec "<<info;
    for(int i = 0; i < values_per_node; ++i)
      std::cout<<" "<<std::sqrt(norms[i]);
    std::cout<<std::endl;
  }
}

int m3dc1_mesh_write(const char* filename, int *option)
{
  if(*option==0 ||*option==3)
  {
    apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
    apf::MeshEntity* e;
    int dim=2, num_ent;
    m3dc1_mesh_getnument (&dim, &num_ent);
    vector<double> geoId (num_ent);
    apf:: MeshIterator* it = mesh->begin(dim);
    while ((e = mesh->iterate(it)))
    {
      int ent_id = getMdsIndex(m3dc1_mesh::instance()->mesh, e);
      int geom_class_dim,geom_class_id;
      m3dc1_ent_getgeomclass (&dim, &ent_id, &geom_class_dim, &geom_class_id);
      geoId.at(ent_id)=geom_class_id;
    }
    mesh->end(it);
    apf::writeVtkFiles(filename,m3dc1_mesh::instance()->mesh);
    int one=1;
    if(*option==3) output_face_data (&one, &geoId[0], "geoId");
  }
  else
  {
    char filename_buff[256];
    sprintf(filename_buff, "%s.smb",filename);
    int fieldID=12;
    double dofBuff[1024];
    m3dc1_field * mf = (*(m3dc1_mesh::instance()->field_container))[fieldID];
    apf::Mesh2* mesh = m3dc1_mesh::instance()->mesh;
    apf::Field* f = mf ->get_field();
    int numDof = countComponents(f);
    apf::MeshTag* tag = mesh->createDoubleTag("field12", numDof);
    apf::MeshEntity* e;
    const int dim=0;
    apf:: MeshIterator* it = mesh->begin(dim);
    while ((e = mesh->iterate(it)))
    {
      apf::getComponents(f, e, 0, dofBuff);
      mesh->setDoubleTag(e,tag, dofBuff);
    }
    mesh->end(it);
    m3dc1_mesh::instance()->mesh->writeNative(filename_buff);
    apf::removeTagFromDimension(mesh, tag, dim);
    mesh->destroyTag(tag);
  }
}

int sum_edge_data (double * data, int* size)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  int num_edge=0, edg_dim=1;
  m3dc1_mesh_getnument(&edg_dim, &num_edge);
  PCU_Comm_Begin();
  for(int i=0; i<num_edge; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, edg_dim, i);
    int own_partid=get_ent_ownpartid(m, e);
    apf::MeshEntity* own_e = get_ent_owncopy(m, e);
    if (own_partid==PCU_Comm_Self()) continue;
    PCU_COMM_PACK(own_partid, own_e);
    PCU_Comm_Pack(own_partid,&(data[i*(*size)]),(*size)*sizeof(double));
  }

  PCU_Comm_Send();
  double* receive_buff = new double [*size];
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* edge;
      PCU_COMM_UNPACK(edge);
      PCU_Comm_Unpack(receive_buff, (*size)*sizeof(double));
      int iedge = getMdsIndex(m, edge);
      for (int i = 0; i < *size; i++)
        data[iedge*(*size)+i]+=receive_buff[i];
    }

  PCU_Comm_Begin();
  for(int i=0; i<num_edge; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, edg_dim, i);
    if (!is_ent_original(m,e) || !m->isShared(e))
      continue;

    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&(data[i*(*size)]),(*size)*sizeof(double));
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* edge;
      PCU_COMM_UNPACK(edge);
      PCU_Comm_Unpack(receive_buff, (*size)*sizeof(double));
      int iedge = getMdsIndex(m, edge);
      for (int i = 0; i < *size; i++)
        data[iedge*(*size)+i]=receive_buff[i];
    }
   delete []receive_buff;
}

int get_node_error_from_elm (double * elm_data, int * size, double* nod_data)
{
  apf::Mesh2* m = m3dc1_mesh::instance()->mesh;
  int num_node=0, num_elm=0, nod_dim=0, elm_dim=2;
  m3dc1_mesh_getnument(&nod_dim, &num_node);
  m3dc1_mesh_getnument(&elm_dim, &num_elm);
  PCU_Comm_Begin();
  double* buff = new double[*size];
  double* area = new double[num_node];
  for(int i=0; i<num_node; i++)
    area[i]=0.;
  for(int i=0; i<(*size)*num_elm; i++ )
    for(int j=0; j<*size; j++)
    {
      assert(elm_data[(*size)*i+j]==elm_data[(*size)*i+j]);
      assert(elm_data[(*size)*i+j]>=0);
    }
  for(int i=0; i<num_node; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, nod_dim, i);
    int own_partid=get_ent_ownpartid(m, e);
    apf::MeshEntity* own_e = get_ent_owncopy(m, e);
    apf::Adjacent adjacent;
    m->getAdjacent(e,2,adjacent);
    for(int j=0; j<adjacent.getSize(); j++)
    {
       apf::MeshElement* me = createMeshElement(m, adjacent[j]);
       double s = apf::measure(me);
       int ielm = getMdsIndex(m, adjacent[j]);
       assert(ielm>=0 &&ielm<num_elm);
       for(int k=0; k<*size; k++)
          nod_data[i*(*size)+k]+=s*elm_data[(*size)*ielm+k];
       area[i]+=s;
       destroyMeshElement(me);
    }

    if (own_partid==PCU_Comm_Self()) continue;
    PCU_COMM_PACK(own_partid, own_e);
    PCU_COMM_PACK(own_partid, area[i]);
    PCU_Comm_Pack(own_partid, &(nod_data[(*size)*i]), sizeof(double)*(*size));
  }

  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* node;
      PCU_COMM_UNPACK(node);
      int inode = getMdsIndex(m, node);
      double s;
      PCU_COMM_UNPACK(s);
      area[inode]+=s;
      PCU_Comm_Unpack(buff, (*size)*sizeof(double));
      for (int i = 0; i < *size; i++)
        nod_data[inode*(*size)+i]+=buff[i];
    }

  for(int i=0; i<num_node; i++)
  {
    for(int j=0; j<*size; j++)
      nod_data[i*(*size)+j]/=area[i];
  }
  PCU_Comm_Begin();
  for(int i=0; i<num_node; i++)
  {
    apf::MeshEntity* e = getMdsEntity(m, nod_dim, i);
    if (!is_ent_original(m,e) || !m->isShared(e))
      continue;
    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&(nod_data[i*(*size)]),(*size)*sizeof(double));
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* node;
      PCU_COMM_UNPACK(node);
      PCU_Comm_Unpack(buff, (*size)*sizeof(double));
      int inode = getMdsIndex(m, node);
      for (int i = 0; i < *size; i++)
        nod_data[inode*(*size)+i]=buff[i];
    }
  delete []buff;
  delete []area;  
}
int m3dc1_field_max (FieldID* field_id, double * max_val, double * min_val)
{
  m3dc1_field* mf = (*(m3dc1_mesh::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  if(mf->get_value_type()) dofPerEnt *= 2;
  if (dofPerEnt==0) return M3DC1_FAILURE;

  int num_vtx=0;
  int vertex_type=0;
  m3dc1_mesh_getnument (&vertex_type, &num_vtx);
  std::vector<double> maxVal(dofPerEnt, -1e30), minVal(dofPerEnt,1e30), dofs(dofPerEnt);
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofs[0]);
    for(int i=0; i<dofPerEnt; i++)
    {
      if(maxVal[i]<dofs[i]) maxVal[i]=dofs[i];
      if(minVal[i]>dofs[i]) minVal[i]=dofs[i];
    }
  }
  MPI_Allreduce(&(maxVal[0]), max_val, dofPerEnt, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&(minVal[0]), min_val, dofPerEnt, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  return M3DC1_SUCCESS;
}

#ifdef M3DC1_TRILINOS
#include <Epetra_MultiVector.h>
#include <AztecOO.h>
#include <Epetra_Version.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include "m3dc1_ls.h"

void verifyFieldEpetraVector(apf::Field* f, Epetra_MultiVector* x)
{
  double* field_data =getArrayData(f);
  assert(countComponents(f)*m3dc1_mesh::instance()->num_local_ent[0]==x->MyLength());

  for(int i=0; i<x->MyLength(); ++i)
  { 
    assert(!value_is_nan((*x)[0][i]) && !value_is_nan(field_data[i]));

    if (!(m3dc1_double_isequal((*x)[0][i], field_data[i])))
      std::cout<<"[p"<<PCU_Comm_Self()<<"] x["<<i<<"]="<<(*x)[0][i]
                <<", field_data["<<i<<"]="<<field_data[i]<<"\n";
      assert(m3dc1_double_isequal((*x)[0][i], field_data[i]));
  }
}

void copyField2EpetraVector(apf::Field* f, Epetra_MultiVector *x, bool own_only)
{
  int num_dof = countComponents(f);
  double* dof_data= new double[num_dof];
  apf::MeshEntity* e;
  apf::MeshIterator* it = m3dc1_mesh::instance()->mesh->begin(0);
  while ((e = m3dc1_mesh::instance()->mesh->iterate(it)))
  {
    if (own_only && m3dc1_mesh::instance()->mesh->getOwner(e)!=PCU_Comm_Self()) continue;
    int global_id = getNumber(get_global_numbering(), e, 0, 0);
    getComponents(f, e, 0, dof_data);
    for(int i=0; i<num_dof; ++i)
    {
      int ierr = x->ReplaceGlobalValue(global_id*num_dof+i, 0, 0, dof_data[i]);
      assert(!ierr);
    }
  }
  m3dc1_mesh::instance()->mesh->end(it);
  delete [] dof_data;

  double* in_field_data = getArrayData(f);
  for( int i=0 ; i<x->MyLength() ; ++i) 
  {
    e = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i/num_dof);    
    if (own_only && m3dc1_mesh::instance()->mesh->getOwner(e)!=PCU_Comm_Self()) continue;
    if (!(m3dc1_double_isequal((*x)[0][i], in_field_data[i])))
        std::cout<<"[p"<<PCU_Comm_Self()<<"] x["<<i<<"]="<<(*x)[0][i]
                 <<", in_field_data["<<i<<"]="<<in_field_data[i]<<"\n";
    assert(m3dc1_double_isequal((*x)[0][i], in_field_data[i]));
  }
}
#endif

int m3dc1_epetra_create(int* matrix_id, int* matrix_type, int* scalar_type, FieldID* field_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);

  if (mat)
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: matrix with id "<<*matrix_id<<" already created\n";
    return M3DC1_FAILURE; 
  }
  // check field exists
  if (!m3dc1_mesh::instance()->field_container || !m3dc1_mesh::instance()->field_container->count(*field_id))
  {
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<" failed: field with id "<<*field_id<<" doesn't exist\n";
    return M3DC1_FAILURE; 
  }
  m3dc1_ls::instance()->add_matrix(*matrix_id, new m3dc1_epetra(*matrix_id, *matrix_type, *scalar_type, *field_id));
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_delete(int* matrix_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  typedef std::map<int, m3dc1_matrix*> matrix_container_map;
  m3dc1_ls::instance()->matrix_container->erase(matrix_container_map::key_type(*matrix_id));
  delete mat;
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_insert(int* matrix_id, int* row, int* col, int* scalar_type, double* val)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  assert(*scalar_type==M3DC1_REAL);

  int err = mat->epetra_mat->ReplaceGlobalValues(*row, 1, val, col);
  if (err) {
    err =mat->epetra_mat->InsertGlobalValues(*row, 1, val, col);
    assert(err == 0);
  }
  return M3DC1_SUCCESS;
#endif
}

void print_elem (int elem_id)
{
  int ielm_dim = (m3dc1_mesh::instance()->mesh->getDimension()==2)? 2:3;
  apf::MeshEntity* e = getMdsEntity(m3dc1_mesh::instance()->mesh, ielm_dim, elem_id);

  apf::Downward downward;
  int num_node_per_element =  m3dc1_mesh::instance()->mesh->getDownward(e, 0, downward);
  switch (num_node_per_element)
  { 
    case 3: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
              <<getNumber(get_global_numbering(), downward[0], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[1], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[2], 0, 0)<<"\n";
            break;
    case 4: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
              <<getNumber(get_global_numbering(), downward[0], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[1], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[2], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[3], 0, 0)<<"\n";
            break;
    case 5: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
              <<getNumber(get_global_numbering(), downward[0], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[1], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[2], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[3], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[4], 0, 0)<<"\n";
            break;
    case 6: std::cout <<"["<<PCU_Comm_Self()<<"] elem "<<elem_id<<": nodes "
              <<getNumber(get_global_numbering(), downward[0], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[1], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[2], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[3], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[4], 0, 0)<<" "
              <<getNumber(get_global_numbering(), downward[5], 0, 0)<<"\n";
            break;
    default: break;
  }
}

#ifdef M3DC1_TRILINOS
// equivalent to Petsc::MatSetValues(*A, rsize, rows, csize, columns, &petscValues[0], ADD_VALUES);
void epetra_add_values(Epetra_CrsMatrix* mat, int rsize, int * rows, int csize, int * columns, double* values)
{
  double val[1];
  int col[1];
  assert(!mat->IndicesAreLocal() && !mat->IndicesAreContiguous());

  for (int i=0; i<rsize; i++)
  {
    for(int j=0; j<csize; j++)
    {
      col[0] = columns[j];
      val[0] = values[i*csize+j];
      int ierr = mat->SumIntoGlobalValues(rows[i], 1, val, col);
      if (ierr) 
        ierr =mat->InsertGlobalValues(rows[i], 1, val, col);
      assert(!ierr);
    }
  } // for i
}

// seol -- this does weird thing so shouldn't be used
// equivalent to Petsc::MatSetValues(*A, rsize, rows, csize, columns, &petscValues[0], ADD_VALUES);
void epetra_add_values_wrong(Epetra_CrsMatrix* mat, int rsize, int * rows, int csize, int * columns, double* values)
{
  assert(!mat->IndicesAreLocal() && !mat->IndicesAreContiguous());

  for (int i=0; i<rsize; i++)
  {
    int ierr = mat->SumIntoGlobalValues(rows[i], csize, &values[i*csize], columns);
    if (ierr) 
      ierr =mat->InsertGlobalValues(rows[i], csize, &values[i*csize], columns);
    assert(!ierr);
  } // for i
}

#endif


int m3dc1_epetra_addblock(int* matrix_id, int * ielm, int* rowVarIdx, int * columnVarIdx, double * values)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);

  if (!mat)
    return M3DC1_FAILURE;

  int field = mat->get_field_id();
  // need to change later, should get the value from field calls ...
  int dofPerVar = 6;
  char field_name[256];
  int num_values, value_type, total_num_dof; 
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  dofPerVar=total_num_dof/num_values;
  int nodes[6];
  int ent_dim=0;
  int ielm_dim = 2;
  int nodes_per_element=sizeof(nodes)/sizeof(int), nodes_per_element_get;

  if (m3dc1_mesh::instance()->mesh->getDimension()==3) ielm_dim =3;
  m3dc1_ent_getadj (&ielm_dim, ielm, &ent_dim, nodes, &nodes_per_element, &nodes_per_element_get);
  nodes_per_element=nodes_per_element_get;
  int start_global_dof_id,end_global_dof_id_plus_one;
  int start_global_dof,end_global_dof_id;
  // need to change later, should get the value from field calls ...
  int scalar_type = mat->get_scalar_type();
  assert(scalar_type==value_type);
  int numDofs = total_num_dof;
  int numVar = numDofs/dofPerVar;
  assert(*rowVarIdx<numVar && *columnVarIdx<numVar);
  int* rows = new int[dofPerVar*nodes_per_element];
  int* columns = new int[dofPerVar*nodes_per_element];

  if (mat->matrix_type==M3DC1_MULTIPLY)
  {
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      for(int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowVarIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnVarIdx)*dofPerVar+i;
      }
    }
    //FIXME: mmat->add_values(dofPerVar*nodes_per_element, rows,dofPerVar*nodes_per_element, columns, values);
    epetra_add_values(mat->epetra_mat, dofPerVar*nodes_per_element, 
                      rows,dofPerVar*nodes_per_element, columns, values);     
  }
  else //M3DC1_SOLVE
  {
    int nodeOwner[6];
    int columns_bloc[6], rows_bloc[6];
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      m3dc1_ent_getownpartid (&ent_dim, nodes+inode, nodeOwner+inode);
      m3dc1_ent_getglobaldofid (&ent_dim, nodes+inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
      rows_bloc[inode]=nodes[inode]*numVar+*rowVarIdx;
      columns_bloc[inode]=nodes[inode]*numVar+*columnVarIdx;
      for(int i=0; i<dofPerVar; i++)
      {
        rows[inode*dofPerVar+i]=start_global_dof_id+(*rowVarIdx)*dofPerVar+i;
        columns[inode*dofPerVar+i]=start_global_dof_id+(*columnVarIdx)*dofPerVar+i;
      }
    }
    int numValuesNode = dofPerVar*dofPerVar*nodes_per_element*(1+scalar_type);
    int offset=0;
    for(int inode=0; inode<nodes_per_element; inode++)
    {
      // FIXME: smat->add_values(dofPerVar, rows+dofPerVar*inode, dofPerVar*nodes_per_element, columns, values+offset);
      epetra_add_values(mat->epetra_mat, dofPerVar, rows+dofPerVar*inode, 
                       dofPerVar*nodes_per_element, columns, values+offset);
      offset += numValuesNode;
    }
  }
  delete [] rows;
  delete [] columns;

  return M3DC1_SUCCESS;
#endif
}


int m3dc1_epetra_setbc(int* matrix_id, int* row)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  int field = mat->get_field_id();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
#ifdef DEBUG
  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
#endif
  global_ordinal_type row_g = start_global_dof_id+*row%total_num_dof;
  global_ordinal_type col[1]; col[0] = row_g;
  double val[1]; val[0]=1.0; 
  // MatSetValue(*A, row, row, 1.0, ADD_VALUES);
  int err = mat->epetra_mat->SumIntoGlobalValues(row_g, 1, val, col);
  if (err) 
    err =mat->epetra_mat->InsertGlobalValues(row_g, 1, val, col);
  assert(err == 0);

  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_setlaplacebc (int * matrix_id, int *row, int * numVals, int *columns, double * values)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  std::vector <global_ordinal_type> columns_g(*numVals);
  int field = mat->get_field_id();
  int num_values, value_type, total_num_dof;
  char field_name[256];
  m3dc1_field_getinfo(&field, field_name, &num_values, &value_type, &total_num_dof);
  int inode = *row/total_num_dof;
  int ent_dim=0, start_global_dof_id, end_global_dof_id_plus_one;
  m3dc1_ent_getglobaldofid (&ent_dim, &inode, &field, &start_global_dof_id, &end_global_dof_id_plus_one);
#ifdef DEBUG
  int start_dof_id, end_dof_id_plus_one;
  m3dc1_ent_getlocaldofid (&ent_dim, &inode, &field, &start_dof_id, &end_dof_id_plus_one);
  assert(*row>=start_dof_id&&*row<end_dof_id_plus_one);
  for (int i=0; i<*numVals; i++)
    assert(columns[i]>=start_dof_id&&columns[i]<end_dof_id_plus_one);
#endif
  global_ordinal_type row_g = start_global_dof_id+*row%total_num_dof;
  for(int i=0; i<*numVals; i++)
    columns_g.at(i) = start_global_dof_id+columns[i]%total_num_dof;
//  (dynamic_cast<matrix_solve*>(mat))->set_row(row_g, *numVals, &columns_g[0], values);
  int err = mat->epetra_mat->SumIntoGlobalValues(row_g, *numVals, values, &columns_g[0]);
  if (err) 
    err =mat->epetra_mat->InsertGlobalValues(row_g, *numVals, values, &columns_g[0]);
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_print(int* matrix_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;

  // assemble matrix
  Epetra_Export exporter(/*target*/*(mat->_overlap_map),/*source*/*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  A.Export(*(mat->epetra_mat),exporter,Add);
  A.FillComplete();
  A.Print(cout);
  return M3DC1_SUCCESS;
#endif
}

#ifdef M3DC1_TRILINOS
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
#endif

int m3dc1_solver_aztec(int* matrix_id, FieldID* x_fieldid, FieldID* b_fieldid)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }
  else
    if (!PCU_Comm_Self())
      std::cout <<"[M3D-C1 INFO] "<<__func__<<": matrix "<<* matrix_id<<", field "<<* x_fieldid<<"\n";

  // assemble matrix
  Epetra_Export exporter(/*target*/*(mat->_overlap_map),/*source*/*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  A.Export(*(mat->epetra_mat),exporter,Add);
  A.FillComplete();

  // copy field to vec  
  apf::Field* b_field = (*(m3dc1_mesh::instance()->field_container))[*b_fieldid]->get_field();
  double* b_field_data = getArrayData(b_field);

  Epetra_MultiVector b_field_vec(*(mat->_overlap_map), 1);

  for (int i=0; i<b_field_vec.MyLength(); ++i)
    b_field_vec[0][i] = b_field_data[i];

  Epetra_MultiVector b(*(mat->_owned_map), 1);
  b.Export(b_field_vec,exporter,Insert);

  // vector for solution
  Epetra_MultiVector x(*(mat->_owned_map), 1);

  Epetra_LinearProblem problem(&A,&x,&b);
  AztecOO solver(problem);
  // FIXME: this crashes
  //solver.SetAztecOption(AZ_precond,AZ_Jacobi);
  //solver.SetAztecOption(AZ_output,AZ_none);

  solver.Iterate(100,1.0E-8);
  mat->num_solver_iter = solver.NumIters();

  // print residual

/* FIXME: this crashes
  double* norm_data;
  x.Norm2(norm_data);
  if (!PCU_Comm_Self()) 
    std::cout <<"[M3D-C1 INFO] "<<__func__<<" completed - #iter="<<mat->num_solver_iter
              <<", norm true resudeal="<<solver.TrueResidual()<<", Norm2 of Ax - b = " << norm_data[0] << std::endl;
*/
  // get solution
  Epetra_Import importer(*(mat->_overlap_map),*(mat->_owned_map));
  Epetra_MultiVector sol_x(*(mat->_overlap_map),1);
  sol_x.Import(x, importer, Add);

  double** s;
  sol_x.ExtractView(&s);
  apf::Field* x_field = (*(m3dc1_mesh::instance()->field_container))[*x_fieldid]->get_field();
  double* x_field_data = getArrayData(x_field);
  for (int i=0; i<sol_x.MyLength(); ++i)
    x_field_data[i] =s[0][i];

  return M3DC1_SUCCESS;
#endif
}

// solve Ax=b
#ifdef M3DC1_TRILINOS
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
#endif

int m3dc1_solver_amesos(int* matrix_id, FieldID* x_fieldid, FieldID* b_fieldid, const char* solver_name)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_SOLVE)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  Epetra_Export exporter(*(mat->_overlap_map),*(mat->_owned_map));
  Epetra_CrsMatrix A(Copy, *(mat->_owned_map), mat->nge);
  
  apf::Field* b_field = (*(m3dc1_mesh::instance()->field_container))[*b_fieldid]->get_field();
  double* b_field_data = getArrayData(b_field);

  Epetra_MultiVector b_field_vec(*(mat->_overlap_map), 1);
  // copy field to vec
  for (int i=0; i<b_field_vec.MyLength(); ++i)
    b_field_vec[0][i] = b_field_data[i];

  Epetra_MultiVector b(*(mat->_owned_map), 1 );

  A.Export(*(mat->epetra_mat),exporter,Add);
  b.Export(b_field_vec,exporter,Insert);

  A.FillComplete();

// using SuperLUDIST

  // Before we do anything, check that the solver is enabled
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  if( !Amesos2::query(solver_name) )
  {
    if (!PCU_Comm_Self()) 
      std::cout <<"[M3D-C1 ERROR] "<<__func__<<": "<< solver_name << " not enabled.  Exiting..." << std::endl;
    return M3DC1_FAILURE;	// Otherwise CTest will pick it up as failure, which it isn't really
  }
  else
    if (!PCU_Comm_Self())  *fos <<__func__<<": "<<Amesos2::version() <<" with "<<solver_name<< std::endl;
  // Constructor from Factory
  typedef Epetra_CrsMatrix MAT;
  typedef Epetra_MultiVector MV;

  Teuchos::RCP<MAT> rcp_A = Teuchos::rcp(&A);
//  A->Print(*fos);
//    *fos <<"["<<PCU_Comm_Self()<<"] GlobalNumNonZero="<<rcp_A->NumGlobalNonzeros()<<"\n";
//  int nrows = A->NumGlobalEntries();
//  if (!PCU_Comm_Self()) 
//      std::cout <<"[M3D-C1 ERROR] "<<__func__<<": nrows="<<nrows<<std::endl;

    int numVecs=1;
    Epetra_MultiVector x(*(mat->_owned_map), 1 );
    Teuchos::RCP<MV> X = Teuchos::rcp(&x);
    Teuchos::RCP<MV> B = Teuchos::rcp(&b);

    // copy field to vec
    for (int i=0; i<B->MyLength(); ++i)
      (*B)[0][i] = b_field_data[i];

    // Solve A*Xhat = B for Xhat using the Superlu solver
    Teuchos::RCP<Amesos2::Solver<MAT,MV> > solver;
    try 
    {
      solver = Amesos2::create<MAT,MV>(solver_name, rcp_A, X, B );
    }
    catch (std::invalid_argument e)
    {
      if (!PCU_Comm_Self()) *fos <<"[M3D-C1 ERROR] "<<__func__<<": "<< e.what() << std::endl;
      return M3DC1_FAILURE;
    }

    solver->symbolicFactorization();
    solver->numericFactorization();
    solver->solve();

    //solver->printTiming(*fos);
    //X.Describe(*fos, Teuchos::VERB_EXTREME);

  //  X->Print(*fos); //, Teuchos::VERB_EXTREME);
  // print residual
//  double norm_data[1];
//  x.Norm2(norm_data);
//  *fos << "["<<PCU_Comm_Self()<<"] Norm2 of Ax - b = " << norm_data[0] << std::endl;

  // get solution
  Epetra_Import importer(*(mat->_overlap_map),*(mat->_owned_map));
  Epetra_MultiVector sol_x(*(mat->_overlap_map),1);
  sol_x.Import(x, importer, Add);
 
  double** s;
  sol_x.ExtractView(&s);
  apf::Field* x_field = (*(m3dc1_mesh::instance()->field_container))[*x_fieldid]->get_field();
  double* x_field_data = getArrayData(x_field);
  for (int i=0; i<sol_x.MyLength(); ++i)
    x_field_data[i] =s[0][i];
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_solver_getnumiter(int* matrix_id, int * num_iter)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  *num_iter = mat->num_solver_iter;
  return M3DC1_SUCCESS;
#endif

}

// local matrix multiplication
// do accumulate(out_field) for global result
int m3dc1_epetra_multiply(int* matrix_id, FieldID* in_fieldid, FieldID* out_fieldid)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat || mat->matrix_type!=M3DC1_MULTIPLY)
  {
    std::cout <<"[M3D-C1 ERROR] "<<__func__<<" matrix not exists or matrix type mismatch (id"<<*matrix_id<<")\n";
    return M3DC1_FAILURE;
  }

  if (!mat->epetra_mat->Filled())
    mat->epetra_mat->FillComplete();
  assert(mat->epetra_mat->Filled());

  apf::Field* in_field = (*(m3dc1_mesh::instance()->field_container))[*in_fieldid]->get_field();
  double* in_field_data =getArrayData(in_field);

  Epetra_MultiVector x(mat->epetra_mat->RowMap(), 1);
  // copy field to vec
  for (int i=0; i<x.MyLength(); ++i)
    x[0][i] = in_field_data[i];
  Epetra_MultiVector b(mat->epetra_mat->RowMap(), 1);
  EPETRA_CHK_ERR(mat->epetra_mat->Multiply(false, x, b));
  apf::Field* out_field = (*(m3dc1_mesh::instance()->field_container))[*out_fieldid]->get_field();
  double* out_field_data =getArrayData(out_field);
  b.ExtractCopy(out_field_data, b.MyLength());
  accumulate(out_field);
  return M3DC1_SUCCESS;
#endif
}

int m3dc1_epetra_freeze(int* matrix_id)
{
#ifndef M3DC1_TRILINOS
  if (!PCU_Comm_Self()) std::cout <<"[M3D-C1 ERROR] "<<__func__<<" not supported: compile the library with \"-DENABLE_TRILINOS=ON\" in config.sh\n";
  return M3DC1_FAILURE;
#else
  m3dc1_epetra* mat = m3dc1_ls::instance()->get_matrix(*matrix_id);
  if (!mat)
    return M3DC1_FAILURE;
  mat->epetra_mat->FillComplete();
  assert(mat->epetra_mat->Filled());
  return M3DC1_SUCCESS;
#endif
}


// Ghosted Mesh Field Functions


//*******************************************************
int m3dc1_ghost_load(int* nlayers)
//*******************************************************
{ 
  if (m3dc1_model::instance()->local_planeid == 0){
    assert(*nlayers>0);
    m3dc1_ghost::instance()->mesh = apf::makeEmptyMdsMesh(m3dc1_model::instance()->model, 2, false);
    m3dc1_ghost::instance()->nlayers = *nlayers;
    m3dc1_ghost::instance()->is_ghosted = true;
    
    // Set up ghosted mesh via omega_h
    osh_t osh_mesh = osh::fromAPF(m3dc1_mesh::instance()->mesh);
    osh_ghost(&osh_mesh, m3dc1_ghost::instance()->nlayers);
    osh::toAPF(osh_mesh, m3dc1_ghost::instance()->mesh);
    osh_free(osh_mesh);
  } else {
    assert(0);
  }
  m3dc1_ghost::instance()->initialize();
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_ghost_getnument (int* /* in*/ ent_dim,
			   int* /* out */ num_ent)
//*******************************************************
{
  if (*ent_dim<0 || *ent_dim > 3)
    return M3DC1_FAILURE;
  *num_ent = m3dc1_ghost::instance()->num_local_ent[*ent_dim];
  return M3DC1_SUCCESS; 
}


//*******************************************************
int m3dc1_gfield_getnewid ( FieldID* /*out*/field_id )
//*******************************************************
{
  *field_id = fieldIdMax+1;
  return M3DC1_SUCCESS;
}

// *scalar_type is either M3DC1_REAL or M3DC1_COMPLEX
int m3dc1_gfield_create (FieldID* /*in*/ field_id,
			 const char* /* in */ field_name,
			 int* /*in*/ num_values,
			 int* /*in*/ scalar_type,
			 int* /*in*/ num_dofs_per_value)
{
  if (!m3dc1_ghost::instance()->field_container)
    m3dc1_ghost::instance()->field_container=new std::map<FieldID, m3dc1_field*>;

  // shape evaluation will be performed outside the APF
  // only need to tell APF all dofs are attached to mesh vertex
  int components = (*num_values)*(*scalar_type+1)*(*num_dofs_per_value);
  apf::Field* f = createPackedField(m3dc1_ghost::instance()->mesh, field_name, components);
  m3dc1_ghost::instance()->field_container->insert(std::map<FieldID, m3dc1_field*>::value_type(*field_id, new m3dc1_field(*field_id, f, *num_values, *scalar_type, *num_dofs_per_value)));
  apf::freeze(f); // switch dof data from tag to array

#ifdef DEBUG
  if (!PCU_Comm_Self()) 
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": ID "<<*field_id<<", #values "<<*num_values<<", #dofs "<<countComponents(f)<<", name "<<field_name<<"\n";
#endif

  if (*field_id>fieldIdMax) fieldIdMax=*field_id;
  double val[2]={0,0};
  m3dc1_gfield_assign(field_id, val, scalar_type);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_delete (FieldID* /*in*/ field_id)
//*******************************************************
{
  if (!m3dc1_ghost::instance()->field_container)
    return M3DC1_FAILURE;
  if (!m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;

  apf::Field* f = (*m3dc1_ghost::instance()->field_container)[*field_id]->get_field();
#ifdef DEBUG
  if (!PCU_Comm_Self()) 
    std::cout<<"[M3D-C1 INFO] "<<__func__<<": ID "<<*field_id<<", name "<<getName(f)<<"\n";
#endif

  destroyField(f);

  // remove f from field container
  delete (*m3dc1_ghost::instance()->field_container)[*field_id];
  m3dc1_ghost::instance()->field_container->erase(*field_id);
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_getinfo(FieldID* /*in*/ field_id, 
			 char* /* out*/ field_name,
			 int* num_values, 
			 int* scalar_type,
			 int* total_num_dof)
//*******************************************************
{
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  apf::Field* f = mf ->get_field();
  strcpy(field_name, getName(f));
  *num_values = mf -> get_num_value();
  *scalar_type = mf ->get_value_type();
  *total_num_dof = countComponents(f);
  if (*scalar_type) *total_num_dof/=2;
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_exist(FieldID* field_id,
		       int * exist)
//*******************************************************
{
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    *exist = 0;
  else
    *exist = 1;
  return M3DC1_SUCCESS;
}

//*******************************************************
void m3dc1_gfield_synchronize(apf::Field* f)
//*******************************************************
{
  apf::Mesh2* m = m3dc1_ghost::instance()->mesh;
  apf::MeshEntity* e;       

  int num_dof, n = countComponents(f);
  double* sender_data = new double[n];
  double* dof_data = new double[n]; 

  PCU_Comm_Begin();

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (!is_ent_original(m,e) || !m->isShared(e))
      continue;
    getComponents(f, e, 0, dof_data);

    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&(dof_data[0]),n*sizeof(double));
    }
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* r;
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&(sender_data[0]),n*sizeof(double));
      setComponents(f, r, 0, sender_data);
    }
  delete [] dof_data;
  delete [] sender_data;
}

//*******************************************************
int m3dc1_gfield_sync (FieldID* /* in */ field_id)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_gfield_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_gfield_synchronize((*m3dc1_ghost::instance()->field_container)[*field_id]->get_field());
#ifdef DEBUG
  m3dc1_gfield_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

// send non-owned copies' dof to owner copy and add them up
//*******************************************************
void m3dc1_gfield_accumulate(apf::Field* f)
//*******************************************************
{
  apf::Mesh2* m = m3dc1_ghost::instance()->mesh;
  apf::MeshEntity* e;       

  int num_dof, own_partid, n = countComponents(f);
  double* dof_data = new double[n];
  double* sender_data = new double[n];
  apf::MeshEntity* own_e;
  apf::MeshEntity* r;
  std::map<apf::MeshEntity*, std::vector<double> > save_map;

  PCU_Comm_Begin();

  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    own_partid=get_ent_ownpartid(m, e);
    if (own_partid==PCU_Comm_Self()) continue;

    own_e = get_ent_owncopy(m, e);

    getComponents(f, e, 0, &(dof_data[0]));
      
    PCU_COMM_PACK(own_partid, own_e);
    PCU_Comm_Pack(own_partid,&(dof_data[0]),n*sizeof(double));
  }
  m->end(it);

  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&(sender_data[0]),n*sizeof(double));
      for (int i = 0; i < n; ++i)
save_map[r].push_back(sender_data[i]);      
    }

  for (std::map<apf::MeshEntity*, std::vector<double> >::iterator mit=save_map.begin(); mit!=save_map.end(); ++mit)
  {
    e = mit->first;
    getComponents(f, e, 0, dof_data);
    int num_data = mit->second.size()/n;
    for (int i=0; i<num_data;++i)
    {
      for (int j=0; j<n; ++j)
dof_data[j] += mit->second[i*n+j];
    }
    setComponents(f, e, 0, dof_data);
  } 
  delete [] dof_data;
  delete [] sender_data;
}


//*******************************************************
int m3dc1_gfield_sum (FieldID* /* in */ field_id)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_gfield_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_gfield_accumulate((*m3dc1_ghost::instance()->field_container)[*field_id]->get_field());
  m3dc1_gfield_synchronize((*m3dc1_ghost::instance()->field_container)[*field_id]->get_field());
#ifdef DEBUG
  m3dc1_gfield_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_gfield_sumsq (FieldID* /* in */ field_id,
			double* /* out */ sum)
//*******************************************************
{
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  
  apf::Field* f = (*(m3dc1_ghost::instance()->field_container))[*field_id]->get_field();

#ifdef DEBUG
  int isnan;
  m3dc1_gfield_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  *sum=0.;
  int num_dof = countComponents(f);

  double* dof_data= new double[num_dof];
  apf::MeshEntity* e;
  apf::MeshIterator* it = m3dc1_ghost::instance()->mesh->begin(0);
  while ((e = m3dc1_ghost::instance()->mesh->iterate(it)))
  {
    if (!is_ent_original(m3dc1_ghost::instance()->mesh,e)) continue;
    getComponents(f, e, 0, dof_data);
    for(int i=0; i<num_dof; ++i)
      *sum+=dof_data[i]*dof_data[i];
  }
  m3dc1_ghost::instance()->mesh->end(it);
  delete [] dof_data;
  return M3DC1_SUCCESS;
}


/** field dof functions */
//*******************************************************
int m3dc1_gfield_getlocaldofid (FieldID* field_id, 
				int* /* out */ start_dof_id,
				int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  *start_dof_id=0;
  *end_dof_id_plus_one=num_dof*m3dc1_ghost::instance()->num_local_ent[0];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_getowndofid (FieldID* field_id, 
			      int* /* out */ start_dof_id,
			      int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();

  int num_own_ent = m3dc1_ghost::instance()->num_own_ent[0];
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  
  int start_id = num_own_ent;
  PCU_Exscan_Ints(&start_id,1);

  *start_dof_id=start_id*num_dof;
  *end_dof_id_plus_one=*start_dof_id+num_own_ent*num_dof;
  return M3DC1_SUCCESS;
}
 
//******************************************************* 
int m3dc1_gfield_getglobaldofid (FieldID* field_id, 
				 int* /* out */ start_dof_id,
				 int* /* out */ end_dof_id_plus_one)
//*******************************************************
{
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
  m3dc1_field * mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  apf::Field* f = mf->get_field();
  int num_dof = (mf->get_value_type())?countComponents(f)/2:countComponents(f);
  assert(mf->get_num_value()*mf->get_dof_per_value()==num_dof);  

  *start_dof_id=0;
  *end_dof_id_plus_one=*start_dof_id+num_dof*m3dc1_ghost::instance()->num_global_ent[0];
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_getnumlocaldof (FieldID* field_id,
				 int* /* out */ num_local_dof)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  m3dc1_field * mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  *num_local_dof = (m3dc1_ghost::instance()->num_local_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_getnumowndof (FieldID* field_id,
			       int* /* out */ num_own_dof)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  m3dc1_field * mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  *num_own_dof = (m3dc1_ghost::instance()->num_own_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_getnumglobaldof (FieldID* field_id,
				  int* /* out */ num_global_dof)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  m3dc1_field * mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  *num_global_dof = (m3dc1_ghost::instance()->num_global_ent[0])*mf->get_num_value()*mf->get_dof_per_value();
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_getdataptr (FieldID* field_id,
			     double** pts)
//*******************************************************
{
#ifdef DEBUG
  if (!m3dc1_ghost::instance()->field_container ||
      !m3dc1_ghost::instance()->field_container->count(*field_id))
    return M3DC1_FAILURE;
#endif
  apf::Field* f = (*(m3dc1_ghost::instance()->field_container))[*field_id]->get_field();
  if (!isFrozen(f)) freeze(f);
  *pts=getArrayData(f);
  return M3DC1_SUCCESS;
}

// add field2 to field1
//*******************************************************
int m3dc1_gfield_add(FieldID* /*inout*/ field_id1,
		     FieldID* /*in*/ field_id2)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_gfield_isnan(field_id1, &isnan);
  assert(isnan==0);
  m3dc1_gfield_isnan(field_id2, &isnan);
  assert(isnan==0);
#endif
  m3dc1_field * mf1 = (*(m3dc1_ghost::instance()->field_container))[*field_id1];
  m3dc1_field * mf2 = (*(m3dc1_ghost::instance()->field_container))[*field_id2];
  int dofPerEnt1 = mf1->get_num_value()*mf1->get_dof_per_value(); 
  int dofPerEnt2 = mf2->get_num_value()*mf2->get_dof_per_value();
  assert(mf1->get_value_type()==mf2->get_value_type());
  std::vector<double> dofs1(dofPerEnt1*(1+mf1->get_value_type())), dofs2(dofPerEnt2*(1+mf2->get_value_type()));
  int dofMin = std::min(dofPerEnt1,dofPerEnt2);
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_ghost_getnument (&vertex_type, &num_vtx);
  int dofPerEntDummy[2];
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id1, dofPerEntDummy, &dofs1[0]);
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id2, dofPerEntDummy+1, &dofs2[0]);
    for(int i=0; i<dofMin*(1+mf1->get_value_type()); i++)
      dofs1.at(i)+=dofs2.at(i);
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id1, &dofPerEnt1, &dofs1[0]);
  } 
#ifdef DEBUG
  m3dc1_gfield_isnan(field_id1, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_mult(FieldID* /*inout*/ field_id,
		      double* fac,
		      int * scalar_type)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_gfield_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_ghost_getnument (&vertex_type, &num_vtx);
  double dofs[FIXSIZEBUFF], dofsNew[FIXSIZEBUFF];
  m3dc1_field* mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  int value_type = mf->get_value_type();
  assert(dofPerEnt<=sizeof(dofs)/sizeof(double)*(1+value_type));
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, dofs);
    if(*scalar_type==0)
    {
      for(int i=0; i<dofPerEnt*(1+value_type); i++)
dofsNew[i]=*fac*dofs[i];
    }
    else
    {
      for(int i=0; i<dofPerEnt; i++)
      {
dofsNew[2*i]=fac[0]*dofs[2*i]-fac[1]*dofs[2*i+1];
dofsNew[2*i+1]=fac[0]*dofs[2*i+1]+fac[1]*dofs[2*i];
      }
    }
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofsNew[0]);
  }
#ifdef DEBUG
  m3dc1_gfield_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_assign(FieldID* /*inout*/ field_id,
			double* fac,
			int * scalar_type)
//*******************************************************
{
  m3dc1_field* mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  if (dofPerEnt==0) return M3DC1_FAILURE;

  int num_vtx=0;
  int vertex_type=0;
  m3dc1_ghost_getnument (&vertex_type, &num_vtx);
  std::vector<double> dofs(dofPerEnt*(1+mf->get_value_type()), fac[0]);
  if(*scalar_type)
    for(int i=0; i<dofPerEnt; i++)
      dofs.at(2*i+1)=fac[1];
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofs[0]);
  }
#ifdef DEBUG
  int isnan;
  m3dc1_gfield_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_copy(FieldID* /* out */ field_id1,
		      FieldID* /* in */ field_id2)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_gfield_isnan(field_id1, &isnan);
  assert(isnan==0);
#endif
  m3dc1_field * mf1 = (*(m3dc1_ghost::instance()->field_container))[*field_id1];
  m3dc1_field * mf2 = (*(m3dc1_ghost::instance()->field_container))[*field_id2];
  apf::Field* f1 =  mf1->get_field();
  apf::Field* f2 =  mf2->get_field();
  int dofPerEnt1 = mf1->get_num_value()*mf1->get_dof_per_value();
  int dofPerEnt2 = mf2->get_num_value()*mf2->get_dof_per_value();
  assert(mf1->get_value_type()==mf2->get_value_type());
  std::vector<double> dofs1(dofPerEnt1*(1+mf1->get_value_type())), dofs2(dofPerEnt2*(1+mf2->get_value_type()));
  int dofMin = std::min(dofPerEnt1,dofPerEnt2);
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_ghost_getnument (&vertex_type, &num_vtx);
  int dofPerEntDummy[2];
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id1, dofPerEntDummy, &dofs1[0]);
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id2, dofPerEntDummy+1, &dofs2[0]);
    for(int i=0; i<dofMin*(1+mf1->get_value_type()); i++)
      dofs1.at(i)=dofs2.at(i);
    m3dc1_ent_setdofdata (&vertex_type, &inode, field_id1, &dofPerEnt1, &dofs1[0]);
  }
#ifdef DEBUG
  m3dc1_gfield_isnan(field_id2, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_retrieve (FieldID* /* in */ field_id,
			   double * /*out*/ data,
			   int * /* in */size)
//*******************************************************
{
#ifdef DEBUG
  int isnan;
  m3dc1_gfield_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  //m3dc1_gfield_print(field_id);
  int num_local_dof, num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_gfield_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
  m3dc1_gfield_getnumlocaldof (field_id, &num_local_dof);

  double* pts=NULL;
  m3dc1_gfield_getdataptr (field_id, &pts);
  memcpy(data, pts, *size*(1+value_type)*sizeof(double));
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_set (FieldID* /* in */ field_id,
		      double * /*in*/ data,
		      int * /* in */size)
//*******************************************************
{
  int num_local_dof, num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_gfield_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
  m3dc1_gfield_getnumlocaldof (field_id, &num_local_dof);

  double* pts=NULL;
  m3dc1_gfield_getdataptr (field_id, &pts);
  memcpy(pts, data, *size*(1+value_type)*sizeof(double));
#ifdef DEBUG
  int isnan;
  m3dc1_gfield_isnan(field_id, &isnan);
  assert(isnan==0);
#endif
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_insert(FieldID* /* in */ field_id,
			int /* in */ * local_dof, 
			int * /* in */ size,
			double* /* in */ values,
			int * type,
			int * op)
//*******************************************************
{
  int num_local_dof, num_values, value_type, total_num_dof;
  char field_name[FIXSIZEBUFF];
  m3dc1_gfield_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);
#ifdef DEBUG
  m3dc1_gfield_getnumlocaldof (field_id, &num_local_dof);
  assert(*local_dof<num_local_dof);
  if(!value_type) assert(!(*type)); // can not insert complex value to real vector
  for(int i=0; i<*size*(1+(*type)); i++)
  {
#ifdef REPLACENANWITHZERO
    if(values[i]!=values[i]) values[i]=0;
#else
    assert(values[i]==values[i]);
#endif
  }
#endif
  std::vector<double> values_convert(*size*(1+value_type),0);
  if(!(*type)&&value_type) // real into complex
  {
    for(int i=0; i<*size; i++)
    {
      values_convert.at(2*i)=values[i];
    }
  }
  else
  {
    for(int i=0; i<*size*(1+value_type); i++)
    {
      values_convert.at(i)=values[i];
    }
  }
  double * dataptr;
  int ibegin=*local_dof*(1+value_type);
  m3dc1_gfield_getdataptr(field_id, &dataptr);
  if(*op==0) // set value
  {
   for(int i=0; i<*size*(1+value_type); i++)
     dataptr[ibegin+i]=values_convert.at(i);
  }
  else
  {
    for(int i=0; i<*size*(1+value_type); i++)
      dataptr[ibegin+i]+=values_convert[i];
  }
  return M3DC1_SUCCESS;
}


//*******************************************************
int m3dc1_gfield_isnan(FieldID* /* in */ field_id,
		       int * isnan)
//*******************************************************
{
  *isnan=0;
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_ghost_getnument (&vertex_type, &num_vtx);
  int dofPerEnt;
  double dofs[FIXSIZEBUFF];
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, dofs);
    for(int i=0; i<dofPerEnt; i++)
      if(value_is_nan(dofs[i])) 
        *isnan=1;
  }
  return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_write(FieldID* field_id,
		       const char* file_name)
//*******************************************************
{ 
  m3dc1_field* mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  if (dofPerEnt==0) return M3DC1_FAILURE;
  FILE * fp =fopen(file_name, "w");
  int num_vtx=0;
  int vertex_type=0;
  m3dc1_ghost_getnument (&vertex_type, &num_vtx);
  std::vector<double> dofs(dofPerEnt*(1+mf->get_value_type()),0);
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofs[0]);
    for(int i=0; i<dofPerEnt; i++)
    {
      for(int j=0; j<1+mf->get_value_type(); j++)
        fprintf(fp, "%lf ", dofs[i*(1+mf->get_value_type())+j]);
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
}


//*******************************************************
int m3dc1_gfield_compare(FieldID* field_id_1,
			 FieldID* field_id_2)
//*******************************************************
{
  apf::Field* f_1 = (*(m3dc1_ghost::instance()->field_container))[*field_id_1]->get_field();
  double* field_data_1 =getArrayData(f_1);

  apf::Field* f_2 = (*(m3dc1_ghost::instance()->field_container))[*field_id_2]->get_field();
  double* field_data_2 =getArrayData(f_2);

  int num_dof_1=countComponents(f_1);
  int num_dof_2=countComponents(f_2);
  if (num_dof_1!=num_dof_2) 
  {
    if (!PCU_Comm_Self()) 
      cout<<"[M3DC1 INFO] "<<__func__<<": #dof mismatch "<<getName(f_1)
          <<"- "<<num_dof_1<<", "<<getName(f_2)<<"- "<<num_dof_2<<"\n";
    return M3DC1_FAILURE;
  }
  int ierr = M3DC1_SUCCESS;
  for (int i=0; i<num_dof_1*m3dc1_ghost::instance()->num_local_ent[0]; ++i)
  {  if (!m3dc1_double_isequal(field_data_1[i], field_data_2[i])) 
    {
     cout<<"[M3DC1 ERROR] "<<__func__<<": "<<getName(f_1)<<"["<<i<<"]="<<field_data_1[i]
          <<", "<<getName(f_2)<<"["<<i<<"]="<<field_data_2[i]<<"\n";
      ierr=M3DC1_FAILURE;
      break;
    }
  }
  int global_ierr;
  MPI_Allreduce(&ierr, &global_ierr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
  if (global_ierr==M3DC1_FAILURE)
  {
    if (!PCU_Comm_Self())
      cout<<"[M3DC1 INFO] "<<__func__<<": dof value mismatch of fields "<<getName(f_1)
          <<" and "<<getName(f_2)<<"\n";
    
    return M3DC1_FAILURE;
  }
  else
    return M3DC1_SUCCESS;
}

//*******************************************************
int m3dc1_gfield_print(FieldID* field_id)
//*******************************************************
{ 
  apf::Field* f = (*(m3dc1_ghost::instance()->field_container))[*field_id]->get_field();
  double* field_data =getArrayData(f);
  
  apf::MeshEntity* e;
  int global_ent_id, num_dof=countComponents(f);
  double* dof_data = new double[num_dof];

  switch (num_dof)
  {
    case 1: {
              int i=0;
              apf::MeshIterator* it = m3dc1_ghost::instance()->mesh->begin(0);
              while ((e = m3dc1_ghost::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_ghost::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<"]\n";
               }
               m3dc1_ghost::instance()->mesh->end(it);
               break;
            }
    case 2: {
              int i=0;
              apf::MeshIterator* it = m3dc1_ghost::instance()->mesh->begin(0);
              while ((e = m3dc1_ghost::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_ghost::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<"]\n";
               }
               m3dc1_ghost::instance()->mesh->end(it);
               break;
            }
    case 3: {
              int i=0;
              apf::MeshIterator* it = m3dc1_ghost::instance()->mesh->begin(0);
              while ((e = m3dc1_ghost::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_ghost::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<"]\n";
               }
               m3dc1_ghost::instance()->mesh->end(it);
               break;
            }
    case 4: {
              int i=0;
              apf::MeshIterator* it = m3dc1_ghost::instance()->mesh->begin(0);
              while ((e = m3dc1_ghost::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_ghost::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<"]\n";
               }
               m3dc1_ghost::instance()->mesh->end(it);
               break;
            }
    case 6: {
              int i=0;
              apf::MeshIterator* it = m3dc1_ghost::instance()->mesh->begin(0);
              while ((e = m3dc1_ghost::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_ghost::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<", "<<dof_data[4]
               <<", "<<dof_data[5]
               <<"]\n";
               }
               m3dc1_ghost::instance()->mesh->end(it);
               break;
          }
    case 8: {
              int i=0;
              apf::MeshIterator* it = m3dc1_ghost::instance()->mesh->begin(0);
              while ((e = m3dc1_ghost::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_ghost::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<", "<<dof_data[4]
               <<", "<<dof_data[5]
               <<", "<<dof_data[6]
               <<", "<<dof_data[7]
               <<"]\n";
               }
               m3dc1_ghost::instance()->mesh->end(it);
               break;
            }
    case 12: {
              int i=0;
              apf::MeshIterator* it = m3dc1_ghost::instance()->mesh->begin(0);
              while ((e = m3dc1_ghost::instance()->mesh->iterate(it)))
  //            for (int i=0; i<m3dc1_mesh::instance()->num_local_ent[0]; ++i)
              { 
                assert(e ==getMdsEntity(m3dc1_ghost::instance()->mesh, 0, i++));
               // if (m3dc1_mesh::instance()->mesh->getOwner(e)!=PCU_Comm_Self()) continue;
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<", "<<dof_data[4]
               <<", "<<dof_data[5]
               <<", "<<dof_data[6]
               <<", "<<dof_data[7]
               <<", "<<dof_data[8]
               <<", "<<dof_data[9]
               <<", "<<dof_data[10]
               <<", "<<dof_data[11]
               <<"]\n";
               }
               m3dc1_ghost::instance()->mesh->end(it);
               break;
             }
    case 18: {
              int i=0;
              apf::MeshIterator* it = m3dc1_ghost::instance()->mesh->begin(0);
              while ((e = m3dc1_ghost::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_ghost::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<", "<<dof_data[4]
               <<", "<<dof_data[5]
               <<", "<<dof_data[6]
               <<", "<<dof_data[7]
               <<", "<<dof_data[8]
               <<", "<<dof_data[9]
               <<", "<<dof_data[10]
               <<", "<<dof_data[11]
               <<", "<<dof_data[12]
               <<", "<<dof_data[13]
               <<", "<<dof_data[14]
               <<", "<<dof_data[15]
               <<", "<<dof_data[16]
               <<", "<<dof_data[17]
               <<"]\n";
               }
               m3dc1_ghost::instance()->mesh->end(it);
               break;
             }
    case 24: {
              int i=0;
              apf::MeshIterator* it = m3dc1_ghost::instance()->mesh->begin(0);
              while ((e = m3dc1_ghost::instance()->mesh->iterate(it)))
              { 
                assert(e ==getMdsEntity(m3dc1_ghost::instance()->mesh, 0, i++));
                global_ent_id = getNumber(get_global_numbering(), e, 0, 0);
                getComponents(f, e, 0, dof_data);
                std::cout<<"[p"<<PCU_Comm_Self()<<"] field "<<getName(f)
               <<"/ent "<<global_ent_id
               <<": ["<<dof_data[0]
               <<", "<<dof_data[1]
               <<", "<<dof_data[2]
               <<", "<<dof_data[3]
               <<", "<<dof_data[4]
               <<", "<<dof_data[5]
               <<", "<<dof_data[6]
               <<", "<<dof_data[7]
               <<", "<<dof_data[8]
               <<", "<<dof_data[9]
               <<", "<<dof_data[10]
               <<", "<<dof_data[11]
               <<", "<<dof_data[12]
               <<", "<<dof_data[13]
               <<", "<<dof_data[14]
               <<", "<<dof_data[15]
               <<", "<<dof_data[16]
               <<", "<<dof_data[17]
               <<", "<<dof_data[18]
               <<", "<<dof_data[19]
               <<", "<<dof_data[20]
               <<", "<<dof_data[21]
               <<", "<<dof_data[22]
               <<", "<<dof_data[23]
               <<"]\n";
               }
               m3dc1_ghost::instance()->mesh->end(it);
               break;
            }
      default: if (!PCU_Comm_Self()) std::cout<<__func__<<" failed for field "<<getName(f)<<": does support "<<num_dof<<" dofs\n";
              return M3DC1_FAILURE;
  } // switch
  return M3DC1_SUCCESS;
}

#ifndef M3DC1_MESHGEN
//*******************************************************
int m3dc1_gfield_sum_plane (FieldID* /* in */ field_id)
//*******************************************************
{
  MPI_Comm icomm= m3dc1_model::instance()->getMPICommPlane();
  int num_vtx,num_dof=0, vertex_type=0;
  m3dc1_gfield_getnumlocaldof(field_id, &num_dof);
  char field_name[256];
  int num_values, value_type, total_num_dof;
  m3dc1_gfield_getinfo(field_id, field_name, &num_values, &value_type, &total_num_dof);

  m3dc1_ghost_getnument (&vertex_type, &num_vtx);

  int data_size=total_num_dof*num_vtx*(1+value_type);
  assert(total_num_dof*num_vtx==num_dof);
  double * thevec = NULL;
  m3dc1_gfield_getdataptr(field_id, &thevec);
  double * sendbuf = new double [data_size];
  m3dc1_gfield_retrieve (field_id, sendbuf, &num_dof);
  MPI_Allreduce (sendbuf,thevec,data_size,MPI_DOUBLE,MPI_SUM,icomm) ;
  m3dc1_gfield_sync(field_id);
  delete []sendbuf;
}

int adapt_by_gfield (int * fieldId,
		     double* psi0,
		     double * psil)
{
  FILE *fp = fopen("sizefieldParam", "r");
  if(!fp)
  {
    std::cout<<" file sizefieldParam not found "<<std::endl;
    throw 1;
  }
  double param[13];
  set<int> field_keep;
  field_keep.insert(*fieldId);
  apf::Mesh2* mesh = m3dc1_ghost::instance()->mesh;
  /*
  for (int i=0; i < mesh->countFields(); ++i)
  {
    Field* f = mesh->getField(i);
    std::cout<<"find field "<<getName(f)<<endl;
  }
  for (int i=0; i < mesh->countNumberings(); ++i)
  {
    Numbering* n = mesh->getNumbering(i);
    std::cout<<"find numbering "<<getName(n)<<endl;
  }*/
  //m3dc1_mesh ::instance ()->clean(field_keep);
  //apf::writeVtkFiles("before",m3dc1_mesh::instance()->mesh);
  //mesh->writeNative("before.smb");
  // read the size field parameters
  for(int i=0; i<13; i++)
    fscanf(fp, "%lf ", &param[i]);
  fclose(fp);
  apf::Field* psiField = (*(m3dc1_ghost::instance()->field_container))[*fieldId]->get_field();

  int node_glb_order=NODE_GLB_ORDER;
  m3dc1_gfield_delete (&node_glb_order);

  // delete all the matrix
  std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
  while (m3dc1_solver::instance()-> matrix_container->size())
  {
    std::map<int, m3dc1_matrix*> :: iterator mat_it = m3dc1_solver::instance()-> matrix_container->begin();
    delete mat_it->second;
    m3dc1_solver::instance()->matrix_container->erase(mat_it);
  }

  int valueType = (*(m3dc1_ghost::instance()->field_container))[*fieldId]->get_value_type();
  SizeFieldPsi sf (psiField, *psi0, *psil, param, valueType);
  double mmax[2], mmin[2];
  m3dc1_model_getmaxcoord(mmax,mmax+1);
  m3dc1_model_getmincoord(mmin,mmin+1);
  //Vortex sfv(mesh, center, mmax[0]-mmin[0]);

  ReducedQuinticImplicit shape;
  vector<apf::Field*> fields;
  std::map<FieldID, m3dc1_field*> :: iterator it=m3dc1_ghost::instance()->field_container->begin();
  while(it!=m3dc1_ghost::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    assert(valueType==complexType);
    if(complexType) group_complex_dof(field, 1);
    if(isFrozen(field)) unfreeze(field);
    if(!PCU_Comm_Self()) std::cout<<"Solution transfer: add field "<<apf::getName(field)<<std::endl;
    fields.push_back(field);
    it++;
  }
  while(mesh->countNumberings())
  {
    apf::Numbering* n = mesh->getNumbering(0);
    if(!PCU_Comm_Self()) std::cout<<"destroying numbering "<<getName(n)<<endl;
    apf::destroyNumbering(n);
  }
  ReducedQuinticTransfer slnTrans(mesh,fields, &shape);
  ma::Input* in = ma::configure(mesh,&sf,&slnTrans);
  //ma::Input* in = ma::configure(mesh,&sfv);
  in->maximumIterations = 9;
  //char filename[256];
  //sprintf(filename,"before%d",adapt_time);
  //apf::writeVtkFiles(filename,mesh);
  //set<int> field_keep;
  //field_keep.insert(*fieldId);
  //m3dc1_mesh ::instance ()->clean(field_keep);
  in->shouldRunPostZoltan = true;
  ma::adapt(in);
  reorderMdsMesh(mesh);
  //mesh->verify();
  
  /*while(mesh->countFields())
  {
    Field* f = mesh->getField(0);
    std::cout<<"destroying field "<<getName(f)<<endl;
    destroyField(f);
  }*/
  /*
  while(mesh->countNumberings())
  {
    Numbering* n = mesh->getNumbering(0);
    std::cout<<"destroying numbering after adapt "<<getName(n)<<endl;
    destroyNumbering(n);
  }*/
  //sprintf(filename,"adapted%d",adapt_time++);
  //apf::writeVtkFiles(filename,mesh);
  //mesh->writeNative("adapted.smb");

  m3dc1_mesh::instance()->initialize();
  it=m3dc1_ghost::instance()->field_container->begin();
  while(it!=m3dc1_ghost::instance()->field_container->end())
  {
    apf::Field* field = it->second->get_field();
    int complexType = it->second->get_value_type();
    if(complexType) group_complex_dof(field, 0);
    if(!isFrozen(field)) freeze(field);
#ifdef DEBUG
    int isnan;
    int fieldId= it->first; 
    m3dc1_gfield_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    m3dc1_gfield_synchronize(field);
    get_global_numbering();
#ifdef DEBUG
    m3dc1_gfield_isnan(&fieldId, &isnan);
    assert(isnan==0);
#endif
    it++;
  }
}
#endif

int m3dc1_gfield_printcompnorm(FieldID* /* in */ field_id,
			       const char* info)
{
  double* pts=NULL;
  m3dc1_gfield_getdataptr (field_id, &pts);
  int num_dof;
  m3dc1_gfield_getnumlocaldof(field_id, &num_dof);
  m3dc1_field * mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];

  apf::Field* f = mf ->get_field();
  int dof_per_node = countComponents(f);
  assert(dof_per_node == num_dof * 6);
  int num_comp=dof_per_node/C1TRIDOFNODE;
  vector<double> norms(dof_per_node/C1TRIDOFNODE);
  int j=0;
  for(int i=0; i<num_dof/C1TRIDOFNODE; i++)
  {
     for(int k=0; k<6; k++)
        norms.at(j)+=pts[i*C1TRIDOFNODE+k]*pts[i*C1TRIDOFNODE+k];
     j++;
     j%=num_comp;
  }
  vector<double> buff=norms;
  MPI_Allreduce(&buff[0],&norms[0], num_comp, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  int psize;
  MPI_Comm_size(MPI_COMM_WORLD,&psize);
  if(PCU_Comm_Self() == psize-1)
  {
    std::cout<< "norm of vec "<<info;
    for(int i=0; i<num_comp; i++)
      std::cout<<" "<<std::sqrt(norms[i]);
    std::cout<<std::endl;
  }
}

int m3dc1_gfield_max (FieldID* field_id,
		      double * max_val,
		      double * min_val)
{
  m3dc1_field* mf = (*(m3dc1_ghost::instance()->field_container))[*field_id];
  int dofPerEnt = mf->get_num_value()*mf->get_dof_per_value();
  if(mf->get_value_type()) dofPerEnt *= 2;
  if (dofPerEnt==0) return M3DC1_FAILURE;

  int num_vtx=0;
  int vertex_type=0;
  m3dc1_ghost_getnument (&vertex_type, &num_vtx);
  std::vector<double> maxVal(dofPerEnt, -1e30), minVal(dofPerEnt,1e30), dofs(dofPerEnt);
  for(int inode=0; inode<num_vtx; inode++)
  {
    m3dc1_ent_getdofdata (&vertex_type, &inode, field_id, &dofPerEnt, &dofs[0]);
    for(int i=0; i<dofPerEnt; i++)
    {
      if(maxVal[i]<dofs[i]) maxVal[i]=dofs[i];
      if(minVal[i]>dofs[i]) minVal[i]=dofs[i];
    }
  }
  MPI_Allreduce(&(maxVal[0]), max_val, dofPerEnt, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&(minVal[0]), min_val, dofPerEnt, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  return M3DC1_SUCCESS;
}



