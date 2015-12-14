#ifndef NOM3DC1
#define NOM3DC1
#endif
#include "m3dc1_scorec.h"
#include <iostream>
#include <assert.h>
#include "apf.h"
#include "petscksp.h"
#include "PCU.h"
#include "m3dc1_mesh.h" // debugging purpose
#include "apfMesh.h" // debugging purpose
#include <apfMDS.h>
using namespace std;
static char help[] = "testing solver functions; \n first do mat-vec product A*b=c; solve Ax=c; compare x and b\n\n";

bool AlmostEqualDoubles(double A, double B,
            double maxDiff, double maxRelDiff);

#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <parma.h>
#if defined(__linux__)
#include <malloc.h>
#else
#include <cstdlib>
#endif

#ifdef __bgq__
#include <spi/include/kernel/memory.h>

static double get_peak()
{
  uint64_t heap;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  return heap;
}

#elif defined (__linux__)

static double get_peak()
{
  return mallinfo().arena;
}

#else

static double get_peak()
{
  if(!PCU_Comm_Self())
    printf("%s:%d: OS Not supported\n", __FILE__, __LINE__);
  return(-1.0);
}

#endif

static void print_stats(const char* name, double value)
{
  double min, max, avg;
  min = value;
  PCU_Min_Doubles(&min, 1);
  max = value;
  PCU_Max_Doubles(&max, 1);
  avg = value;
  PCU_Add_Doubles(&avg, 1);
  avg /= PCU_Comm_Peers();
  double imb = max / avg;
  if (!PCU_Comm_Self())
    printf("%s: min %f max %f avg %f imb %f\n", name, min, max, avg, imb);
}

#if defined(__linux__)

static double get_chunks()
{
  struct mallinfo m = mallinfo();
  return m.uordblks + m.hblkhd;
}

#else
static double get_chunks()
{
  if(!PCU_Comm_Self())
    printf("%s:%d: OS Not supported\n", __FILE__, __LINE__);
  return(-1.0);
}
#endif


static void list_tags(apf::Mesh* m)
{
  if (PCU_Comm_Self())
    return;
  apf::DynamicArray<apf::MeshTag*> tags;
  m->getTags(tags);
  for (size_t i = 0; i < tags.getSize(); ++i)
    printf("tag: \"%s\"\n",m->getTagName(tags[i]));
}

int main( int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  m3dc1_scorec_init();
  int op=0;
  m3dc1_matrix_setassembleoption(&op);
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  PetscLogDouble mem;
  if (argc<4 & !PCU_Comm_Self())
  {
    cout<<"Usage: ./main  model mesh #planes real(0)/complex(1) "<<endl;
    return M3DC1_FAILURE;
  }
  double t1, t2, t3, t4, t5;
  double tf1,tf2;

  int num_plane=1;
  if (argc>3)
  {
    num_plane = atoi(argv[3]);
    if (num_plane>1 && PCU_Comm_Peers()%num_plane==0)
      m3dc1_model_setnumplane (&num_plane);
  }

  m3dc1_model_load(argv[1]);
  m3dc1_model_print();

  int scalar_type=0;
  if(argc>4) scalar_type=atoi(argv[4]);
#ifndef PETSC_USE_COMPLEX
  if(scalar_type!=0)
  {
    if(!PCU_Comm_Self()) cout<<"PETSc is not configured with --with-scalar-type=complex"<<endl;
    return 1;
  } 
;
#else
  if(!PCU_Comm_Self()) cout<<"PETSC_USE_COMPLEX defined"<<endl;
#endif

  m3dc1_mesh_load(argv[2]);
  //int three=3;
  //m3dc1_mesh_write("geoId", &three);

  printStats(m3dc1_mesh::instance()->mesh);
  if (num_plane>1)
  {
    int zero=0;
    m3dc1_mesh_build3d(&zero, &zero, &zero);
  }
 
  if (atoi(argv[3])==0) // model/mesh testing
  {
    apf::MeshTag* tag_loaded = m3dc1_mesh::instance()->mesh->findTag("norm_curv");
    if (tag_loaded)
    {
    double curv;
    double norm[3];
    double norm_curv[3];
    apf::MeshEntity* ment;
    for (int i=0; i<m3dc1_mesh::instance()->mesh->count(0); ++i)
    {
      m3dc1_node_getnormvec(&i, &norm[0]);
      m3dc1_node_getcurv (&i, &curv);
      ment = getMdsEntity(m3dc1_mesh::instance()->mesh, 0, i);
      m3dc1_mesh::instance()->mesh->getDoubleTag(ment, tag_loaded, &norm_curv[0]);
      assert(m3dc1_double_isequal(norm[0], norm_curv[0]));
      assert(m3dc1_double_isequal(norm[1], norm_curv[1]));
      assert(m3dc1_double_isequal(norm[2], 0.0));
      assert(m3dc1_double_isequal(curv, norm_curv[2]));
      } 
    }
  apf::writeVtkFiles("mesh",m3dc1_mesh::instance()->mesh);
  m3dc1_mesh::instance()->mesh->writeNative("mesh.smb");
  print_stats("kernel heap", get_peak());
  print_stats("malloc used", get_chunks());
  print_stats("elements", m3dc1_mesh::instance()->mesh->count(m3dc1_mesh::instance()->mesh->getDimension()));
  print_stats("vertices", m3dc1_mesh::instance()->mesh->count(0));
  Parma_PrintPtnStats(m3dc1_mesh::instance()->mesh, "");
  list_tags(m3dc1_mesh::instance()->mesh);
    PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return M3DC1_SUCCESS;
  }
  // set/get field dof values
  int num_vertex, num_own_vertex, vertex_dim=0;
  int num_elem, elem_dim=2;
  if(num_plane>1) elem_dim=3; // use wedge in 3D
  m3dc1_mesh_getnument (&vertex_dim, &num_vertex);
  m3dc1_mesh_getnumownent (&vertex_dim, &num_own_vertex);
  m3dc1_mesh_getnument(&elem_dim, &num_elem);
#if 0
  for(int i=0; i<num_vertex; i++)
  {
    double xyz[3];
    m3dc1_node_getcoord(&i, xyz);
    int ent_dim=0;
    int geom_class_dim;
    int geom_class_id;
    double normal[3];
    double curv;
    m3dc1_ent_getgeomclass ( &ent_dim, &i, &geom_class_dim, &geom_class_id);
    int bdy;
    m3dc1_node_isongeombdry(&i, &bdy);
    if(bdy) assert(geom_class_dim<2);
    if(geom_class_dim<2)
    {
      assert(bdy);
      m3dc1_node_getnormvec (&i, normal);
      m3dc1_node_getcurv (&i, &curv);
      cout<<" rank "<<PCU_Comm_Self()<<" node "<<i<<"("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<") geo class "<<geom_class_dim<<" "<<geom_class_id<<" normal curv "<<normal[0]<<" "<<normal[1]<<" "<<curv<<endl;
    }
    else
      cout<<" rank "<<PCU_Comm_Self()<<" node "<<i<<"("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<") geo class "<<geom_class_dim<<" "<<geom_class_id<<endl;
    // 2D mesh, z component =0
    if(num_plane==1) assert(AlmostEqualDoubles(xyz[2], 0, 1e-6, 1e-6));
  }
#endif

  int value_type[] = {scalar_type,scalar_type};

  int b_field=1, c_field=2, x_field=3;
#ifndef ASSEMBLEONLY
  int num_values=1;
#else
  int num_values=3;
#endif
  int num_dofs=6;
  if (num_plane>1) num_dofs=12;
  int num_dofs_node = num_values * num_dofs;
  m3dc1_field_create (&b_field, "b_field", &num_values, value_type, &num_dofs);
  m3dc1_field_create (&c_field, "c_field", &num_values, value_type, &num_dofs);
  m3dc1_field_create (&x_field, "x_field", &num_values, value_type, &num_dofs);
  m3dc1_field_printcompnorm(&b_field, "b_field init info");
  PetscMemoryGetCurrentUsage(&mem);
  PetscSynchronizedPrintf(MPI_COMM_WORLD, "process %d mem usage %f M \n ",PCU_Comm_Self(), mem/1e6);
  PetscSynchronizedFlush(MPI_COMM_WORLD, NULL);
  if(!PCU_Comm_Self()) cout<<" set b field ..."<<endl;
  // fill b field
  for(int inode=0; inode<num_vertex; inode++)
  {
    double xyz[3];
    m3dc1_node_getcoord(&inode, xyz);
    // 2D mesh, z component =0
    if(num_plane==1) assert(AlmostEqualDoubles(xyz[2], 0, 1e-6, 1e-6));
    vector<double> dofs(num_dofs_node*(1+scalar_type));
    for(int i=0; i<num_dofs_node*(1+scalar_type); i++)
      dofs.at(i)=xyz[i%3];
    m3dc1_ent_setdofdata(&vertex_dim, &inode, &b_field, &num_dofs_node, &dofs.at(0));
  }
  m3dc1_field_printcompnorm(&b_field, "b_field after set info");
  PetscMemoryGetCurrentUsage(&mem);
  PetscSynchronizedPrintf(MPI_COMM_WORLD, "process %d mem usage %f M \n ",PCU_Comm_Self(), mem/1e6);
  PetscSynchronizedFlush(MPI_COMM_WORLD, NULL);
  if(!PCU_Comm_Self()) cout<<" set matrix ..."<<endl; 
  t1 = MPI_Wtime();
  // fill matrix 
  // the matrix is diagnal dominant; thus should be positive definite
  int matrix_mult=1, matrix_solve=2;
  int matrix_mult_type = M3DC1_MULTIPLY;
  int matrix_solve_type = M3DC1_SOLVE;
#ifndef ASSEMBLEONLY
  m3dc1_matrix_create(&matrix_mult, &matrix_mult_type, value_type, &b_field);
#endif
  m3dc1_matrix_create(&matrix_solve, &matrix_solve_type, value_type, &b_field);
  double diag_value=2.0, off_diag=1.0;
  int node_elm = 3;
  if(num_plane>1) node_elm=6;
  int num_dofs_element = num_dofs_node*node_elm;
  vector<double> block(num_dofs_element*num_dofs_element*(1+scalar_type),0);
  for(int i=0; i<num_dofs_element; i++)
  {
    for(int j=0; j<num_dofs_element; j++)
    {
      double val= (i==j? diag_value: off_diag);
      if(!scalar_type) block.at(i*num_dofs_element+j)=val;
      else
      {
        block.at(2*i*num_dofs_element+2*j)=val;
        block.at(2*i*num_dofs_element+2*j+1)=off_diag;
      }
    }
  }
  //check geo
#ifndef ASSEMBLEONLY
  if(num_plane >1)
  {
    for( int ielm = 0; ielm < num_elem; ielm++)
    {
      int ent_dim=3, adj_dim=0, num_adj_ent, adj_ent_allocated_size=6;
      int nodes[6];
      m3dc1_ent_getadj (&ent_dim, &ielm, &adj_dim, nodes, &adj_ent_allocated_size, &num_adj_ent);
      assert(num_adj_ent==adj_ent_allocated_size);
      double normal1[3], normal2[3], curv1, curv2;
      for(int i=0; i<3; i++)
      {
        int is_bdy=0;
        m3dc1_node_isongeombdry(nodes+i, &is_bdy);
        if(is_bdy)
        {
          int geom_class_dim,geom_class_id;
          m3dc1_ent_getgeomclass (&adj_dim, nodes+i, &geom_class_dim, &geom_class_id);
          m3dc1_node_getnormvec(nodes+i, normal1);
          m3dc1_node_getcurv(nodes+i,&curv1);
          if(geom_class_dim==0)
          {
             cout<<" Node classified on geometric vertex with normal "<<normal1[0]<<" "<<normal1[1]<<" curv "<<curv1<<endl;
          }
          m3dc1_node_isongeombdry(nodes+i+3, &is_bdy);
          assert(is_bdy);
          m3dc1_node_getnormvec(nodes+i+3, normal2);
          m3dc1_node_getcurv(nodes+i+3,&curv2);
          for(int i=0; i<2; i++)
            assert(AlmostEqualDoubles(normal1[i], normal2[i], 1e-6, 1e-6));
          assert(AlmostEqualDoubles(curv1,curv2, 1e-6, 1e-6));
        }
      }
    }
  }
#endif
  for(int ielm = 0; ielm < num_elem; ielm++)
  {
    tf1 = MPI_Wtime();
    tf2 = MPI_Wtime();
    int nodes[256];
    int size_alloc=256;
    int num_nodes=-1;
    m3dc1_ent_getadj (&elem_dim, &ielm, &vertex_dim, nodes, &size_alloc, &num_nodes);
    for(int rowVar=0; rowVar< num_values; rowVar++)
    {
      for(int colVar=0; colVar< num_values; colVar++)
      {
         vector<double> block_tmp = block;
         if(rowVar!=colVar)
         {
           for(int i=0; i<block_tmp.size(); i++) block_tmp.at(i)*=0.5/num_values;
         }
         m3dc1_matrix_insertblock(&matrix_solve, &ielm, &rowVar, &colVar, &block_tmp[0]);
#ifndef ASSEMBLEONLY
         m3dc1_matrix_insertblock(&matrix_mult, &ielm, &rowVar, &colVar, &block_tmp[0]);
#endif
      }
    }
  }
  t2 = MPI_Wtime();
  PetscMemoryGetCurrentUsage(&mem);
  PetscSynchronizedPrintf(MPI_COMM_WORLD, "process %d mem usage %f M \n ",PCU_Comm_Self(), mem/1e6);
  PetscSynchronizedFlush(MPI_COMM_WORLD,NULL);
  if(!PCU_Comm_Self()) cout<<" assemble matrix ..."<<endl;
#ifndef ASSEMBLEONLY
  m3dc1_matrix_freeze(&matrix_mult);
#endif
  m3dc1_matrix_freeze(&matrix_solve); 
  //m3dc1_matrix_write(&matrix_mult, "matrixMult.m");
  //m3dc1_matrix_write(&matrix_solve, "matrixSolve.m");
  PetscMemoryGetCurrentUsage(&mem);
  PetscSynchronizedPrintf(MPI_COMM_WORLD, "process %d mem usage %f M \n ",PCU_Comm_Self(), mem/1e6);
  PetscSynchronizedFlush(MPI_COMM_WORLD, NULL);
  t3 = MPI_Wtime();
  // calculate c field
  //m3dc1_field_print(&b_field);
#ifndef ASSEMBLEONLY
  if(!PCU_Comm_Self()) cout<<" do matrix-vector multiply ..."<<endl;
  m3dc1_matrix_multiply(&matrix_mult, &b_field, &c_field); 
  //m3dc1_field_print(&c_field);
  t4 = MPI_Wtime();
  // let's test field operations here
  m3dc1_field_copy(&x_field, &c_field);
  double val[]={2.};
  int realtype=0;
  m3dc1_field_mult(&x_field, val, &realtype);
  m3dc1_field_add(&x_field, &c_field);
  for(int i=0; i<num_vertex; i++)
  {
    vector<double> dofs1(num_dofs_node*(1+scalar_type)), dofs2(num_dofs_node*(1+scalar_type));
    int num_dofs_t;
    m3dc1_ent_getdofdata(&vertex_dim, &i, &c_field, &num_dofs_t, &dofs1.at(0));
    assert(num_dofs_t==num_dofs_node);
    m3dc1_ent_getdofdata(&vertex_dim, &i, &x_field, &num_dofs_t, &dofs2.at(0));
    assert(num_dofs_t==num_dofs_node);
    for(int i=0; i<num_dofs_node*((1+scalar_type)); i++)
      assert(AlmostEqualDoubles(dofs2.at(i), (val[0]+1)*dofs1.at(i), 1e-6, 1e-6));
  }
  // copy c field to x field
  m3dc1_field_copy(&x_field, &c_field);
  //m3dc1_field_print(&x_field);
  if(!PCU_Comm_Self()) cout<<" solve ..."<<endl;
  // solve Ax=c
  m3dc1_matrix_solve(&matrix_solve, &x_field);
  //m3dc1_field_print(&x_field);
  t5 = MPI_Wtime();
  // verify x=b
  for(int inode=0; inode<num_vertex; inode++)
  {
    vector<double> dofs_x(num_dofs_node*(1+scalar_type)), dofs_b(num_dofs_node*(1+scalar_type));
    m3dc1_ent_getdofdata(&vertex_dim, &inode, &b_field, &num_dofs_node, &dofs_b.at(0));
    m3dc1_ent_getdofdata(&vertex_dim, &inode, &x_field, &num_dofs_node, &dofs_x.at(0));
    for(int idof=0; idof<num_dofs_node*(1+scalar_type); idof++)
      assert(AlmostEqualDoubles(dofs_b.at(idof),dofs_x.at(idof), 1e-3, 1e-3));
  }
#else
  t5=t4 = MPI_Wtime();
#endif
  if(!PCU_Comm_Self())
  {
    cout<<" time: fill matrix "<<t2-t1<<" assemble "<<t3-t2<<" mult "<<t4-t3<<" solve "<<t5-t4<<endl; 
    cout<<" mid flush "<<tf2-tf1<<endl;
  }
  m3dc1_matrix_delete(&matrix_mult);
  m3dc1_matrix_delete(&matrix_solve);
  m3dc1_field_delete(&x_field);
  m3dc1_field_delete(&b_field);
  m3dc1_field_delete(&c_field);

  PetscFinalize();
  m3dc1_scorec_finalize();
  MPI_Finalize();
  return M3DC1_SUCCESS;
}

bool AlmostEqualDoubles(double A, double B,
            double maxDiff, double maxRelDiff)
{
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

