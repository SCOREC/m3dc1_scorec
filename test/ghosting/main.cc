/* Test access of fields on a ghosted mesh obtained via the workflow:
   APF -> omega_h -> ghost -> APF.
*/

#ifndef NOM3DC1
#define NOM3DC1
#endif

#include "m3dc1_ghost.h"
#include "m3dc1_mesh.h"
#include "m3dc1_model.h"
#include "apfOmega_h.h"
#include "m3dc1_scorec.h"
#include "m3dc1_mesh.h"
#include "gmi_mesh.h"
#include "gmi.h"
#include "gmi_null.h"
#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "PCU.h"
#include <cstdlib>
#include <iostream>
#include <assert.h>

using namespace std;

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  m3dc1_scorec_init();
  // PCU_Comm_Init();
  
  int op = 0, scalar_type = 0;
  int value_type[] = {scalar_type, scalar_type};
  int field_1 = 1, field_2 = 2, field_3 = 3;
  int num_values = 3;
  int num_dofs = 6;
  
  int num_dofs_node = num_values * num_dofs;
  
  int num_vertex, num_own_vertex, vertex_dim = 0;
  int num_elem, elem_dim = 2;

  int num_plane = 1;

  if (argc < 4 && !PCU_Comm_Self()) {
    printf("Usage: ./ghost_example model mesh #planes "
	   "real(0)/complex(1)\n");
    return M3DC1_FAILURE;
  }

  if (argc > 3) {
    num_plane = atoi(argv[3]);
    if (num_plane > 1 && PCU_Comm_Peers()%num_plane == 0)
      m3dc1_model_setnumplane (&num_plane);
  }
  
  if(argc > 4)
    scalar_type = atoi(argv[4]);

  if(num_plane>1)
    elem_dim=3; // use wedge in 3D

  m3dc1_model_load(argv[1]);
  // m3dc1_model_print();
  m3dc1_mesh_load(argv[2]);
  // printStats(m3dc1_mesh::instance()->mesh);

  m3dc1_mesh_getnument (&vertex_dim, &num_vertex);
  m3dc1_mesh_getnumownent (&vertex_dim, &num_own_vertex);
  m3dc1_mesh_getnument(&elem_dim, &num_elem);

  // Create and fill fields on APF mesh
  m3dc1_field_create(&field_1,
		     "field_1",
		     &num_values,
		     value_type,
		     &num_dofs);
  m3dc1_field_create(&field_2,
		     "field_2",
		     &num_values,
		     value_type,
		     &num_dofs);
  m3dc1_field_create(&field_3,
		     "field_3",
		     &num_values,
		     value_type,
		     &num_dofs);


  // Load mesh using null model
  gmi_register_null();
  apf::Mesh2* mesh = apf::loadMdsMesh(".null", argv[2]);
  mesh->verify();

  // APF -> omega_h -> APF
  osh_t om = osh::fromAPF(mesh);
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  mesh = apf::makeEmptyMdsMesh(gmi_load(".null"), osh_dim(om), false);
  /*  mesh = apf::makeEmptyMdsMesh(m3dc1_model->instance()->model,
      osh_dim(om), false);*/
  osh_ghost(&om, 1);
  osh::toAPF(om, mesh);
  // mesh->verify();

  // Clean up and finalize
  osh_free(om);
  mesh->destroyNative();
  apf::destroyMesh(mesh);

  m3dc1_field_delete(&field_1);
  m3dc1_field_delete(&field_2);
  m3dc1_field_delete(&field_3);
  
  m3dc1_scorec_finalize();
  // PCU_Comm_Free();
  MPI_Finalize();
  return M3DC1_SUCCESS;}





