/****************************************************************************** 

  (c) 2005-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef M3DC1_GHOST_H
#define M3DC1_GHOST_H
#include "map"
#include <set>
#include "utility"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apf.h"
#include "m3dc1_scorec.h"
#include "m3dc1_field.h"

bool is_ent_original_ghost(apf::Mesh2* mesh, apf::MeshEntity* e);
int get_ent_ownpartid_ghost(apf::Mesh2* mesh, apf::MeshEntity* ent);
apf::MeshEntity* get_ent_owncopy_ghost(apf::Mesh2* mesh, apf::MeshEntity* ent);

class m3dc1_ghost
{
public:
  m3dc1_ghost();
  ~m3dc1_ghost();
  static m3dc1_ghost* instance();
  static void destroy();
  // functions
  void reset();
  void clean( std::set<int> & fields_keep);
  void initialize(); // to be called after creating ghost mesh.

  // data
  apf::Mesh2* mesh;
  bool is_ghosted;
  int nlayers;
  int ordering_opt;

  // local counter for fast info retrieval
  int num_local_ent[4];
  int num_global_ent[4];
  int num_own_ent[4];

  // field container 
  std::map<FieldID, m3dc1_field*>* field_container;

  // Populate fields by synchronizing with m3dc1_mesh's field
  // int field_synchronize();

  // tag for local entity id
  apf::MeshTag* local_entid_tag;

  // tag for owned partid attached to the part bdry entities
  apf::MeshTag* own_partid_tag; 

  // tags for second order adjanceny info
  apf::MeshTag* num_global_adj_node_tag;
  apf::MeshTag* num_own_adj_node_tag;
private:
  void set_node_adj_tag();
  static m3dc1_ghost* _instance;
};
#endif


