/****************************************************************************** 

  (c) 2005-2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef M3DC1_MESHGEN
#ifdef M3DC1_TRILINOS
#ifndef M3DC1_LS_H
#define M3DC1_LS_H
#include <map>
#include "apf.h"
#include "m3dc1_scorec.h"

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
typedef int global_ordinal_type;
#else
//typedef long long global_ordinal_type; -- 64bit should throw compilation error
#endif

// NOTE: all field realted interaction is done through m3dc1 api rather than apf
class m3dc1_epetra
{
public:
  m3dc1_epetra(int i, int t, int s, FieldID field);
  ~m3dc1_epetra();
  int initialize(); // create a matrix and solver object
  int destroy(); // delete a matrix and solver object
  int set_value(int row, int col, int operation, double real_val, double imag_val); //insertion/addition with global numbering
  // values use row-wise, rsize * csize block
  int add_values(int rsize, int * rows, int csize, int * columns, double* values);
  int get_scalar_type() { return scalar_type; }
  int get_field_id() { return field_id;}
  int assemble();
  Epetra_CrsMatrix* epetra_mat;
  Epetra_Map* _overlap_map;
  Epetra_Map* _owned_map;
  global_ordinal_type nge;
  int num_solver_iter;
  int matrix_type; /* 0 for multiplication, 1 for solver*/
protected:
  int id;
  int scalar_type;
  apf::Field* _field;
  int field_id; // the field that provides numbering
};

class m3dc1_ls
{
public:
// functions
  m3dc1_ls(){matrix_container = new std::map<int, m3dc1_epetra*>;}
  ~m3dc1_ls();
  static m3dc1_ls* instance();
  m3dc1_epetra* get_matrix(int matrix_id);
  void add_matrix(int matrix_id, m3dc1_epetra*);
// data
  std::map<int, m3dc1_epetra*>* matrix_container;
private:
  static m3dc1_ls* _instance;
};
#endif
#endif
#endif //#ifndef M3DC1_MESHGEN
