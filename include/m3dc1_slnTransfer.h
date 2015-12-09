/****************************************************************************** 

  Copyright (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef M3DC1_SLNTRANSFER_H
#define M3DC1_SLNTRANSFER_H
#include "maSolutionTransfer.h"
#include <apf.h>
#include <apfVector.h>
#include <apfNumbering.h>
#include <iostream>
#include <vector>
#include <ReducedQuintic.h>
using namespace std;

class ReducedQuinticTransfer : public ma::SolutionTransfer
{
  public:
    ReducedQuinticTransfer(apf::Mesh* m, vector<apf::Field*>& fds, ReducedQuintic* shape): mesh(m), fields(fds), thecase(shape)
    {
      int maxComp=0;
      for(int i=0; i<fds.size(); i++)
        maxComp=max(maxComp,apf::countComponents(fds.at(i)));
      value.allocate(6*maxComp);
    }
    virtual void onVertex(
        apf:: MeshElement* parent,
        ma:: Vector const& xi,
        ma:: Entity* vert);
    virtual bool hasNodesOn(int dimension)
    {
      return dimension==0;
    }
  protected:
    vector<apf::Field*> fields;
    apf::Mesh* mesh;
    apf::NewArray<double> value;
    ReducedQuintic* thecase;
    static int dofNode;
};
#endif
