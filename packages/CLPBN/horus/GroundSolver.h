#ifndef HORUS_GROUNDSOLVER_H
#define HORUS_GROUNDSOLVER_H

#include <iomanip>

#include "FactorGraph.h"
#include "Horus.h"


using namespace std;

class GroundSolver
{
  public:
    GroundSolver (const FactorGraph& factorGraph) : fg(factorGraph) { }

    virtual ~GroundSolver() { } // ensure that subclass destructor is called

    virtual Params solveQuery (VarIds queryVids) = 0;

    virtual void printSolverFlags (void) const = 0;

    void printAnswer (const VarIds& vids);

    void printAllPosterioris (void);

    static Params getJointByConditioning (GroundSolverType,
        FactorGraph, const VarIds& jointVarIds);

  protected:
    const FactorGraph& fg;

    DISALLOW_COPY_AND_ASSIGN (GroundSolver);
};

#endif // HORUS_GROUNDSOLVER_H

