#include <cassert>

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <sstream>

#include <YapInterface.h>

#include "ParfactorList.h"
#include "FactorGraph.h"
#include "LiftedOperations.h"
#include "LiftedVe.h"
#include "VarElim.h"
#include "LiftedBp.h"
#include "CountingBp.h"
#include "BeliefProp.h"
#include "LiftedKc.h"
#include "ElimGraph.h"
#include "BayesBall.h"


namespace Horus {

namespace {

Parfactor* readParfactor (YAP_Term);

ObservedFormulas* readLiftedEvidence (YAP_Term);

std::vector<unsigned> readUnsignedList (YAP_Term);

Params readParameters (YAP_Term);

YAP_Term fillSolutionList (const std::vector<Params>&);

}

typedef std::pair<ParfactorList*, ObservedFormulas*> LiftedNetwork;



int
createLiftedNetwork()
{
  Parfactors parfactors;
  YAP_Term parfactorList = YAP_ARG1;
  while (parfactorList != YAP_TermNil()) {
    YAP_Term pfTerm = YAP_HeadOfTerm (parfactorList);
    //Parfactor * pf=
    //readParfactor (pfTerm);
     //std::cout<< " vvv\n";std::cout.flush();

     //parfactors.push_back (pf);
    parfactors.push_back (readParfactor (pfTerm));
     //std::cout<< " wwww\n";std::cout.flush();

    parfactorList = YAP_TailOfTerm (parfactorList);
  }

  // LiftedUtils::printSymbolDictionary();
  if (Globals::verbosity > 2) {
    Util::printHeader ("INITIAL PARFACTORS");
    for (size_t i = 0; i < parfactors.size(); i++) {
      parfactors[i]->print();
      std::cout << std::endl;
    }
  }

  ParfactorList* pfList = new ParfactorList (parfactors);

  if (Globals::verbosity > 2) {
    Util::printHeader ("SHATTERED PARFACTORS");
    pfList->print();
  }
  
  // read evidence
  ObservedFormulas* obsFormulas = readLiftedEvidence (YAP_ARG2);

  LiftedNetwork* network = new LiftedNetwork (pfList, obsFormulas);

  YAP_Int p = (YAP_Int) (network);
  return YAP_Unify (YAP_MkIntTerm (p), YAP_ARG3);
}



int
createGroundNetwork()
{
  std::string factorsType ((char*) YAP_AtomName (YAP_AtomOfTerm (YAP_ARG1)));
  FactorGraph* fg = new FactorGraph();
  if (factorsType == "bayes") {
    fg->setFactorsAsBayesian();
  }
  YAP_Term factorList = YAP_ARG2;
  while (factorList != YAP_TermNil()) {
    YAP_Term factor = YAP_HeadOfTerm (factorList);
    // read the var ids
    VarIds varIds = readUnsignedList (YAP_ArgOfTerm (1, factor));
    // read the ranges
    Ranges ranges = readUnsignedList (YAP_ArgOfTerm (2, factor));
    // read the parameters
    Params params = readParameters (YAP_ArgOfTerm (3, factor));
    // read dist id
    unsigned distId = (unsigned) YAP_IntOfTerm (YAP_ArgOfTerm (4, factor));
    fg->addFactor (Factor (varIds, ranges, params, distId));
    factorList = YAP_TailOfTerm (factorList);
  }
  unsigned nrObservedVars = 0;
  YAP_Term evidenceList = YAP_ARG3;
  while (evidenceList != YAP_TermNil()) {
    YAP_Term evTerm = YAP_HeadOfTerm (evidenceList);
    unsigned vid = (unsigned) YAP_IntOfTerm ((YAP_ArgOfTerm (1, evTerm)));
    unsigned ev  = (unsigned) YAP_IntOfTerm ((YAP_ArgOfTerm (2, evTerm)));
    assert (fg->getVarNode (vid));
    fg->getVarNode (vid)->setEvidence (ev);
    evidenceList = YAP_TailOfTerm (evidenceList);
    nrObservedVars ++;
  }
  if (FactorGraph::exportToLibDai()) {
    fg->exportToLibDai ("model.fg");
  }
  if (FactorGraph::exportToUai()) {
    fg->exportToUai ("model.uai");
  }
  if (FactorGraph::exportGraphViz()) {
    fg->exportToGraphViz ("model.dot");
  }
  if (FactorGraph::printFactorGraph()) {
    fg->print();
  }
  if (Globals::verbosity > 0) {
    std::cout << "factor graph contains " ;
    std::cout << fg->nrVarNodes() << " variables and " ;
    std::cout << fg->nrFacNodes() << " factors " << std::endl;
  }
  YAP_Int p = (YAP_Int) (fg);
  return YAP_Unify (YAP_MkIntTerm (p), YAP_ARG4);
}



int
runLiftedSolver()
{
  LiftedNetwork* network = (LiftedNetwork*) YAP_IntOfTerm (YAP_ARG1);
  ParfactorList copy (*network->first);
  LiftedOperations::absorveEvidence (copy, *network->second);

  LiftedSolver* solver = 0;
  switch (Globals::liftedSolver) {
    case LiftedSolverType::lveSolver: solver = new LiftedVe (copy); break;
    case LiftedSolverType::lbpSolver: solver = new LiftedBp (copy); break;
    case LiftedSolverType::lkcSolver: solver = new LiftedKc (copy); break;
  }

  if (Globals::verbosity > 0) {
    solver->printSolverFlags();
    std::cout << std::endl;
  }

  YAP_Term taskList = YAP_ARG2;
  std::vector<Params> results;
  while (taskList != YAP_TermNil()) {
    Grounds queryVars;
    YAP_Term jointList = YAP_HeadOfTerm (taskList);
    while (jointList != YAP_TermNil()) {
      YAP_Term ground = YAP_HeadOfTerm (jointList);
      if (YAP_IsAtomTerm (ground)) {
        std::string name ((char*) YAP_AtomName (YAP_AtomOfTerm (ground)));
        queryVars.push_back (Ground (LiftedUtils::getSymbol (name)));
      } else {
        assert (YAP_IsApplTerm (ground));
        YAP_Functor yapFunctor = YAP_FunctorOfTerm (ground);
        std::string name ((char*) (YAP_AtomName (
            YAP_NameOfFunctor (yapFunctor))));
        unsigned arity = (unsigned) YAP_ArityOfFunctor (yapFunctor);
        Symbol functor = LiftedUtils::getSymbol (name);
        Symbols args;
        for (unsigned i = 1; i <= arity; i++) {
          YAP_Term ti = YAP_ArgOfTerm (i, ground);
          assert (YAP_IsAtomTerm (ti));
          std::string arg ((char *) YAP_AtomName (YAP_AtomOfTerm (ti)));
          args.push_back (LiftedUtils::getSymbol (arg));
        }
        queryVars.push_back (Ground (functor, args));
      }
      jointList = YAP_TailOfTerm (jointList);
    }
    results.push_back (solver->solveQuery (queryVars));
    taskList = YAP_TailOfTerm (taskList);
  }

  delete solver;

  return YAP_Unify (fillSolutionList (results), YAP_ARG3);
}



int
runGroundSolver()
{
  FactorGraph* fg = (FactorGraph*) YAP_IntOfTerm (YAP_ARG1);

  std::vector<VarIds> tasks;
  YAP_Term taskList = YAP_ARG2;
  while (taskList != YAP_TermNil()) {
    tasks.push_back (readUnsignedList (YAP_HeadOfTerm (taskList)));
    taskList = YAP_TailOfTerm (taskList);
  }

  FactorGraph* mfg = fg;
  if (fg->bayesianFactors()) {
    std::set<VarId> vids;
    for (size_t i = 0; i < tasks.size(); i++) {
      Util::addToSet (vids, tasks[i]);
    }
    mfg = BayesBall::getMinimalFactorGraph (
        *fg, VarIds (vids.begin(), vids.end()));
  }

  GroundSolver* solver = 0;
  CountingBp::setFindIdenticalFactorsFlag (false);
  switch (Globals::groundSolver) {
    case GroundSolverType::veSolver:  solver = new VarElim    (*mfg); break;
    case GroundSolverType::bpSolver:  solver = new BeliefProp (*mfg); break;
    case GroundSolverType::CbpSolver: solver = new CountingBp (*mfg); break;
  }

  if (Globals::verbosity > 0) {
    solver->printSolverFlags();
    std::cout << std::endl;
  }

  std::vector<Params> results;
  results.reserve (tasks.size());
  for (size_t i = 0; i < tasks.size(); i++) {
    results.push_back (solver->solveQuery (tasks[i]));
  }

  delete solver;
  if (fg->bayesianFactors()) {
    delete mfg;
  }

  return YAP_Unify (fillSolutionList (results), YAP_ARG3);
}



int
setParfactorsParams()
{
  LiftedNetwork* network = (LiftedNetwork*) YAP_IntOfTerm (YAP_ARG1);
  ParfactorList* pfList = network->first;
  YAP_Term distIdsList = YAP_ARG2;
  YAP_Term paramsList  = YAP_ARG3;
  std::unordered_map<unsigned, Params> paramsMap;
  while (distIdsList != YAP_TermNil()) {
    unsigned distId = (unsigned) YAP_IntOfTerm (
        YAP_HeadOfTerm (distIdsList));
    assert (Util::contains (paramsMap, distId) == false);
    paramsMap[distId] = readParameters (YAP_HeadOfTerm (paramsList));
    distIdsList = YAP_TailOfTerm (distIdsList);
    paramsList  = YAP_TailOfTerm (paramsList);
  }
  ParfactorList::iterator it = pfList->begin();
  while (it != pfList->end()) {
    assert (Util::contains (paramsMap, (*it)->distId()));
    (*it)->setParams (paramsMap[(*it)->distId()]);
    ++ it;
  }
  return TRUE;
}



int
setFactorsParams()
{
  FactorGraph* fg = (FactorGraph*) YAP_IntOfTerm (YAP_ARG1);
  YAP_Term distIdsList = YAP_ARG2;
  YAP_Term paramsList  = YAP_ARG3;
  std::unordered_map<unsigned, Params> paramsMap;
  while (distIdsList != YAP_TermNil()) {
    unsigned distId = (unsigned) YAP_IntOfTerm (
        YAP_HeadOfTerm (distIdsList));
    assert (Util::contains (paramsMap, distId) == false);
    paramsMap[distId] = readParameters (YAP_HeadOfTerm (paramsList));
    distIdsList = YAP_TailOfTerm (distIdsList);
    paramsList  = YAP_TailOfTerm (paramsList);
  }
  const FacNodes& facNodes = fg->facNodes();
  for (size_t i = 0; i < facNodes.size(); i++) {
    unsigned distId = facNodes[i]->factor().distId();
    assert (Util::contains (paramsMap, distId));
    facNodes[i]->factor().setParams (paramsMap[distId]);
  }
  return TRUE;
}



int
setVarsInformation()
{
  Var::clearVarsInfo();
  std::vector<std::string> labels;
  YAP_Term labelsL = YAP_ARG1;
  while (labelsL != YAP_TermNil()) {
    YAP_Atom atom = YAP_AtomOfTerm (YAP_HeadOfTerm (labelsL));
    labels.push_back ((char*) YAP_AtomName (atom));
    labelsL = YAP_TailOfTerm (labelsL);
  }
  unsigned count = 0;
  YAP_Term stateNamesL = YAP_ARG2;
  while (stateNamesL != YAP_TermNil()) {
    States states;
    YAP_Term namesL = YAP_HeadOfTerm (stateNamesL);
    while (namesL != YAP_TermNil()) {
      YAP_Atom atom = YAP_AtomOfTerm (YAP_HeadOfTerm (namesL));
      states.push_back ((char*) YAP_AtomName (atom));
      namesL = YAP_TailOfTerm (namesL);
    }
    Var::addVarInfo (count, labels[count], states);
    count ++;
    stateNamesL = YAP_TailOfTerm (stateNamesL);
  }
  return TRUE;
}



int
setHorusFlag()
{
  std::string option ((char*) YAP_AtomName (YAP_AtomOfTerm (YAP_ARG1)));
  std::string value;
  if (option == "verbosity") {
    std::stringstream ss;
    ss << (int) YAP_IntOfTerm (YAP_ARG2);
    ss >> value;
  } else if (option == "bp_accuracy") {
    std::stringstream ss;
    ss << (float) YAP_FloatOfTerm (YAP_ARG2);
    ss >> value;
  } else if (option == "bp_max_iter") {
    std::stringstream ss;
    ss << (int) YAP_IntOfTerm (YAP_ARG2);
    ss >> value;
  } else {
    value = ((char*) YAP_AtomName (YAP_AtomOfTerm (YAP_ARG2)));
  }
  return Util::setHorusFlag (option, value);
}



int
freeGroundNetwork()
{
  delete (FactorGraph*) YAP_IntOfTerm (YAP_ARG1);
  return TRUE;
}



int
freeLiftedNetwork()
{
  LiftedNetwork* network = (LiftedNetwork*) YAP_IntOfTerm (YAP_ARG1);
  delete network->first;
  delete network->second;
  delete network;
  return TRUE;
}



namespace {

Parfactor*
readParfactor (YAP_Term pfTerm)
{
  // read dist id
  unsigned distId = YAP_IntOfTerm (YAP_ArgOfTerm (2, pfTerm));
  size_t indexConv=-1;
  std::vector<bool> conv;
  // read the ranges
  Ranges ranges;
  YAP_Term rangeList = YAP_ArgOfTerm (4, pfTerm);
  while (rangeList != YAP_TermNil()) {
    unsigned range = (unsigned) YAP_IntOfTerm (YAP_HeadOfTerm (rangeList));
    ranges.push_back (range);
    rangeList = YAP_TailOfTerm (rangeList);
  }
  
  YAP_Term typeTerm = YAP_ArgOfTerm (1, pfTerm);
  std::string nameType;
  
  if (YAP_IsAtomTerm(typeTerm)) {
      YAP_Atom yapF = YAP_AtomOfTerm (typeTerm);
      nameType = ((char*) YAP_AtomName (yapF));
      
  } else {
        YAP_Functor yapF = YAP_FunctorOfTerm (typeTerm);
        nameType = ((char*) YAP_AtomName (YAP_NameOfFunctor(yapF)));
  }
  
  //std::cout<<"name: "<<nameType<<"\n";
  bool dep =  !nameType.compare("deputy");
  //std::cout<<"deputy? " << dep <<"\n";
  if (dep)
  {
    indexConv= (size_t) 1;
    //index=(size_t) YAP_IntOfTerm (YAP_ArgOfTerm(2,typeTerm));
    //std::cout<<" -> indexConv: " <<indexConv <<"\n";
  }
  
  bool nor= !YAP_IsAtomTerm (typeTerm);
  //std::cout<<"noisy_or? " << nor <<"\n";
  
  if (nor)
  {
    indexConv= (size_t) YAP_IntOfTerm (YAP_ArgOfTerm(1,typeTerm));
    //index=(size_t) YAP_IntOfTerm (YAP_ArgOfTerm(2,typeTerm));
    //std::cout<<" -> indexConv: " <<indexConv <<"\n";
  }
  //if(dep)  indexConv=1;
  // read parametric random vars
  ProbFormulas formulas;
  unsigned count = 0;
  std::unordered_map<YAP_Term, LogVar> lvMap;
  YAP_Term pvList = YAP_ArgOfTerm (3, pfTerm);
  while (pvList != YAP_TermNil()) {
    //std::cout << " entering cycle\n"; std::cout.flush();

    YAP_Term formulaTerm = YAP_HeadOfTerm (pvList);
    if (YAP_IsAtomTerm (formulaTerm)) {
      std::string name ((char*) YAP_AtomName (YAP_AtomOfTerm (formulaTerm)));
      //std::cout << "  first if (atomic): " << name; std::cout<< "\n";
      Symbol functor = LiftedUtils::getSymbol (name);
      //std::cout<< "   nor, dep, count, indexConv " << nor << " " << dep << " " << count << " " << (int) indexConv << "\n";
      if ((nor || dep) && (count==indexConv))
         conv.push_back(true);
      //formulas.push_back (ProbFormula (functor, ranges[count]));
      
      else
        conv.push_back(false);
      
      //std::cout << "   conv : " << conv << "\n";
      formulas.push_back(ProbFormula (functor, ranges[count]));
      //std::cout<< "   ac\n";

    } else {
      LogVars logVars,orLogVars;
      YAP_Functor yapFunctor = YAP_FunctorOfTerm (formulaTerm);
      std::string name ((char*) YAP_AtomName (
      YAP_NameOfFunctor (yapFunctor)));
      //std::cout << "  second if (not_atomic): " << name;  std::cout<< "\n"; std::cout.flush();
/*
    if (nor && id==index)
    {
	std::cout<< "or"; std::cout.flush();
    YAP_Term varsList = YAP_ArgOfTerm (2, formulaTerm);
    
    j = YAP_IntOfTerm(YAP_ArgOfTerm (1, formulaTerm));
    std::cout<<" "<<j; std::cout.flush();

    YAP_Term prv = YAP_ArgOfTerm (3, formulaTerm);
    YAP_Functor yapFunctorprv = YAP_FunctorOfTerm (prv);
	std::string nameprv ((char*) YAP_AtomName (
          YAP_NameOfFunctor (yapFunctorprv)));
	Symbol functor = LiftedUtils::getSymbol (nameprv);
      unsigned arity = (unsigned) YAP_ArityOfFunctor (yapFunctorprv);
      std::cout<< "arity"<< arity<< "\n";
      for (unsigned i = 1; i <= arity; i++) {
        YAP_Term ti = YAP_ArgOfTerm (i, prv);
        std::unordered_map<YAP_Term, LogVar>::iterator it = lvMap.find (ti);
        if (it != lvMap.end()) {
          logVars.push_back (it->second);
        } else {
         unsigned newLv = lvMap.size();
         lvMap[ti] = newLv;
         logVars.push_back (newLv);
        }}
	std::cout<<"log v "<<logVars<<"\n";
     while (varsList !=YAP_TermNil())
    {
      YAP_Term ti=YAP_HeadOfTerm(varsList);
        std::unordered_map<YAP_Term, LogVar>::iterator it = lvMap.find (ti);
        if (it != lvMap.end()) {
          orLogVars.push_back (it->second);
        } else {
         unsigned newLv = lvMap.size();
         lvMap[ti] = newLv;
         orLogVars.push_back (newLv);
        }
    varsList=YAP_TailOfTerm(varsList);
    }
      std::cout<< " "<< j<<" j\n";std::cout.flush();
      //ProbFormula pf=ProbFormula (functor, logVars,orLogVars, ranges[count],j);
      formulas.push_back (ProbFormula (functor, logVars,orLogVars, ranges[count],j));
//        formulas.push_back (ProbFormula (functor, logVars, ranges[count]));
	std::cout<< " pf\n";std::cout.flush();

//        formulas.push_back(pf);
//	std::cout<< " pf1\n";std::cout.flush();
	
    }

    
    else 
    
    {*/
      Symbol functor = LiftedUtils::getSymbol (name);
      unsigned arity = (unsigned) YAP_ArityOfFunctor (yapFunctor);
      //std::cout <<"    -> arity: " << arity << "\n";
      for (unsigned i = 1; i <= arity; i++) {
        YAP_Term ti = YAP_ArgOfTerm (i, formulaTerm);
        std::unordered_map<YAP_Term, LogVar>::iterator it = lvMap.find (ti);
        if (it != lvMap.end()) {
          logVars.push_back (it->second);
        } else {
         unsigned newLv = lvMap.size();
         lvMap[ti] = newLv;
         logVars.push_back (newLv);
        }
      }
            
      //std::cout<< "   nor, dep, count, indexConv " << nor << " " << dep << " " << count << " " << indexConv << "\n";
      if ((nor || dep) && (count==indexConv))
        conv.push_back(true);
      //formulas.push_back (ProbFormula (functor, ranges[count]));
      
      else
        conv.push_back(false);
 
      //std::cout << "   conv: " << conv << "\n";
/*       if (nor && (count==indexConv))
      { //std::cout<< "indexconv "<<indexConv<<"\n";
      formulas.push_back (ProbFormula (functor, logVars, ranges[count],indexConv));
      }
      else*/
      formulas.push_back (ProbFormula (functor, logVars, ranges[count]));
    //}
    }
    count ++;
    pvList = YAP_TailOfTerm (pvList);
  }
 //std::cout<< " pppp\n";std::cout.flush();
 
  // read the parameters
  Params params = readParameters (YAP_ArgOfTerm (5, pfTerm));

  // read the constraint
  Tuples tuples;
  if (lvMap.size() >= 1) {
    YAP_Term tupleList = YAP_ArgOfTerm (6, pfTerm);
    while (tupleList != YAP_TermNil()) {
      YAP_Term term = YAP_HeadOfTerm (tupleList);
      assert (YAP_IsApplTerm (term));
      YAP_Functor yapFunctor = YAP_FunctorOfTerm (term);
      unsigned arity = (unsigned) YAP_ArityOfFunctor (yapFunctor);
      assert (lvMap.size() == arity);
      Tuple tuple (arity);
      for (unsigned i = 1; i <= arity; i++) {
        YAP_Term ti = YAP_ArgOfTerm (i, term);
        if (YAP_IsAtomTerm (ti) == false) {
          std::cerr << "Error: the constraint contains free variables." ;
          std::cerr << std::endl;
          exit (EXIT_FAILURE);
        }
        std::string name ((char*) YAP_AtomName (YAP_AtomOfTerm (ti)));
	//std::cout<< name;std::cout<<" tuple";
        tuple[i - 1] = LiftedUtils::getSymbol (name);
      }
      tuples.push_back (tuple);
      tupleList = YAP_TailOfTerm (tupleList);
    }
  }
  //std::cout<< "returning new Parfactor\n";std::cout.flush();

  return new Parfactor (formulas, params, tuples, distId, conv , nor);
  //std::cout<<"pf lvs "<<pf->logVarSet()<<"\n";
}



ObservedFormulas*
readLiftedEvidence (YAP_Term observedList)
{
  ObservedFormulas* obsFormulas = new ObservedFormulas();
  while (observedList != YAP_TermNil()) {
    YAP_Term pair = YAP_HeadOfTerm (observedList);
    YAP_Term ground = YAP_ArgOfTerm (1, pair);
    Symbol functor;
    Symbols args;
    if (YAP_IsAtomTerm (ground)) {
      std::string name ((char*) YAP_AtomName (YAP_AtomOfTerm (ground)));
      functor = LiftedUtils::getSymbol (name);
    } else {
      assert (YAP_IsApplTerm (ground));
      YAP_Functor yapFunctor = YAP_FunctorOfTerm (ground);
      std::string name ((char*) (YAP_AtomName (
          YAP_NameOfFunctor (yapFunctor))));
      functor = LiftedUtils::getSymbol (name);
      unsigned arity = (unsigned) YAP_ArityOfFunctor (yapFunctor);
      for (unsigned i = 1; i <= arity; i++) {
        YAP_Term ti = YAP_ArgOfTerm (i, ground);
        assert (YAP_IsAtomTerm (ti));
        std::string arg ((char *) YAP_AtomName (YAP_AtomOfTerm (ti)));
        args.push_back (LiftedUtils::getSymbol (arg));
      }
    }
    unsigned evidence = (unsigned) YAP_IntOfTerm (YAP_ArgOfTerm (2, pair));
    bool found = false;
    for (size_t i = 0; i < obsFormulas->size(); i++) {
      if ((*obsFormulas)[i].functor()  == functor &&
          (*obsFormulas)[i].arity()    == args.size() &&
          (*obsFormulas)[i].evidence() == evidence) {
        (*obsFormulas)[i].addTuple (args);
        found = true;
      }
    }
    if (found == false) {
      obsFormulas->push_back (ObservedFormula (functor, evidence, args));
    }
    observedList = YAP_TailOfTerm (observedList);
  }
  return obsFormulas;
}



std::vector<unsigned>
readUnsignedList (YAP_Term list)
{
  std::vector<unsigned> vec;
  while (list != YAP_TermNil()) {
    vec.push_back ((unsigned) YAP_IntOfTerm (YAP_HeadOfTerm (list)));
    list = YAP_TailOfTerm (list);
  }
  return vec;
}



Params
readParameters (YAP_Term paramL)
{
  Params params;
  assert (YAP_IsPairTerm (paramL));
  while (paramL != YAP_TermNil()) {
    YAP_Term hd = YAP_HeadOfTerm (paramL);
    if (YAP_IsFloatTerm (hd)) {
      params.push_back ((double) YAP_FloatOfTerm (hd));
    } else {
      params.push_back ((double) YAP_IntOfTerm (hd));
    }
    paramL = YAP_TailOfTerm (paramL);
  }
  if (Globals::logDomain) {
    Util::log (params);
  }
  return params;
}



YAP_Term
fillSolutionList (const std::vector<Params>& results)
{
  YAP_Term list = YAP_TermNil();
  for (size_t i = results.size(); i-- > 0; ) {
    const Params& beliefs = results[i];
    YAP_Term queryBeliefsL = YAP_TermNil();
    for (size_t j = beliefs.size(); j-- > 0; ) {
      YAP_Int sl       = YAP_InitSlot (list);
      YAP_Term belief  = YAP_MkFloatTerm (beliefs[j]);
      queryBeliefsL    = YAP_MkPairTerm (belief, queryBeliefsL);
      list             = YAP_GetFromSlot (sl);
      YAP_RecoverSlots (1);
    }
    list = YAP_MkPairTerm (queryBeliefsL, list);
  }
  return list;
}

}



extern "C" void
init_predicates()
{
  YAP_UserCPredicate ("cpp_create_lifted_network",
      createLiftedNetwork, 3);

  YAP_UserCPredicate ("cpp_create_ground_network",
      createGroundNetwork, 4);

  YAP_UserCPredicate ("cpp_run_lifted_solver",
      runLiftedSolver, 3);

  YAP_UserCPredicate ("cpp_run_ground_solver",
      runGroundSolver, 3);

  YAP_UserCPredicate ("cpp_set_parfactors_params",
      setParfactorsParams, 3);

  YAP_UserCPredicate ("cpp_set_factors_params",
      setFactorsParams, 3);

  YAP_UserCPredicate ("cpp_set_vars_information",
      setVarsInformation, 2);

  YAP_UserCPredicate ("cpp_set_horus_flag",
      setHorusFlag, 2);

  YAP_UserCPredicate ("cpp_free_lifted_network",
      freeLiftedNetwork, 1);

  YAP_UserCPredicate ("cpp_free_ground_network",
     freeGroundNetwork, 1);
}

}  // namespace Horus

