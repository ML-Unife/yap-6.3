/*************************************************************************
*									 *
*	 YAP Prolog 							 *
*									 *
*	Yap Prolog was developed at NCCUP - Universidade do Porto	 *
*									 *
* Copyright L.Damas, V.S.Costa and Universidade do Porto 1985-1997	 *
*									 *
**************************************************************************
*									 *
* File:		eval.c							 *
* Last rev:								 *
* mods:									 *
* comments:	arithmetical expression evaluation			 *
*									 *
*************************************************************************/
#ifdef SCCS
static char     SccsId[] = "%W% %G%";
#endif

/*
 * This file implements arithmetic operations 
 *
 */
#include "Yap.h"
#include "Yatom.h"
#include "Heap.h"
#include "eval.h"

yap_error_number _YAP_matherror = YAP_NO_ERROR;

#define E_FUNC   blob_type
#define E_ARGS   arith_retptr o
#define USE_E_ARGS   o
#define RBIG(v)     (o)->big = v; return(big_int_e)

#define RINT(v)     (o)->Int = v; return(long_int_e)
#define RFLOAT(v)   (o)->dbl = v; return(double_e)
#define RERROR()    return(db_ref_e)

static Term
EvalToTerm(blob_type bt, union arith_ret *res)
{
  switch (bt) {
  case long_int_e:
    return(MkIntegerTerm(res->Int));
  case double_e:
    return(MkFloatTerm(res->dbl));
#ifdef USE_GMP
  case big_int_e:
    return(_YAP_MkBigIntTerm(res->big));
#endif
  default:
    return(TermNil);
  }
}

static E_FUNC
Eval(Term t, E_ARGS)
{
  if (IsVarTerm(t)) {
    _YAP_Error(INSTANTIATION_ERROR,TermNil,"in arithmetic");
    P = (yamop *)FAILCODE;
    RERROR();
  }
  if (IsApplTerm(t)) {
    Functor fun = FunctorOfTerm(t);
    switch ((CELL)fun) {
    case (CELL)FunctorLongInt:
      RINT(LongIntOfTerm(t));
    case (CELL)FunctorDouble:
      RFLOAT(FloatOfTerm(t));
#ifdef USE_GMP
    case (CELL)FunctorBigInt:
      RBIG(_YAP_BigIntOfTerm(t));
#endif
    default:
      {
	Int n = ArityOfFunctor(fun);
	Atom name  = NameOfFunctor(fun);
	ExpEntry *p;

	if (EndOfPAEntr(p = RepExpProp(_YAP_GetExpProp(name, n)))) {
	  Term ti[2];

	  /* error */
	  ti[0] = t;
	  ti[1] = MkIntegerTerm(n);
	  t = _YAP_MkApplTerm(_YAP_MkFunctor(_YAP_LookupAtom("/"),2), 2, ti);
	  _YAP_Error(TYPE_ERROR_EVALUABLE, t,
		"functor %s/%d for arithmetic expression",
		RepAtom(name)->StrOfAE,n);
	  P = (yamop *)FAILCODE;
	  RERROR();
	}
	if (n == 1)
	  return(p->FOfEE.unary(ArgOfTerm(1,t), USE_E_ARGS));
	return(p->FOfEE.binary(ArgOfTerm(1,t),ArgOfTerm(2,t), USE_E_ARGS));
      }
    }
  } else if (IsPairTerm(t)) {
    return(Eval(HeadOfTerm(t), USE_E_ARGS));
  } else if (IsIntTerm(t)) {
    RINT(IntOfTerm(t));
  } else {
    Atom name = AtomOfTerm(t);
    ExpEntry *p;

    if (EndOfPAEntr(p = RepExpProp(_YAP_GetExpProp(name, 0)))) {
      /* error */
      _YAP_Error(TYPE_ERROR_EVALUABLE, t,
	    "atom %s for arithmetic expression",
	    RepAtom(name)->StrOfAE);
      P = (yamop *)FAILCODE;
      RERROR();
    }
    return(p->FOfEE.constant(USE_E_ARGS));
  }
}

E_FUNC
_YAP_Eval(Term t, E_ARGS)
{
  if (IsVarTerm(t)) {
    _YAP_Error(INSTANTIATION_ERROR,TermNil,"in arithmetic");
    P = (yamop *)FAILCODE;
    RERROR();
  }
  if (IsApplTerm(t)) {
    Functor fun = FunctorOfTerm(t);
    switch ((CELL)fun) {
    case (CELL)FunctorLongInt:
      RINT(LongIntOfTerm(t));
    case (CELL)FunctorDouble:
      RFLOAT(FloatOfTerm(t));
#ifdef USE_GMP
    case (CELL)FunctorBigInt:
      RBIG(_YAP_BigIntOfTerm(t));
#endif
    default:
      {
	Int n = ArityOfFunctor(fun);
	Atom name  = NameOfFunctor(fun);
	ExpEntry *p;

	if (EndOfPAEntr(p = RepExpProp(_YAP_GetExpProp(name, n)))) {
	  Term ti[2];

	  /* error */
	  ti[0] = t;
	  ti[1] = MkIntegerTerm(n);
	  t = _YAP_MkApplTerm(_YAP_MkFunctor(_YAP_LookupAtom("/"),2), 2, ti);
	  _YAP_Error(TYPE_ERROR_EVALUABLE, t,
		"functor %s/%d for arithmetic expression",
		RepAtom(name)->StrOfAE,n);
	  P = (yamop *)FAILCODE;
	  RERROR();
	}
	if (n == 1)
	  return(p->FOfEE.unary(ArgOfTerm(1,t), USE_E_ARGS));
	return(p->FOfEE.binary(ArgOfTerm(1,t),ArgOfTerm(2,t), USE_E_ARGS));
      }
    }
  } else if (IsPairTerm(t)) {
    return(Eval(HeadOfTerm(t), USE_E_ARGS));
  } else if (IsIntTerm(t)) {
    RINT(IntOfTerm(t));
  } else {
    Atom name = AtomOfTerm(t);
    ExpEntry *p;

    if (EndOfPAEntr(p = RepExpProp(_YAP_GetExpProp(name, 0)))) {
      /* error */
      _YAP_Error(TYPE_ERROR_EVALUABLE, t,
	    "atom %s for arithmetic expression",
	    RepAtom(name)->StrOfAE);
      P = (yamop *)FAILCODE;
      RERROR();
    }
    return(p->FOfEE.constant(USE_E_ARGS));
  }
}

static Int 
p_is(void)
{				/* X is Y	 */
  union arith_ret res;
  blob_type bt;

  bt = Eval(Deref(ARG2), &res);
  return (_YAP_unify_constant(ARG1,EvalToTerm(bt,&res)));
}

void
_YAP_InitEval(void)
{
  /* here are the arithmetical predicates */
  _YAP_InitConstExps();
  _YAP_InitUnaryExps();
  _YAP_InitBinaryExps();
  _YAP_InitCPred("is", 2, p_is, TestPredFlag | SafePredFlag);
}

