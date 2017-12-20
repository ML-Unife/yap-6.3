%
% This module defines PFL, the Prolog Factor Language.
%
%

:- module(pfl,
		[op(550,yfx,@),
		 op(550,yfx,::),
		 op(1150,fx,bayes),
		 op(1150,fx,het),
		 op(1150,fx,deputy),
		 op(1150,fx,markov),
		 factor/6,
		 skolem/2,
		 defined_in_factor/2,
		 defined_in_factor/3,
		 evidence/2,
		 get_pfl_cpt/5, % given id and keys,  return new keys and cpt
		 get_pfl_parameters/3, % given id return par factor parameter
		 new_pfl_parameters/3, % given id  set new parameters
		 get_first_pvariable/2, % given id get firt pvar (useful in bayesian)
		 get_factor_pvariable/2, % given id get any pvar
		 add_ground_factor/5    %add a new bayesian variable (for now)
		]).

:- reexport(library(clpbn),
		[clpbn_flag/2 as pfl_flag,
		 set_clpbn_flag/2 as set_pfl_flag,
		 set_solver/1,
		 set_em_solver/1,
		 conditional_probability/3,
		 pfl_init_solver/5,
		 pfl_run_solver/3
		]).

:- reexport(library(clpbn/aggregates),
		[avg_factors/5]).

:- reexport('clpbn/horus',
		[set_horus_flag/2]).

:- ( % if clp(bn) has done loading, we're top-level
	predicate_property(set_pfl_flag(_,_), imported_from(clpbn))
   ->
	% we're using factor language
	% set appropriate flag
	set_pfl_flag(use_factors,on)
   ;
	% we're within clp(bn), no need to do anything
	true
   ).

:- use_module(library(atts)).

:- use_module(library(lists),
		[nth0/3,
		 append/3,
		 member/2
		]).

:- dynamic factor/6, skolem_in/2, skolem/2, preprocess/3, evidence/2, id/1.

user:term_expansion( bayes((Formula ; Phi ; Constraints)), pfl:factor(bayes,Id,FList,FV,Phi,Constraints)) :-
	!,
	term_variables(Formula, FreeVars),
	FV =.. [''|FreeVars],
	new_id(Id),
	process_args(Formula, Id, 0, _, FList, []).
user:term_expansion( het((Formula ;  Phi;Constraints)), pfl:factor(het(0),Id,FList,FV,Phi,Constraints)) :-
	!,
	term_variables(Formula, FreeVars),
	FV =.. [''|FreeVars],
	new_id(Id),
	process_args(Formula, Id, 0, _, FList, []).
user:term_expansion( deputy((Formula ; Constraints)), pfl:factor(deputy,Id,FList,FV,[1.0, 0.0, 0.0, 1.0],Constraints)) :-
	!,
	term_variables(Formula, FreeVars),
	FV =.. [''|FreeVars],
	new_id(Id),
	process_args(Formula, Id, 0, _, FList, []).
user:term_expansion( markov((Formula ; Phi ; Constraints)), pfl:factor(markov,Id,FList,FV,Phi,Constraints)) :-
	!,
	term_variables(Formula, FreeVars),
	FV =.. [''|FreeVars],
	new_id(Id),
	process_args(Formula, Id, 0, _, FList, []).
user:term_expansion( Id@N, L ) :-
	atom(Id), number(N), !,
	N1 is N + 1,
	findall(G,generate_entity(1, N1, Id, G), L).
user:term_expansion( Goal, [] ) :-
	preprocess(Goal, Sk,Var), !,
	(ground(Goal) -> true ; throw(error('non ground evidence',Goal))),
%	prolog_load_context(module, M),
	assert(pfl:evidence(Sk,Var)).
user:term_expansion( Goal, [] ) :-
	skolem( Goal, Dom),
	( Dom == [f,t] -> true ; throw(error('evidence for what value?',Goal))),
	(ground(Goal) -> true ; throw(error('non ground evidence',Goal))),
%	prolog_load_context(module, M),
	assert(pfl:evidence(Goal,1)).

Id@N :-
	generate_entity(0, N, Id, G),
	assert_static(user:G),
	fail.
_Id@_N.

add_ground_factor(bayes, Domain, Vars, CPT, Id) :-
	Vars = [K|_],
	( skolem(K,_Domain) -> true ; assert(skolem(K, Domain)) ),
	new_id(Id),
	asserta(skolem_in(K, Id)),
	assert(factor(bayes, Id, Vars, [], CPT, [])).

skolem(_Id:Key,Dom) :- skolem(Key, Dom).

defined_in_factor(Key, Factor) :-
	skolem_in(Key, Id),
	factor(bayes, Id, [Key|FList], FV, Phi, Constraints), !,
	Factor = factor(bayes, Id, [Key|FList], FV, Phi, Constraints).
defined_in_factor(Key, Factor) :-
	skolem_in(Key, Id),
	factor(markov, Id, FList, FV, Phi, Constraints),
	member(Key, FList),
	Factor = factor(markov, Id, FList, FV, Phi, Constraints).


defined_in_factor(Key, Id, 0) :-
	skolem_in(Key, Id),
	factor(bayes, Id, [Key|_FList], _FV, _Phi, _Constraints), !.
defined_in_factor(Key, Id, I) :-
	skolem_in(Key, Id),
	factor(markov, Id, FList, _FV, _Phi, _Constraints),
	nth0(I, FList, Key).


generate_entity(N, N, _, _) :- !.
generate_entity(I0, _N, Id, T) :-
	atomic_concat(p, I0, P),
	T =.. [Id, P].
generate_entity(I0, N, Id, T) :-
	I is I0+1,
	generate_entity(I, N, Id, T).

id(0).

new_id(Id) :-
	retract(id(Id0)),
	Id is Id0+1,
	assert(id(Id)).

process_args(V, _Id, _I0, _I ) --> { var(V) }, !,
	{ throw(error(instantiation_error,pfl:process_args)) }.
process_args((Arg1,V), Id, I0, I ) --> { var(V) }, !,
	{ I is I0+1 },
	process_arg(Arg1, Id, I),
	[V].
process_args((Arg1,Arg2), Id, I0, I ) --> !,
	process_args(Arg1, Id, I0, I1),
	process_args(Arg2, Id, I1, I).
process_args(Arg1, Id, I0, I ) -->
	{ I is I0+1 },
	process_arg(Arg1, Id, I).

process_arg(Sk::D, Id, _I) -->
	!,
	{
	  new_skolem(Sk,D),
	  assert(skolem_in(Sk, Id))
	},
	[Sk].
process_arg(or(In,Var,Sk), Id, _I) -->
	!,
	{
	  % if :: been used before for this skolem
	  % just keep on using it,
	  % otherwise, assume it is t,f
	  ( \+ \+ skolem(or(In,Var,Sk),_D) -> true ; new_skolem(or(In,Var,Sk),[none,one]) ),
	  assert(skolem_in(or(In,Var,Sk), Id))
	},
	[or(In,Var,Sk)].


process_arg(Sk, Id, _I) -->
	!,
	{
	  % if :: been used before for this skolem
	  % just keep on using it,
	  % otherwise, assume it is t,f
	  ( \+ \+ skolem(Sk,_D) -> true ; new_skolem(Sk,[f,t]) ),
	  assert(skolem_in(Sk, Id))
	},
	[Sk].

%
% redefinition
%
new_skolem(Sk, D) :-
	copy_term(Sk, Sk1),
	skolem(Sk1, D1),
	functor(Sk1, N, A),
	functor(Sk , N, A), !,
	( D1 = D -> true ; throw(pfl(permission_error(redefining_domain(Sk),D:D1)))).
%
%
% create interface and skolem descriptor
%
new_skolem(Sk, D) :-
	functor(Sk, N, A),
	functor(NSk, N, A),
	% [f,t] is special for evidence
	( D = [f,t] -> assert((evidence(NSk, 1) :- user:NSk)) ; true ),
	interface_predicate(NSk),
	assert(skolem(NSk, D)).

interface_predicate(Sk) :-
	Sk =.. SKAs,
	append(SKAs, [Var], ESKAs),
	ESk =.. ESKAs,
	assert(preprocess(ESk, Sk, Var)),
	% transform from PFL to CLP(BN) call
	assert_static((user:ESk :-
	  evidence(Sk,Ev) -> Ev = Var;
	  var(Var) -> insert_atts(Var,Sk) ;
	  add_evidence(Sk,Var)
	)).

insert_atts(Var,Sk) :-
	clpbn:put_atts(Var,[key(Sk)]).

add_evidence(Sk,Var) :-
	skolem(Sk,D),
	once(nth0(E,D,Var)),
	clpbn:put_atts(V,[key(Sk),evidence(E)]),
	( catch(b_getval(pfl_evidence, Vs), _, fail) ->
	    b_setval(pfl_evidence, [V|Vs])
	;
	    b_setval(pfl_evidence, [V])
	).


get_pfl_cpt(Id, Keys, Ev, NewKeys, Out) :-
 	factor(_Type,Id,[Key|_],_FV,avg,_Constraints), !,
 	Keys = [Key|Parents],
 	avg_factors(Key, Parents, 0.0, Ev, NewKeys, Out).
get_pfl_cpt(Id, Keys, _, Keys, Out) :-
 	factor(_Type,Id,Keys,_FV,Phi,_Constraints),
	( Phi = [_|_] -> Phi = Out ; call(user:Phi, Out) ).

get_pfl_parameters(Id, Keys, Out) :-
 	factor(_Type,Id,Keys,_FV,Phi,_Constraints),
	( Phi = [_|_] -> Phi = Out ; call(user:Phi, Out) ).


new_pfl_parameters(Id, Keys, NewPhi) :-
	retract(factor(Type,Id,Keys,FV,_Phi,Constraints)),
	assert(factor(Type,Id,Keys,FV,NewPhi,Constraints)),
	fail.
new_pfl_parameters(_Id, _Keys, _NewPhi).

get_pfl_factor_sizes(Id, DSizes) :-
	factor(_Type, Id, FList, _FV, _Phi, _Constraints),
	get_sizes(FList, DSizes).

get_sizes([], []).
get_sizes(Key.FList, Sz.DSizes) :-
	skolem(Key, Domain),
	length(Domain, Sz),
	get_sizes(FList, DSizes).

% only makes sense for bayesian networks
get_first_pvariable(Id,Var) :-
	factor(_Type, Id,Var._FList,_FV,_Phi,_Constraints).

% only makes sense for bayesian networks
get_factor_pvariable(Id,Var) :-
	factor(_Type, Id,FList,_FV,_Phi,_Constraints),
	member(Var, FList).

