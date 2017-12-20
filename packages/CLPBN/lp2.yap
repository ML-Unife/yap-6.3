/*
 cpl2pfl_translator 	module
 Program that translate a 'cplint file' in a 'pfl file'
 Author Zese Riccardo
 */

:- module(cpl2pfl_translator,[
	translate/2,
	translate/1,
	op(550,yfx,::)
	%constraints/1
	%op(1150,fx,constraints)
	]).

:- use_module(library(lists)).
:- use_module(library(pfl)).
%:- [-'create_table.pl'].

:-set_prolog_flag(single_var_warnings,on).



% Cplint default settings
setting(ground_body,true).		% if not true, the correct behavior is not garantee
setting(debug_mode,false).		% if true, the translator will write on standard input



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%													%
%		TABLE CREATION										%
%													%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*
 * Creates table for OR parfactor
 */
build_or_table(Terms,P,T):-
	create_table(P,Terms,T,'OR').

/*
 * Creates table for AND parfactor
 */	
build_and_table(Terms,P,T):-
	create_table(P,Terms,T,'AND').

/*
 * Create table for general parfactor
 */
create_table(P,Terms,[NP,P],_):-
	length(Terms,1),!,
	NP is 1 - P.
	
create_table(P,[\+_|L],Table,Mode):-
	!,
	length(L,Length),
	NC is 2^(Length),
	compute_column_value(L,Val),
	NP is 1 - P,
	build_row(NC,Val,[P,NP],Row1,Mode),
	build_row(NC,Val,[NP,P],Row2,Mode),
	append(Row1,Row2,Table),
	!.
	
create_table(P,[_H|L],Table,Mode):-
	length(L,Length),
	NC is 2^(Length),
	compute_column_value(L,Val),
	NP is 1.0 - P,
	build_row(NC,Val,[NP,P],Row1,Mode),
	build_row(NC,Val,[P,NP],Row2,Mode),
	append(Row1,Row2,Table),
	!.
	
/*
 * Create a row for the table
 */
build_row(0,_,_,[],_):-!.

build_row(NC,Val,[P,NP],[Prob|Valore1],Mode):-
	compute_value(Val,Valore,Mode),
	(Valore == 't' -> Prob = P ; Prob = NP),
	next_column(Val,Val1),
	N is NC - 1,
	build_row(N,Val1,[P,NP],Valore1,Mode).
	
/*
 * Computes the truth value of a column label
 */	
compute_value(['t'],'t',_).
compute_value(['f'],'f',_).
compute_value([H|T],Val,Mode):-
	compute_value(T,V,Mode),
	(Mode == 'AND' -> and(H,V,Val) ; or(H,V,Val)).
	
/*
 * compute AND and OR
 */	
and('t','t','t').
and('f','t','f').
and('t','f','f').
and('f','f','f').

or('t','t','t').
or('f','t','t').
or('t','f','t').
or('f','f','f').

/*
 * Compute the next column by updating true/false value of the column label
 */	
next_column(Val,Val1):-
	next_column(Val,Val1,_).

next_column(['t'],['f'],1).

next_column(['f'],['t'],1).

next_column(['t'|T],[V1,V|T1],C1):-
	next_column(T,[V|T1],C),
	(C==1 ->
		( V=='t' -> 
			V1='t', C1=0 
	 	  ; 
	 		V1='f', C1=1)
	 ;
	 	V1='t', C1=0
	 ).

next_column(['f'|T],[V1,V|T1],C1):-
	next_column(T,[V|T1],C),
	(C==1 ->
		( V=='t' -> 
			V1='f', C1=0 
	 	  ; 
	 		V1='t', C1=1)
	 ;
	 	V1='f', C1=0
	 ).
	
	
	
/*
 * Initialize the label of the first column
 */	
compute_column_value([],[]).

compute_column_value([\+_|T],['t'|T1]):-
	!,
	compute_column_value(T,T1).
	
compute_column_value([_|T],['f'|T1]):-
	compute_column_value(T,T1).
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	


% :- constraints(L) management
check_contraints([H/N]):-
	functor(A,H,N),
	assert(cons(A)).

check_contraints([H/N|T]):-
	functor(A,H,N),
	assert(cons(A)),
	check_contraints(T).


% creation of lists of terms for type TermsH;Table;TermsB
build_terms([H|TermsH],TermsB,H,B):-
	build_terms(TermsH,TermsB,B).

build_terms(TermsH,[H|TermsB],(H,T)):-
	cons(H),!,
	build_terms(TermsH,TermsB,T).
	
build_terms([H|TermsH],TermsB,(H,T)):-
	build_terms(TermsH,TermsB,T).

build_terms([],[T],T):-
	cons(T),!.	
	
build_terms([T],[],T).	
	

 
% given a list of terms return the list of all variables that occur in the terms
retrieve_term_list_variables((Term,L),Vars):-
	term_variables(Term,VarsT),
	retrieve_term_list_variables(L,VarsL),
	append(VarsT,VarsL,VarsDup),
	remove_duplicates(VarsDup,Vars).

retrieve_term_list_variables(Term,Vars):-
	term_variables(Term,Vars).
	



%%%%%%%%%%%%%%%%%%%%%%

% check for deciding if a term should translate in a noisy or operand


not_same(Var,VarsH):-
	\+ contains_var(Var,VarsH).
	
contains_var(Var,VarsH):-
	member(Var2,VarsH),
	Var == Var2.

is_noisy_or(H,B):-
	term_variables(H,VarsH),
	retrieve_term_list_variables(B,VarsB),!,
	member(Var,VarsB),
	not_same(Var,VarsH),!.
	
%%%%%%%%%%%%%%%%%%%%%%


% creation of a new term for deputy and or parfactor
create_intermediate_term(H,B,C):-
	ch_n(H,N),!,
	retract(ch_n(H,N)),
	N1 is N + 1,
	assert(ch_n(H,N1)),
	create_n_intermediate_term(H,B,C,N).
	
create_intermediate_term(H,B,C):-
	assert(ch_n(H,2)),
	create_n_intermediate_term(H,B,C,1).


create_intermediate_term(H,B):-
	create_intermediate_term(H,B,[112]).

create_n_intermediate_term(H,B,C,N):-
	functor(H,F,_),atom_codes(F,LF),
	atom_codes(N,L),append([LF,L,C],Chars),name(Name,Chars),
	term_variables(H,Vars),
	flatten([Name,Vars],Interm_term),
	B=..Interm_term.
	

% translate a list of predicates in a comma-separated set of predicates
list2term([End],End):-
	length([End],1).
	
list2term([H|TL],(H,TT)):-
	list2term(TL,TT).
	

% update the number of terms that share the same head
check_head(H):-
	functor(H,F,A),
	head_usg(F,A,N),!,
	retract(head_usg(F,A,N)),
	N1 is N + 1,
	assert(head_usg(F,A,N1)).
	
check_head(H):-
	functor(H,F,A),
	assert(head_usg(F,A,1)).
	

% prepare the list of terms to put in the goals list of a parfactor
check_deputy_termsB(TermsB,H,TermB0):-
	term_variables(H,Vars),
	check_deputy_variable_term(TermsB,Vars,TermB0).
	

check_deputy_variable_term(_,[],[]):-!.

check_deputy_variable_term([],_,[]).

check_deputy_variable_term([H|T],Vars,[H|T1]):-
	term_variables(H,VarsH),
	member(Var1,Vars),
	member(Var2,VarsH),
	Var1 == Var2, !,
	check_deputy_variable_term(T,Vars,T1).

check_deputy_variable_term([_H|T],Vars,T1):-
	check_deputy_variable_term(T,Vars,T1).



% substitute the negated terms with positive terms
make_all_positive([],[]).

make_all_positive([\+A|T],[A|T1]):-
 	make_all_positive(T,T1).
	
make_all_positive([A|T],[A|T1]):-
	make_all_positive(T,T1).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Management of the terms in the cpl file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% managements of the contraints
update_clause((:-constraints(L)),[]):-
	!,
	check_contraints(L).

% non ground probabilistic fact	
update_clause((H:P:-B),L):-
	float(P),!,
	%write(H),write(':'),write(P),write(':-'),write(B),nl,		%debug
	check_head(H),
	build_terms(TermsH,_TermsB,H,(B)),
	( length(TermsH,1) ->
			( build_and_table(TermsH,P,Table),
				L = (bayes H;Table;[B]));
			( write('Error in clause '),write(H),write(':'),write(P),write(':-'),write(B),nl,
				write('Non-contraint in the body. I quit!'),nl,
				abort) ).
	%write('	 -> '),write(L),nl,				%debug
	


% non ground probabilistic fact	
update_clause((P::H:-B),L):-
	float(P),!,
	%write(P),write('::'),write(H),write(':-'),write(B),nl,		%debug
	check_head(H),
	build_terms(TermsH,_TermsB,H,(B)),
	( length(TermsH,1) ->
			( build_and_table(TermsH,P,Table),
				L = (bayes H;Table;[B]));
			( write('Error in clause '),write(H),write(':'),write(P),write(':-'),write(B),nl,
				write('Non-contraint in the body. I quit!'),nl,
				abort) ).
	%write('	 -> '),write(L),nl,				%debug


% ground probabilistic fact
update_clause((H:P),(bayes H;Table;[])):-
	!,
	check_head(H),
	build_and_table([H],P,Table).


% ground probabilistic fact
update_clause((P::H),(bayes H;Table;[])):-
	!,
	check_head(H),
	build_and_table([H],P,Table).


% constraint rule
update_clause((H:-B),(H:-B)):-
	 cons(H),!.

% constraint fact
update_clause(H,H):-
	 cons(H),!.


% non-prob rule
update_clause((H:-B),L):-
	%write(H),write(':-'),write(B),nl,				%debug
	check_head(H),
	(is_noisy_or(H,B) ->
			 build_terms(TermsH,TermsB0,H,(B)),
			 build_and_table([H],1.0,Table),
			 create_intermediate_term(H,Interm_term),
			 delete(TermsH,H,Interm_termsH),
			 append([Interm_term],Interm_termsH,Ch_termsH),
			 build_and_table(Ch_termsH,1.0,Ch_table),
			 make_all_positive(Ch_termsH,Ch_termsH0),
			 list2term(Ch_termsH0,Ch_termsHT),
			 %list2term(TermsB0,TermsBT),
			 check_deputy_termsB(TermsB0,H,TermsBD0),
			 
			 make_all_positive(TermsB0,TermsB),
			 make_all_positive(TermsBD0,TermsBD),
			 L = [(deputy (H,Interm_term);TermsBD),(het Ch_termsHT;Ch_table;TermsB)]	
	 ;
			 build_terms(TermsH,TermsB0,H,(B)),
			 build_and_table(TermsH,1.0,Table),
			 make_all_positive(TermsH,TermsH0),
			 list2term(TermsH0,TermsHT),
			 
			 make_all_positive(TermsB0,TermsB),
			 %list2term(TermsB,TermsBT),
			 L = (bayes TermsHT;Table;TermsB)
	).
	
	%write('	 -> '),write(L),nl,				%debug
	
	
update_clause(end_of_file,[]):-!.
	
update_clause(H,H):-
	%write('H '),write(H),nl,					%debug
	check_head(H).
	
/*
% non-prob fact
update_clause(H,Factors):-
	H\= (:- check_contraints(_)),
	!,
	....
*/


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



/*
 Translate
 Argument: 
	FileName is the 'cplint file' name
	Note: The output file name will be FileName concatenated with '_2pfl'
*/
translate(FileName):-
	name(FileName,FNL),
	delete_extension(FNL,FN),
	name(FileNameWE,FN),
	write(FileNameWE),
	atom_concat(FileNameWE,'.pfl',FileNameOut),
	translate(FileName,FileNameOut).

/*
 Translate
 Argument: 
	FileName is the 'cplint file' name
	FileNameOut is the output file name
*/
translate(FileName,FileNameOut):-	

	%%%%%%%%%%%%%%%%%%
	% Startup system %
	%%%%%%%%%%%%%%%%%%

	% clean db
	retractall(ch_n(_,_)),
	retractall(cons(_)),
	retractall(vars_in_head(_,_)),
	retractall(head_usg(_,_,_)),
	assert(ch_n('',1)),
	assert(cons('')),
	assert(term_list([])),
	assert(head_usg('',0,0)),
	
	
	(setting(debug_mode,true) ->
		%atom_concat(FileName,'.cpl',FileNameCpl),
		open(FileName,read,StreamCpl),
		% read clauses file *.cpl
		read_clauses_cpl(StreamCpl,ClausesCpl,VarsName0),
		close(StreamCpl),
		write('ClausesCpl: '),nl,write(ClausesCpl),nl,			%debug
		update_clauses(ClausesCpl,NewClausesCpl0),
		write('NewClausesCpl0: '),nl,write(NewClausesCpl0),nl,		%debug
		look_for_or(NewClausesCpl0,NewClausesCpl,VarsName0,VarsName),
		write('NewClausesCpl: '),nl,write(NewClausesCpl),nl,		%debug
		testwritel(FileNameOut, NewClausesCpl,VarsName)
	 ;
	 	% create translation's files	
		create_translations_file(FileNameOut),
		%atom_concat(FileName,'.cpl',FileNameCpl),
		open(FileName,read,StreamCpl),
		% read clauses file *.cpl
		read_clauses_cpl(StreamCpl,ClausesCpl,VarsName0),
		close(StreamCpl),
		%write('ClausesCpl: '),nl,write(ClausesCpl),nl,			%debug
		update_clauses(ClausesCpl,NewClausesCpl0),
		%write('NewClausesCpl0: '),nl,write(NewClausesCpl0),nl,		%debug
		look_for_or(NewClausesCpl0,NewClausesCpl,VarsName0,VarsName),
		%write('NewClausesCpl: '),nl,write(NewClausesCpl),nl,		%debug
		writel(FileNameOut, NewClausesCpl,VarsName)
	).



/*
 * remove the extension from the file name
 */
delete_extension([],[]).

delete_extension([46|_],[]):- !.

delete_extension([H|T],[H|T1]):-
	delete_extension(T,T1).



/*
 * Creates and initializes the output file.
 */	
create_translations_file(FileName):-
	%atom_concat(FileName,'.pfl',RulesFile),
	open(FileName,write,S),
	current_output(O),
	set_output(S),	
	write('% solver settings'),nl,
	write(':- use_module(library(pfl)).'),nl,
	write(':- set_solver(lve).'),nl,nl,nl,
	write('% flags settings'),nl,
	write(':- yap_flag(unknown,fail).'),nl,
	write(':- set_pfl_flag(use_logarithms,false).'),nl,nl,nl,
	write('% set verbosity level'),nl,
	write('%:- set_pfl_flag(verbosity,3).'),nl,nl,nl,
	set_output(O),
	close(S).
	
/*
 * Transform clauses of cpl file into pfl parfactors.
 */	
update_clauses([],[]).

update_clauses([C|T],T1N):-
	update_clause(C,NC),
	(is_list(NC) ->
		!,
		update_clauses(T,T1),
		append(NC,T1,T1N)
	;
		(NC == [] -> 
			!,
			update_clauses(T,T1N)
		;
			!,
			update_clauses(T,T1),
			append([NC],T1,T1N)
		)
	). 


/*
 * Searches and updates parfactors in order to make OR parfactor.
 */	
look_for_or(NewClausesCpl0,NewClausesCpl,NameVars0,NameVars):-
	retractall(ch_n(_,_)),
	%assert(ch_n('',0)),
	update_or(NewClausesCpl0,NewClausesCpl1),
	retract(head_usg('',0,0)),
	findall(head_usg(F,A,N),head_usg(F,A,N),OrList),
	insert_new_or(OrList,OrToAdd,NameVars0,NameVars),
	reorder_clauses(NewClausesCpl1,OrToAdd,NewClausesCpl).
	
/*
 * Create OR parfactors.
 */	
insert_new_or([],[],VarsName,VarsName).

insert_new_or([head_usg(_F,_A,1)|T],T1,VarsName0,VarsName):- 
	!,
	insert_new_or(T,T1,VarsName0,VarsName).
	
insert_new_or([head_usg(F,A,N)|T],[(bayes TermsHT;Table;B)|T1],VarsNameOld,[Vars|VarsName]):-
	functor(H,F,A),
	create_all_or_intermediate_term(H,N,OrTerms),
	append([H],OrTerms,TermsH),
	build_or_table(TermsH,1.0,Table),
	list2term(TermsH,TermsHT),
	insert_new_or(T,T1,VarsNameOld,VarsName),
	b_for_h(H,B),
	vars_in_head(H,Vars).

	
/*
 * Create terms for a new OR parfactor.
 */	
create_all_or_intermediate_term(H,N,OrTerms):-
	create_all_or_intermediate_term(H,1,N,OrTerms).

create_all_or_intermediate_term(H,Lim,Lim,[H1]):-
	create_n_intermediate_term(H,H1,[],Lim).
	
create_all_or_intermediate_term(H,N,Lim,[H1|T1]):-
	create_n_intermediate_term(H,H1,[],N),
	N1 is N + 1,
	create_all_or_intermediate_term(H,N1,Lim,T1).


/*
 * Sort the parfactor, merging the initial parfactor list with the list containing 
 * the parfactors which model OR between the parfactors of the first list.
 */	
reorder_clauses(H,[],H).

reorder_clauses([bayes(((H,T1);BT))|T],[bayes(((T0,H,T2);BT2))|OrT],
		[bayes(((T0,H,T2);BT2)),bayes(((H,T1);BT))|ResT]):-
	!,
	%write('bayes'),nl,write(bayes(((H,T1);BT))),nl,write(bayes(((T0,H,T2);BT2))),nl,nl,	%debug
	reorder_clauses(T,OrT,ResT).
reorder_clauses([deputy(((H,T1);BT))|T],[bayes(((T0,H,T2);BT2))|OrT],
		[bayes(((T0,H,T2);BT2)),deputy(((H,T1);BT))|ResT]):-
	!,
	%write('deputy'),nl,write(deputy(((H,T1);BT))),nl,write(bayes(((T0,H,T2);BT2))),nl,nl,	%debug
	reorder_clauses(T,OrT,ResT).
reorder_clauses([H|T],OrT,[H|ResT]):-
	reorder_clauses(T,OrT,ResT).	

/*
 * Look for parfactor with same head and rename the heads for making a parfactor
 * that models the OR between the parfactors with same head.
 */	
update_or([],[]).

update_or([bayes(((H,T1);Tbl;T2))| T],[bayes(((H1,T1);Tbl;T2))|TT]):-
	functor(H,F,N),
	head_usg(F,N,Usg),
	assert(b_for_h(H,T2)),
	Usg @> 1,!,
	create_intermediate_term(H,H1,[]),
	update_or(T,TT).
	
update_or([deputy(((H,T1);T2))| T],[deputy(((H1,T1);T2))|TT]):-
	functor(H,F,N),
	head_usg(F,N,Usg),
	assert(b_for_h(H,T2)),
	Usg @> 1,!,
	%retract(ch_n(H,_)),
	create_intermediate_term(H,H1,[]),
	update_or(T,TT).
	
update_or([H| T],[H|TT]):-
	update_or(T,TT).




/*
 * Write parfactors on standard output.
 */		
testwritel(_,[],_).

testwritel(FileNameOut,[bayes((H;Tbl;B))|T],VarsName):-
	!,	
	write('bayes '),writeTerm(H,VarsName,UsedVars),write(';'),write(Tbl),write(';'),writeTermList(B,VarsName,UsedVars),write('.'),nl,
	testwritel(FileNameOut,T,VarsName).
	
testwritel(FileNameOut,[het((H;Tbl;B))|T],VarsName):-
	!,
	write('het '),writeTerm(H,VarsName,UsedVars),write(';'),write(Tbl),write(';'),writeTermList(B,VarsName,UsedVars),write('.'),nl,
	testwritel(FileNameOut,T,VarsName).
	
testwritel(FileNameOut,[deputy((H;B))|T],VarsName):-
	!,	
	write('deputy '),writeTerm(H,VarsName,UsedVars),write(';'),writeTermList(B,VarsName,UsedVars),write('.'),nl,
	testwritel(FileNameOut,T,VarsName).
	
testwritel(FileNameOut,[H|T],VarsName):-
	!,	
	write(H),write('.'),nl,
	testwritel(FileNameOut,T,VarsName).
	
	
/*
 * Write parfactors in file.
 */	
writel(_,[],_).

writel(FileNameOut,[bayes((H;Tbl;B))|T],VarsName):-
	!,	
	%atom_concat(FileNameOut,'.pfl',File),
	open(FileNameOut,append,S),
	current_output(O),
	set_output(S),
	write('bayes '),writeTerm(H,VarsName,UsedVars),write(';'),write(Tbl),write(';'),writeTermList(B,VarsName,UsedVars),write('.'),nl,
	set_output(O),
	close(S),
	writel(FileNameOut,T,VarsName).
	
writel(FileNameOut,[het((H;Tbl;B))|T],VarsName):-
	!,	
	%atom_concat(FileNameOut,'.pfl',File),
	open(FileNameOut,append,S),
	current_output(O),
	set_output(S),
	write('het '),writeTerm(H,VarsName,UsedVars),write(';'),write(Tbl),write(';'),writeTermList(B,VarsName,UsedVars),write('.'),nl,
	set_output(O),
	close(S),
	writel(FileNameOut,T,VarsName).
	
writel(FileNameOut,[deputy((H;B))|T],VarsName):-
	!,	
	%atom_concat(FileNameOut,'.pfl',File),
	open(FileNameOut,append,S),
	current_output(O),
	set_output(S),
	write('deputy '),writeTerm(H,VarsName,UsedVars),write(';'),writeTermList(B,VarsName,UsedVars),write('.'),nl,
	set_output(O),
	close(S),
	writel(FileNameOut,T,VarsName).
	
writel(FileNameOut,[H|T],VarsName):-
	!,	
	%atom_concat(FileNameOut,'.pfl',File),
	open(FileNameOut,append,S),
	current_output(O),
	set_output(S),
	write(H),write('.'),nl,
	set_output(O),
	close(S),
	writel(FileNameOut,T,VarsName).
	
	
	
	
/*
 * Write the list containing goal terms
 */	
writeTermList(Terms,VarsName,UsedVars):-
	write('['),
	writeTermList_int(Terms,VarsName,UsedVars),
	write(']').

writeTermList_int([],_,_).

writeTermList_int([H|T],VarsName,UsedVars):-
	length(T,L), L == 0,
	writeTermWithUsedVars(H,VarsName,UsedVars).

writeTermList_int([H|T],VarsName,UsedVars):-
	writeTermWithUsedVars(H,VarsName,UsedVars),
	write(','),
	writeTermList_int(T,VarsName,UsedVars).
	

/*
 * Write the goal terms.
 */
writeTermWithUsedVars((H,T),VarsName,UsedVars):-
	functor(H,F,_),
	term_variables(H,Vars),
	write(F),
	( length(Vars,L), L@>0 ->
		write('('),
		writeVariablesWithUsedVars(Vars,VarsName,UsedVars),
		write('),')
		;
			write(',')
	),
	writeTermWithUsedVars(T,VarsName,UsedVars).
	
writeTermWithUsedVars((H),VarsName,UsedVars):-
	functor(H,F,_),
	term_variables(H,Vars),
	length(Vars,L),
	L@>0,!,
	write(F),
	write('('),
	writeVariablesWithUsedVars(Vars,VarsName,UsedVars),
	write(')').
	
writeTermWithUsedVars((H),_,_):-
	functor(H,F,_),
	write(F).


/*
 * Write the variables of the goal terms
 */
writeVariablesWithUsedVars([],_,[]).

writeVariablesWithUsedVars([Var|T],VarsName,UsedVars):-
	length(T,0),!,
	( member(VarR,UsedVars), VarR == Var ->
		var2name(N,Var,VarsName),
		write(N)
		;
			write('_')
	).
writeVariablesWithUsedVars([Var|T],VarsName,UsedVars):-
	( member(VarR,UsedVars), VarR == Var ->
		var2name(N,Var,VarsName),
		write(N)
		;
			write('_')
	),
	write(','),
	writeVariablesWithUsedVars(T,VarsName,UsedVars).

/*
 * Write terms' parfactors.
 */
writeTerm((H,T),VarsName,UsedVars):-
	functor(H,F,_),
	term_variables(H,Vars),
	write(F),
	( length(Vars,L), L@>0 ->
		write('('),
		writeVariables(Vars,VarsName,UsedVars0),
		write('),')
		;
			write(','),
			UsedVars0 = []
	),
	writeTerm(T,VarsName,UsedVars1),
	flatten([UsedVars0|UsedVars1],UsedVars0T),
	remove_duplicates(UsedVars0T,UsedVars).
	
writeTerm((H),VarsName,UsedVars):-
	functor(H,F,_),
	term_variables(H,Vars),
	length(Vars,L),
	L@>0,!,
	write(F),
	write('('),
	writeVariables(Vars,VarsName,UsedVars),
	write(')').

writeTerm((H),_,[]):-
	functor(H,F,_),
	write(F).
	

/*
 * Write the variable's name
 */	
writeVariables([],_,[]).

writeVariables([Var|T],VarsName,[Var]):-
	length(T,0),!,
	var2name(N,Var,VarsName),
	write(N).
	
writeVariables([Var|T],VarsName,[Var|UsedVars]):-
	var2name(N,Var,VarsName),
	write(N),write(','),
	writeVariables(T,VarsName,UsedVars).
		
/*
 * Return te name of the variable in the parfactor. 
 */
var2name(Var,Var,[]):-
	write('Problem during variable translation!').
	

var2name(N,Var,[H|_]):-
	member((N=V),H),
	V==Var,!.
	
var2name(N,Var,[_|VarsName]):-
	var2name(N,Var,VarsName).
	
	
	



/*
 * Read clauses in a 'cplint file' 
 */
read_clauses_cpl(S,Clauses,VariablesName):-
	(setting(ground_body,true)->
		read_clauses_ground_body(S,Clauses,VariablesName)
	;
		read_clauses_exist_body(S,Clauses,VariablesName)
	).


/*
 * type of Reading clauses cplint file
 */
read_clauses_ground_body(S,[Cl|Out],[VN|Out2]):-
	read_term(S,Cl,[variable_names(VN)]),
	(Cl=end_of_file->
		Out=[],Out2=[],VN=[]
	;
		read_clauses_ground_body(S,Out,Out2),
		assert_cl_variables_head(Cl,VN)
	).


read_clauses_exist_body(S,[Cl|Out],[VN|Out2]):-
	read_term(S,Cl,[variable_names(VN)]),
	(Cl=end_of_file->
		Out=[],Out2=[],VN=[]
	;
		read_clauses_exist_body(S,Out,Out2)
	).


assert_cl_variables_head((H:_:-_),VN):-
	!,
	%write(vars_in_head(H,VN)),		%debug
	assert(vars_in_head(H,VN)).

assert_cl_variables_head((_::H:-_),VN):-
	!,
	%write(vars_in_head(H,VN)),		%debug
	assert(vars_in_head(H,VN)).

assert_cl_variables_head((H:-_),VN):-
	!,
	%write(vars_in_head(H,VN)),		%debug
	assert(vars_in_head(H,VN)).
	
assert_cl_variables_head(_,_).

/*
 * control if a list contains a predicate.
 */
member_eq(A,[H|_T]):-
	A==H.
	
member_eq(A,[_H|T]):-
	member_eq(A,T).
