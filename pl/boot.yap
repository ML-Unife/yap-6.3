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
* File:		boot.yap						 *
* Last rev:	8/2/88							 *
* mods:									 *
* comments:	boot file for Prolog					 *
*									 *
*************************************************************************/

% This one should come first so that disjunctions and long distance
% cuts are compiled right with co-routining.
%

true :- true.

'$live' :-
	'$init_system',
        '$do_live'.

'$do_live' :-
	repeat,
		'$current_module'(Module),
		( Module==user ->
		    '$compile_mode'(_,0)
		;
		    format(user_error,'[~w]~n', [Module])
		),
		'$system_catch'('$enter_top_level',Module,Error,user:'$Error'(Error)).

'$init_system' :-
	set_stream(user_input,alias('$loop_stream')),
        % do catch as early as possible
	(
	 '$access_yap_flags'(15, 0),
	 '$access_yap_flags'(22, 0),
	 \+ '$uncaught_throw'
	->
	  '$version'
	;
	  true
	),
	(
	 '$access_yap_flags'(22, 0) ->
	 set_value('$verbose',on)
	;
	 set_value('$verbose',off)
	),
	(
	 retractall(user:library_directory(_)),
	 '$system_library_directories'(D),
	 assertz(user:library_directory(D)),
	 fail
	;
	 true
	),
	'$enter_system_mode',
	'$init_globals',
	set_value(fileerrors,1),
	set_value('$gc',on),
	('$exit_undefp' -> true ; true),
	prompt1(' ?- '),
	'$debug_on'(false),
	% simple trick to find out if this is we are booting from Prolog.
	get_value('$user_module',V),
	(
	  V == []
	->
	  '$current_module'(_,prolog)
	  ;
	  '$current_module'(_,V), '$compile_mode'(_,0),
	  ('$access_yap_flags'(16,0) ->
	      ( exists('~/.yaprc') -> load_files('~/.yaprc', []) ; true ),
	      ( exists('~/.prologrc') -> load_files('~/.prologrc', []) ; true ),
	      ( exists('~/prolog.ini') -> load_files('~/prolog.ini', []) ; true )
	  ;
	      true
	  )
	),
	'$db_clean_queues'(0),
% this must be executed from C-code.
%	'$startup_saved_state',
	'$startup_reconsult',
	'$startup_goals',
	set_input(user_input),
	set_output(user_output),
	'$init_or_threads',
	'$run_at_thread_start'.


'$init_globals' :-
	'$init_consult',
	nb_setval('$chr_toplevel_show_store',false),
	nb_setval('$break',0),
	% '$set_read_error_handler'(error), let the user do that
	nb_setval('$open_expands_filename',true),
	nb_setval('$trace',off),
	nb_setval('$system_mode',off),
	nb_setval('$chr_toplevel_show_store',false),
	nb_setval('$assert_all',off),
	nb_setval('$if_skip_mode',no_skip),
	b_setval('$spy_glist',[]),
	nb_setval('$spy_gn',1),
	nb_setval('$debug_run',off),
	nb_setval('$debug_jump',off).

'$init_consult' :-
	set_value('$lf_verbose',informational),
	nb_setval('$if_level',0),
	nb_setval('$endif',off),
	nb_setval('$consulting_file',[]),
	nb_setval('$initialization_goals',off),
	nb_setval('$consulting',false),
	nb_setval('$included_file',[]).
	
'$init_or_threads' :-
	'$c_yapor_threads'(W), !,
	'$start_orp_threads'(W).
'$init_or_threads'.

'$start_orp_threads'(1) :- !.
'$start_orp_threads'(W) :-
	thread_create('$c_worker',_,[detached(true)]),
	W1 is W-1,
	'$start_orp_threads'(W1).


% Start file for yap

/*		I/O predicates						*/

/* meaning of flags for '$write' is
	  1	quote illegal atoms
	  2	ignore operator declarations
	  4	output '$VAR'(N) terms as A, B, C, ...
	  8	use portray(_)
*/

/* main execution loop							*/
'$read_vars'(user_input, Goal, Mod, Pos, Bindings) :-
	get_value('$readline',true), !,
	read_history(h, '!h',
                         [trace, end_of_file],
                         ' ?- ', Goal, Bindings),
	(nonvar(Err) ->
	 print_message(error,Err), fail
	;
	 true
	).
'$read_vars'(Stream,T,Mod,Pos,V) :-
	'$read'(true,T,Mod,V,Pos,Err,Stream),
	(nonvar(Err) ->
	 print_message(error,Err), fail
	;
	 true
	).

% reset alarms when entering top-level.
'$enter_top_level' :-
	'$alarm'(0, 0, _, _),
	fail.
'$enter_top_level' :-
	'$clean_up_dead_clauses',
	fail.
'$enter_top_level' :-
	'$nb_getval'('$break',BreakLevel,fail),
	 '$debug_on'(DBON),
	(
	 '$nb_getval'('$trace', on, fail)
	->
	 TraceDebug = trace
	;
	 DBON == true
	->
	 TraceDebug = debug
	;
	 true
	),
	print_message(informational,prompt(BreakLevel,TraceDebug)),
	fail.
'$enter_top_level' :-
	get_value('$top_level_goal',GA), GA \= [], !,
	set_value('$top_level_goal',[]),
	'$run_atom_goal'(GA),
	set_value('$live','$false').
'$enter_top_level' :-
	'$disable_docreep',
	prompt(_,'| '),
	prompt1(' ?- '),
	'$run_toplevel_hooks',
	'$read_vars'(user_input,Command,_,Pos,Varnames),
	nb_setval('$spy_gn',1),
		% stop at spy-points if debugging is on.

	nb_setval('$debug_run',off),
	nb_setval('$debug_jump',off),
	prompt(_,'|: '),
	'$command'(Command,Varnames,Pos,top),
	'$sync_mmapped_arrays',
	set_value('$live','$false').

%
% first, recover what we need from the saved state...
%
'$startup_saved_state' :-
	'$init_path_extensions',
	fail.
% use if we come from a save_program and we have SWI's shlib
'$startup_saved_state' :-
	recorded('$reload_foreign_libraries',G,R),
	erase(R),
	shlib:reload_foreign_libraries,
	fail.
% use if we come from a save_program and we have a goal to execute
'$startup_saved_state' :-
	recorded('$restore_goal',G,R),
	erase(R),
	prompt(_,'| '),
	'$system_catch'('$do_yes_no'((G->true),user),user,Error,user:'$Error'(Error)),
	fail.
'$startup_saved_state'.

% then recover program.
'$startup_reconsult' :-
	get_value('$consult_on_boot',X), X \= [], !,
	set_value('$consult_on_boot',[]),
	'$do_startup_reconsult'(X).
'$startup_reconsult'.


% then we can execute the programs.
'$startup_goals' :-
	recorded('$startup_goal',G,_),
	'$current_module'(Module),
	'$system_catch'('$query'(once(G), []),Module,Error,user:'$Error'(Error)),
	fail.
'$startup_goals' :-
	get_value('$init_goal',GA),
	GA \= [],
	set_value('$init_goal',[]),
	'$run_atom_goal'(GA),
	fail.
'$startup_goals' :-
	get_value('$myddas_goal',GA), GA \= [],
	set_value('$myddas_goal',[]),
	get_value('$myddas_user',User), User \= [],
	set_value('$myddas_user',[]),
	get_value('$myddas_db',Db), Db \= [],
	set_value('$myddas_db',[]),
	get_value('$myddas_host',HostT),
	( HostT \= [] ->
	  Host = HostT,
	  set_value('$myddas_host',[])
	;
	  Host = localhost
	),
	get_value('$myddas_pass',PassT),
	( PassT \= [] ->
	  Pass = PassT,
	  set_value('$myddas_pass',[])
	;
	  Pass = ''
	),
	use_module(library(myddas)),
	call(db_open(mysql,myddas,Host/Db,User,Pass)),
	'$myddas_import_all',
	fail.
'$startup_goals'.

 %
 % MYDDAS: Import all the tables from one database
 %

 '$myddas_import_all':-
	 call(db_my_show_tables(myddas,table(Table))),
	 call(db_import(myddas,Table,Table)),
	 fail.
 '$myddas_import_all'.
	 


 '$erase_sets' :- 
		 eraseall('$'),
		 eraseall('$$set'),
		 eraseall('$$one'), 
		 eraseall('$reconsulted'), fail.
 '$erase_sets' :- \+ recorded('$path',_,_), recorda('$path',"",_).
 '$erase_sets'.

 '$version' :- 
	 get_value('$version_name',VersionName),
	 print_message(help, version(VersionName)),
	 get_value('$myddas_version_name',MYDDASVersionName),
	 MYDDASVersionName \== [],
	 print_message(help, myddas_version(MYDDASVersionName)),
	 fail.
 '$version' :-
	 recorded('$version',VersionName,_),
	 print_message(help, VersionName),
	 fail.
 '$version'.

 repeat :- '$repeat'.

 '$repeat'.
 '$repeat'.
 '$repeat'.
 '$repeat'.
 '$repeat'.
 '$repeat'.
 '$repeat'.
 '$repeat'.
 '$repeat'.
 '$repeat' :- '$repeat'.

'$start_corouts' :-
	recorded('$corout','$corout'(Name,_,_),R),
	Name \= main,
	finish_corout(R),
	fail.
'$start_corouts' :- 
	eraseall('$corout'),
	eraseall('$result'),
	eraseall('$actual'),
	fail.
'$start_corouts' :- recorda('$actual',main,_),
	recordz('$corout','$corout'(main,main,'$corout'([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[])),_Ref),
	recorda('$result',going,_).

'$command'(C,VL,Pos,Con) :-
	'$access_yap_flags'(9,1), !,
	 '$execute_command'(C,VL,Pos,Con,C).
'$command'(C,VL,Pos,Con) :-
	( (Con = top ; var(C) ; C = [_|_])  ->  
	  '$execute_command'(C,VL,Pos,Con,C), ! ;
	  % do term expansion
	  expand_term(C, EC),
	  % execute a list of commands
	  '$execute_commands'(EC,VL,Pos,Con,C),
	  % succeed only if the *original* was at end of file.
	  C == end_of_file
	).

 %
 % Hack in case expand_term has created a list of commands.
 %
 '$execute_commands'(V,_,_,_,Source) :- var(V), !,
	 '$do_error'(instantiation_error,meta_call(Source)).
 '$execute_commands'([],_,_,_,_) :- !.
 '$execute_commands'([C|Cs],VL,Pos,Con,Source) :- !,
	 (
	   '$execute_command'(C,VL,Pos,Con,Source),
	   fail	
	 ;
	   '$execute_commands'(Cs,VL,Pos,Con,Source)
	 ).
 '$execute_commands'(C,VL,Pos,Con,Source) :-
	 '$execute_command'(C,VL,Pos,Con,Source).



				%
 %
 %

 '$execute_command'(C,_,_,top,Source) :- var(C), !,
	 '$do_error'(instantiation_error,meta_call(Source)).
 '$execute_command'(C,_,_,top,Source) :- number(C), !,
	 '$do_error'(type_error(callable,C),meta_call(Source)).
 '$execute_command'(R,_,_,top,Source) :- db_reference(R), !,
	 '$do_error'(type_error(callable,R),meta_call(Source)).
 '$execute_command'(end_of_file,_,_,_,_) :- !.
 '$execute_command'(Command,_,_,_,_) :-
	 '$nb_getval'('$if_skip_mode', skip, fail),
	 \+ '$if_directive'(Command),
	 !.
 '$execute_command'((:-G),_,_,Option,_) :-
%          !,
	 Option \= top, !,
	 '$current_module'(M),
	 % allow user expansion
	 expand_term((:- G), O),
         O = (:- G1),
	 '$process_directive'(G1, Option, M).
 '$execute_command'((?-G), V, Pos, Option, Source) :-
	 Option \= top, !,
	 '$execute_command'(G, V, Pos, top, Source).
 '$execute_command'(G, V, Pos, Option, Source) :-
	 '$continue_with_command'(Option, V, Pos, G, Source).

 %
 % This command is very different depending on the language mode we are in.
 %
 % ISO only wants directives in files
 % SICStus accepts everything in files
 % YAP accepts everything everywhere
 % 
 '$process_directive'(G, top, M) :-
	 '$access_yap_flags'(8, 0), !, % YAP mode, go in and do it,
	 '$process_directive'(G, consult, M).
 '$process_directive'(G, top, _) :- !,
	 '$do_error'(context_error((:- G),clause),query).
 %
 % allow modules
 %
 '$process_directive'(M:G, Mode, _) :- !,
	 '$process_directive'(G, Mode, M).
 %
 % default case
 %
 '$process_directive'(Gs, Mode, M) :-
	 '$all_directives'(Gs), !,
	 '$exec_directives'(Gs, Mode, M).

 %
 % ISO does not allow goals (use initialization).
 %
 '$process_directive'(D, _, M) :-
	 '$access_yap_flags'(8, 1), !, % ISO Prolog mode, go in and do it,
	 '$do_error'(context_error((:- M:D),query),directive).
 %
 % but YAP and SICStus does.
 %
 '$process_directive'(G, _, M) :-
	 '$exit_system_mode',
	 ( '$notrace'(M:G) -> true ; format(user_error,':- ~w:~w failed.~n',[M,G]) ),
	 '$enter_system_mode'.

 '$continue_with_command'(Where,V,'$stream_position'(C,_P,A1,A2,A3),'$source_location'(_F,L):G,Source) :- !,
	  '$continue_with_command'(Where,V,'$stream_position'(C,L,A1,A2,A3),G,Source).
 '$continue_with_command'(reconsult,V,Pos,G,Source) :-
	 '$go_compile_clause'(G,V,Pos,5,Source),
	 fail.
 '$continue_with_command'(consult,V,Pos,G,Source) :-
	 '$go_compile_clause'(G,V,Pos,13,Source),
	 fail.
 '$continue_with_command'(top,V,_,G,_) :-
	 '$query'(G,V).

 %
 % not 100% compatible with SICStus Prolog, as SICStus Prolog would put
 % module prefixes all over the place, although unnecessarily so.
 %
 '$go_compile_clause'(G,V,Pos,N,Source) :-
	 '$current_module'(Mod),
	 '$go_compile_clause'(G,V,Pos,N,Mod,Mod,Source).
 
'$go_compile_clause'(G,_,_,_,_,_,Source) :-
	var(G), !,
	'$do_error'(instantiation_error,assert(Source)).	
'$go_compile_clause'((G:-_),_,_,_,_,_,Source) :-
	var(G), !,
	'$do_error'(instantiation_error,assert(Source)).	
'$go_compile_clause'(M:G,V,Pos,N,_,_,Source) :- !,
	  '$go_compile_clause'(G,V,Pos,N,M,M,Source).
'$go_compile_clause'((M:H :- B),V,Pos,N,_,BodyMod,Source) :- !,
	  '$go_compile_clause'((H :- B),V,Pos,N,M,BodyMod,Source).
'$go_compile_clause'(G,V,Pos,N,HeadMod,BodyMod,Source) :- !,
	 '$prepare_term'(G, V, Pos, G0, G1, BodyMod, HeadMod, Source),
	 '$$compile'(G1, G0, N, HeadMod).

 '$prepare_term'(G, V, Pos, G0, G1, BodyMod, SourceMod, Source) :-
	 (
	     get_value('$syntaxcheckflag',on)
          ->
	     '$check_term'(Source, V, Pos, BodyMod)
	 ;
	     true 
	 ),
	 '$precompile_term'(G, G0, G1, BodyMod, SourceMod).

 % process an input clause
 '$$compile'(G, G0, L, Mod) :-
	 '$head_and_body'(G,H,_),
	 '$flags'(H, Mod, Fl, Fl),
	 is(NFl, /\, Fl, 0x00002000),
	 (
	  NFl \= 0
	 ->
	  '$assertz_dynamic'(L,G,G0,Mod)
	 ;
	  nb_getval('$assert_all',on)
	 ->
	  functor(H,N,A),
	  '$dynamic'(N/A,Mod),
	  '$assertz_dynamic'(L,G,G0,Mod)
	 ;
	  '$not_imported'(H, Mod),
	  '$compile'(G, L, G0, Mod)
	 ).

'$not_imported'(H, Mod) :-
	recorded('$import','$import'(NM,Mod,NH,H,_,_),_),
	NM \= Mod, !,
	functor(NH,N,Ar),
	'$do_error'(permission_error(modify, static_procedure, NM:N/Ar), consult).
'$not_imported'(_, _).


'$check_if_reconsulted'(N,A) :- 
         once(recorded('$reconsulted',N/A,_)), 
	 recorded('$reconsulted',X,_), 
	 ( X = N/A , !; 
	   X = '$', !, fail; 
	   fail 
	 ). 

'$inform_as_reconsulted'(N,A) :-
	 recorda('$reconsulted',N/A,_).

'$clear_reconsulting' :-
	recorded('$reconsulted',X,Ref),
	erase(Ref),
	X == '$', !,
	( recorded('$reconsulting',_,R) -> erase(R) ).

'$prompt_alternatives_on'(determinism).

/* Executing a query */

'$query'(end_of_file,_).

 % ***************************
 % * -------- YAPOR -------- *
 % ***************************

'$query'(G,V) :-
	 \+ '$undefined'('$c_yapor_on', prolog),
	 '$c_yapor_on',
	 \+ '$undefined'('$c_start_yapor', prolog),
	 '$parallelizable'(G), !,
	 '$parallel_query'(G,V),
	 fail.

 % end of YAPOR

'$query'(G,[]) :-
	 '$prompt_alternatives_on'(OPT),
	 ( OPT = groundness ; OPT = determinism), !,
	 '$yes_no'(G,(?-)).
'$query'(G,V) :-
	 (
	   '$exit_system_mode',
	  yap_hacks:current_choice_point(CP),
	  '$execute'(G),
	  yap_hacks:current_choice_point(NCP),
	  ( '$enter_system_mode' ; '$exit_system_mode', fail),
	  '$delayed_goals'(G, V, NV, LGs),
	  '$write_answer'(NV, LGs, Written),
	  '$write_query_answer_true'(Written),
	  (
	   '$prompt_alternatives_on'(determinism), CP = NCP ->
	   nl(user_error),
	   !
	  ;
	   '$another',
	   !
	  ),
	  fail	 
	 ;
	  '$enter_system_mode',
	  '$out_neg_answer'
	 ).

 '$yes_no'(G,C) :-
	 '$current_module'(M),
	 '$do_yes_no'(G,M),
	 '$delayed_goals'(G, [], NV, LGs),
	 '$write_answer'(NV, LGs, Written),
	 ( Written = [] ->
	 !,'$present_answer'(C, yes);
	 '$another', !
	 ),
	 fail.
 '$yes_no'(_,_) :-
	 '$enter_system_mode',
	 '$out_neg_answer'.

'$add_env_and_fail' :- fail.

'$delayed_goals'(G, V, NV, LGs) :-
	'$attributes':delayed_goals(G, V, NV, LGs), !.
'$delayed_goals'(_, V, NV, []) :-
	copy_term_nat(V, NV).

'$out_neg_answer' :-
	 ( '$undefined'(print_message(_,_),prolog) -> 
	    '$present_answer'(user_error,"no~n", [])
	 ;
	    print_message(help,no)
	 ),
	 fail.

'$do_yes_no'([X|L], M) :- !, '$csult'([X|L], M).
'$do_yes_no'(G, M) :-
	'$exit_system_mode',
	'$execute'(M:G),
	( '$enter_system_mode' ; '$exit_system_mode', fail).

'$write_query_answer_true'([]) :- !,
	format(user_error,'~ntrue',[]).
'$write_query_answer_true'(_).


%
% present_answer has three components. First it flushes the streams,
% then it presents the goals, and last it shows any goals frozen on
% the arguments.
%
'$present_answer'(_,_):-
        flush_output,
	fail.
'$present_answer'((?-), Answ) :-
	nb_getval('$break',BL),
	( BL \= 0 -> 	format(user_error, '[~p] ',[BL]) ;
			true ),
        ( recorded('$print_options','$toplevel'(Opts),_) ->
	   write_term(user_error,Answ,Opts) ;
	   format(user_error,'~w',[Answ])
        ),
	format(user_error,'~n', []).

'$another' :-
	format(user_error,' ? ',[]),
	get0(user_input,C),
	'$do_another'(C).

'$do_another'(C) :-
	(   C== 0'; ->  skip(user_input,10), %'
	    '$add_nl_outside_console',
	    fail
	;
	    C== 10 -> '$add_nl_outside_console',
		( '$undefined'(print_message(_,_),prolog) -> 
			format(user_error,'yes~n', [])
	        ;
		   print_message(help,yes)
		)
	;
	    C== 13 -> 
	    get0(user_input,NC),
	    '$do_another'(NC)	    
	;
	    C== -1 -> halt
	;
	    skip(user_input,10), '$ask_again_for_another'
	).

%'$add_nl_outside_console' :-
%	'$is_same_tty'(user_input, user_error), !.
'$add_nl_outside_console' :-
	format(user_error,'~n',[]).

'$ask_again_for_another' :-
	format(user_error,'Action (\";\" for more choices, <return> for exit)', []),
	'$another'.

'$write_answer'(_,_,_) :-
        flush_output,
	fail.
'$write_answer'(Vs, LBlk, FLAnsw) :-
	'$purge_dontcares'(Vs,IVs),
	'$sort'(IVs, NVs),
	'$prep_answer_var_by_var'(NVs, LAnsw, LBlk),
	'$name_vars_in_goals'(LAnsw, Vs, NLAnsw),
        '$write_vars_and_goals'(NLAnsw, first, FLAnsw).

'$purge_dontcares'([],[]).
'$purge_dontcares'([[[95|_]|_]|Vs],NVs) :- !,
	'$purge_dontcares'(Vs,NVs).
'$purge_dontcares'([V|Vs],[V|NVs]) :-
	'$purge_dontcares'(Vs,NVs).


'$prep_answer_var_by_var'([], L, L).
'$prep_answer_var_by_var'([[Name|Value]|L], LF, L0) :- 
	'$delete_identical_answers'(L, Value, NL, Names),
	'$prep_answer_var'([Name|Names], Value, LF, LI),
	'$prep_answer_var_by_var'(NL, LI, L0).

% fetch all cases that have the same solution.
'$delete_identical_answers'([], _, [], []).
'$delete_identical_answers'([[Name|Value]|L], Value0, FL, [Name|Names]) :-
	Value == Value0, !,
	'$delete_identical_answers'(L, Value0, FL, Names).
'$delete_identical_answers'([VV|L], Value0, [VV|FL], Names) :-
	'$delete_identical_answers'(L, Value0, FL, Names).

% now create a list of pairs that will look like goals.
'$prep_answer_var'(Names, Value, LF, L0) :- var(Value), !,
	'$prep_answer_unbound_var'(Names, LF, L0).
'$prep_answer_var'(Names, Value, [nonvar(Names,Value)|L0], L0).

% ignore unbound variables
'$prep_answer_unbound_var'([_], L, L) :- !.
'$prep_answer_unbound_var'(Names, [var(Names)|L0], L0).

'$gen_name_string'(I,L,[C|L]) :- I < 26, !, C is I+65.
'$gen_name_string'(I,L0,LF) :-
	I1 is I mod 26,
	I2 is I // 26,
	C is I1+65,
	'$gen_name_string'(I2,[C|L0],LF).

'$write_vars_and_goals'([], _, []).
'$write_vars_and_goals'([nl,G1|LG], First, NG) :- !,
	nl(user_error),
	'$write_goal_output'(G1, First, NG, Next, IG),
	'$write_vars_and_goals'(LG, Next, IG).
'$write_vars_and_goals'([G1|LG], First, NG) :-
	'$write_goal_output'(G1, First, NG, Next, IG),
	'$write_vars_and_goals'(LG, Next, IG).

'$goal_to_string'(Format, G, String) :-
	format(codes(String),Format,G).

'$write_goal_output'(var([V|VL]), First, [var([V|VL])|L], next, L) :- !,
        ( First = first -> true ; format(user_error,',~n',[]) ),
	format(user_error,'~s',[V]),
	'$write_output_vars'(VL).
'$write_goal_output'(nonvar([V|VL],B), First, [nonvar([V|VL],B)|L], next, L) :- !,
        ( First = first -> true ; format(user_error,',~n',[]) ),
	format(user_error,'~s',[V]),
	'$write_output_vars'(VL),
	format(user_error,' = ', []),
        ( recorded('$print_options','$toplevel'(Opts),_) ->
	   write_term(user_error,B,[priority(699)|Opts]) ;
	   write_term(user_error,B,[priority(699)])
        ).
'$write_goal_output'(nl, First, NG, First, NG) :- !,
	format(user_error,'~n',[]).
'$write_goal_output'(Format-G, First, NG, Next, IG) :- !,
	G = [_|_], !,
	% dump on string first so that we can check whether we actually
	% had any output from the solver.
	'$goal_to_string'(Format, G, String),
	( String == [] ->
	    % we didn't
	    IG = NG, First = Next
	;
	    % we did
	    ( First = first -> true ; format(user_error,',~n',[]) ),
	    format(user_error, '~s', [String]),
	    NG = [G|IG]
	).
'$write_goal_output'(_-G, First, [G|NG], next, NG) :- !,
        ( First = first -> true ; format(user_error,',~n',[]) ),
        ( recorded('$print_options','$toplevel'(Opts),_) ->
	   write_term(user_error,G,Opts) ;
	   format(user_error,'~w',[G])
        ).
'$write_goal_output'(_M:G, First, [G|NG], next, NG) :- !,
        ( First = first -> true ; format(user_error,',~n',[]) ),
        ( recorded('$print_options','$toplevel'(Opts),_) ->
	   write_term(user_error,G,Opts) ;
	   format(user_error,'~w',[G])
        ).
'$write_goal_output'(G, First, [M:G|NG], next, NG) :-
	'$current_module'(M),
        ( First = first -> true ; format(user_error,',~n',[]) ),
        ( recorded('$print_options','$toplevel'(Opts),_) ->
	   write_term(user_error,G,Opts) ;
	   format(user_error,'~w',[G])
        ).

'$name_vars_in_goals'(G, VL0, G) :-
	'$name_well_known_vars'(VL0),
	'$variables_in_term'(G, [], GVL),
	'$name_vars_in_goals1'(GVL, 0, _).

'$name_well_known_vars'([]).
'$name_well_known_vars'([[Name|V]|NVL0]) :-
	var(V), !,
	V = '$VAR'(Name),
	'$name_well_known_vars'(NVL0).
'$name_well_known_vars'([_|NVL0]) :-
	'$name_well_known_vars'(NVL0).

'$name_vars_in_goals1'([], I, I).
'$name_vars_in_goals1'(['$VAR'([95|Name])|NGVL], I0, IF) :-
	I is I0+1,
	'$gen_name_string'(I0,[],Name), !,
	'$name_vars_in_goals1'(NGVL, I, IF).
'$name_vars_in_goals1'([NV|NGVL], I0, IF) :-
	nonvar(NV),
	'$name_vars_in_goals1'(NGVL, I0, IF).

'$write_output_vars'([]).
'$write_output_vars'([V|VL]) :-
	format(user_error,' = ~s',[V]),
	'$write_output_vars'(VL).

call(G) :- '$execute'(G).

incore(G) :- '$execute'(G).

%
% standard meta-call, called if $execute could not do everything.
%
'$meta_call'(G, M) :-
	yap_hacks:current_choice_point(CP),
	'$call'(G, CP, G, M).


','(X,Y) :-
	yap_hacks:env_choice_point(CP),
	'$current_module'(M),
        '$call'(X,CP,(X,Y),M),
        '$call'(Y,CP,(X,Y),M).
';'((X->A),Y) :- !,
	yap_hacks:env_choice_point(CP),
	'$current_module'(M),
        ( '$execute'(X)
	->
	  '$call'(A,CP,(X->A;Y),M)
	;
	  '$call'(Y,CP,(X->A;Y),M)
	).
';'((X*->A),Y) :- !,
	yap_hacks:env_choice_point(CP),
	'$current_module'(M),
	(
	 yap_hacks:current_choicepoint(DCP),
	 '$execute'(X),
	 yap_hacks:cut_at(DCP),
	 '$call'(A,CP,((X*->A),Y),M)
        ;
	 '$call'(Y,CP,((X*->A),Y),M)
	).
';'(X,Y) :-
	yap_hacks:env_choice_point(CP),
	'$current_module'(M),
        ( '$call'(X,CP,(X;Y),M) ; '$call'(Y,CP,(X;Y),M) ).
'|'(X,Y) :-
	yap_hacks:env_choice_point(CP),
	'$current_module'(M),
        ( '$call'(X,CP,(X|Y),M) ; '$call'(Y,CP,(X|Y),M) ).
'->'(X,Y) :-
	yap_hacks:env_choice_point(CP),
	'$current_module'(M),
        ( '$call'(X,CP,(X->Y),M) -> '$call'(Y,CP,(X->Y),M) ).
'*->'(X,Y) :-
	yap_hacks:env_choice_point(CP),
	'$current_module'(M),
        ( '$call'(X,CP,(X*->Y),M), '$call'(Y,CP,(X*->Y),M) ).
\+(G) :-     \+ '$execute'(G).
not(G) :-    \+ '$execute'(G).

'$cut_by'(CP) :- '$$cut_by'(CP).

%
% do it in ISO mode.
%
'$meta_call'(G,_ISO,M) :-
	'$iso_check_goal'(G,G),
	yap_hacks:current_choice_point(CP),
	'$call'(G, CP, G, M).

'$meta_call'(G, CP, G0, M) :-
	'$call'(G, CP, G0, M).

'$call'(G, CP, G0, _, M) :-  /* iso version */
	'$iso_check_goal'(G,G0),
	'$call'(G, CP, G0, M).


'$call'(M:_,_,G0,_) :- var(M), !,
	'$do_error'(instantiation_error,call(G0)).
'$call'(M:G,CP,G0,_) :- !,
        '$call'(G,CP,G0,M).
'$call'((X,Y),CP,G0,M) :- !,
        '$call'(X,CP,G0,M),
        '$call'(Y,CP,G0,M).
'$call'((X->Y),CP,G0,M) :- !,
	(
	 '$call'(X,CP,G0,M)
          ->
	 '$call'(Y,CP,G0,M)
	).
'$call'((X*->Y),CP,G0,M) :- !,
	'$call'(X,CP,G0,M),
	'$call'(Y,CP,G0,M).
'$call'((X->Y; Z),CP,G0,M) :- !,
	(
	    '$call'(X,CP,G0,M)
         ->
	    '$call'(Y,CP,G0,M)
        ;
	    '$call'(Z,CP,G0,M)
	).
'$call'((X*->Y; Z),CP,G0,M) :- !,
	(
	 yap_hacks:current_choicepoint(DCP),
	 '$call'(X,CP,G0,M),
	 yap_hacks:cut_at(DCP),
	 '$call'(Y,CP,G0,M)
        ;
	 '$call'(Z,CP,G0,M)
	).
'$call'((A;B),CP,G0,M) :- !,
	(
	    '$call'(A,CP,G0,M)
        ;
	    '$call'(B,CP,G0,M)
	).
'$call'((X->Y| Z),CP,G0,M) :- !,
	(
	    '$call'(X,CP,G0,M)
         ->
	 '$call'(Y,CP,G0,M)
        ;
	'$call'(Z,CP,G0,M)
	).
'$call'((X*->Y| Z),CP,G0,M) :- !,
	(
	 yap_hacks:current_choicepoint(DCP),
	 '$call'(X,CP,G0,M),
	 yap_hacks:cut_at(DCP),
	 '$call'(Y,CP,G0,M)
        ;
	 '$call'(Z,CP,G0,M)
	).
'$call'((A|B),CP, G0,M) :- !,
	(
	    '$call'(A,CP,G0,M)
        ;
	    '$call'(B,CP,G0,M)
	).
'$call'(\+ X, _CP, _G0, M) :- !,
	yap_hacks:current_choicepoint(CP),
	\+  '$call'(X,CP,G0,M).
'$call'(not(X), _CP, _G0, M) :- !,
	\+  '$call'(X,CP,G0,M).
'$call'(!, CP, _,_) :- !,
	'$$cut_by'(CP).
'$call'([A|B], _, _, M) :- !,
	'$csult'([A|B], M).
'$call'(G, CP, G0, CurMod) :-
	( '$is_expand_goal_or_meta_predicate'(G,CurMod) ->
	   (
	     '$notrace'(user:goal_expansion(G, CurMod, NG)) ->
	       '$call'(NG, CP, G0,CurMod)
	     ;
	       % repeat other code.
             '$is_metapredicate'(G,CurMod) ->
	       (
	         '$meta_expansion'(G,CurMod,CurMod,CurMod,NG,[]) ->
	         '$execute0'(NG, CurMod)
	       ;
	         '$execute0'(G, CurMod)
	       )
	   ;
	     '$execute0'(G, CurMod)
	   )
	;
	  '$execute0'(G, CurMod)
	).

'$check_callable'(V,G) :- var(V), !,
	'$do_error'(instantiation_error,G).
'$check_callable'(M:_G1,G) :- var(M), !,
	'$do_error'(instantiation_error,G).
'$check_callable'(_:G1,G) :- !,
	'$check_callable'(G1,G).
'$check_callable'(A,G) :- number(A), !,
	'$do_error'(type_error(callable,A),G).
'$check_callable'(R,G) :- db_reference(R), !,
	'$do_error'(type_error(callable,R),G).
'$check_callable'(_,_).

% Called by the abstract machine, if no clauses exist for a predicate
'$undefp'([M|G]) :-
	'$find_goal_definition'(M, G, NM, NG),
	'$execute0'(NG, NM).

'$find_goal_definition'(M, G, NM, NG) :-
	% make sure we do not loop on undefined predicates
        % for undefined_predicates.
	'$enter_undefp',
	(
	 '$get_undefined_pred'(G, M, Goal, NM)
	->
	 '$exit_undefp'
	;
	 once('$find_undefp_handler'(G, M, Goal, NM))
	),
	!,
	Goal \= fail,
	'$complete_goal'(M, Goal, NM, G, NG).

'$complete_goal'(M, G, CurMod, G0, NG) :-
	  (
	   '$is_metapredicate'(G,CurMod)
	  ->
	   '$meta_expansion'(G, CurMod, M, M, NG,[])
	  ;
	   NG = G
	  ).

'$find_undefp_handler'(G,M,NG,user) :-
	functor(G, Na, Ar),
	user:exception(undefined_predicate,M:Na/Ar,Action), !,
	'$exit_undefp',
	(
	 Action == fail
	->
	 NG = fail
	;
	 Action == retry
	->
	 NG = G
	;
	 Action = error
	->
	 '$unknown_error'(M:G)
	;
	 '$do_error'(type_error(atom, Action),M:G)
	).
'$find_undefp_handler'(G,M,NG,user) :-
	\+ '$undefined'(unknown_predicate_handler(_,_,_), user),
	'$system_catch'(unknown_predicate_handler(G,M,NG), user, Error, '$leave_undefp'(Error)), !,
	'$exit_undefp'.
'$find_undefp_handler'(G,M,US,user) :-
	recorded('$unknown','$unknown'(M:G,US),_), !,
	'$exit_undefp'.
'$find_undefp_handler'(_,_,_,_) :-
	'$exit_undefp',
	fail.

'$leave_undefp'(Ball) :-
	'$exit_undefp',
	throw(Ball).


/* This is the break predicate,
	it saves the importante data about current streams and
	debugger state */

break :-
	nb_getval('$system_mode',SystemMode),
	nb_getval('$trace',Trace),
	nb_setval('$trace',off),
	nb_getval('$debug_jump',Jump),
	nb_getval('$debug_run',Run),
	'$debug_on'(Debug),
	'$debug_on'(false),
	nb_getval('$break',BL), NBL is BL+1,
	nb_getval('$spy_gn',SPY_GN),
	b_getval('$spy_glist',GList),
	b_setval('$spy_glist',[]),
	nb_setval('$break',NBL),
	current_output(OutStream), current_input(InpStream),
	format(user_error, '% Break (level ~w)~n', [NBL]),
	'$do_live',
	!,
	set_value('$live','$true'),
	b_setval('$spy_glist',GList),
	nb_setval('$spy_gn',SPY_GN),
	set_input(InpStream), 
	set_output(OutStream),
	'$debug_on'(Debug),
	nb_setval('$debug_jump',Jump),
	nb_setval('$debug_run',Run),
	nb_setval('$trace',Trace),
	nb_setval('$break',BL),
	nb_setval('$system_mode',SystemMode).

'$silent_bootstrap'(F) :-
	'$init_globals',
	nb_setval('$if_level',0),
	get_value('$lf_verbose',OldSilent),
	set_value('$lf_verbose',silent),
	set_stream(user_input,alias('$loop_stream')),
	bootstrap(F),
	% -p option must be processed after initializing the system
	'$init_path_extensions',
	set_value('$lf_verbose', OldSilent).

bootstrap(F) :-
%	'$open'(F, '$csult', Stream, 0, 0, F),
%	'$file_name'(Stream,File),
	open(F, read, Stream),
	stream_property(Stream, file_name(File)),
	'$start_consult'(consult, File, LC),
	file_directory_name(File, Dir),
	getcwd(OldD),
	cd(Dir),
	(
	  get_value('$lf_verbose',silent)
	->
	  true
	;
	  H0 is heapused, '$cputime'(T0,_),
	  format(user_error, '~*|% consulting ~w...~n', [LC,F])
	),
	'$loop'(Stream,consult),
	cd(OldD),
	'$end_consult',
	(
	  get_value('$lf_verbose',silent)
	->
	  true
	;
	  H is heapused-H0, '$cputime'(TF,_), T is TF-T0,
	  format(user_error, '~*|% ~w consulted ~w bytes in ~d msecs~n', [LC,F,H,T])
	),
	!,
	close(Stream).


'$init_path_extensions' :-
	get_value('$extend_file_search_path',P), !,
	P \= [],
	set_value('$extend_file_search_path',[]),
	'$extend_file_search_path'(P).
'$init_path_extensions'.

'$loop'(Stream,Status) :-
	repeat,
%VSC		( '$current_stream'(_,_,Stream) -> true
%VSC		 ; '$abort_loop'(Stream)
%VSC		),
		prompt1('|     '), prompt(_,'| '),
		'$current_module'(OldModule),
		'$system_catch'('$enter_command'(Stream,Status), OldModule, Error,
			 user:'$LoopError'(Error, Status)),
	!.

'$enter_command'(Stream,Status) :-
	'$read_vars'(Stream,Command,_,Pos,Vars),
	'$command'(Command,Vars,Pos,Status).

'$abort_loop'(Stream) :-
	'$do_error'(permission_error(input,closed_stream,Stream), loop).

/* General purpose predicates				*/

'$head_and_body'((H:-B),H,B) :- !.
'$head_and_body'(H,H,true).

%
% split head and body, generate an error if body is unbound.
%
'$check_head_and_body'((H:-B),H,B,P) :- !,
	'$check_head'(H,P).
'$check_head_and_body'(H,H,true,P) :-
	'$check_head'(H,P).

'$check_head'(H,P) :- var(H), !,
	'$do_error'(instantiation_error,P).
'$check_head'(H,P) :- number(H), !,
	'$do_error'(type_error(callable,H),P).
'$check_head'(H,P) :- db_reference(H), !,
	'$do_error'(type_error(callable,H),P).
'$check_head'(_,_).

% term expansion
%
% return two arguments: Expanded0 is the term after "USER" expansion.
%                       Expanded is the final expanded term.
%
'$precompile_term'(Term, Expanded0, Expanded, BodyMod, SourceMod) :-
	'$module_expansion'(Term, Expanded0, ExpandedI, BodyMod, SourceMod), !,
	(
	 '$access_yap_flags'(9,1)      /* strict_iso on */
        ->
	 Expanded = ExpandedI,
	 '$check_iso_strict_clause'(Expanded0)
        ;
	 '$expand_array_accesses_in_term'(ExpandedI,Expanded)
	).
'$precompile_term'(Term, Term, Term, _, _).
	

expand_term(Term,Expanded) :-
	( '$current_module'(Mod), \+ '$undefined'(term_expansion(_,_), Mod),
	  '$notrace'(Mod:term_expansion(Term,Expanded))
        ; \+ '$undefined'(term_expansion(_,_), system),
	  '$notrace'(system:term_expansion(Term,Expanded))
        ; \+ '$undefined'(term_expansion(_,_), user),
	  '$notrace'(user:term_expansion(Term,Expanded))
        ;
	  '$expand_term_grammar'(Term,Expanded)
	),
	!.


%
% Grammar Rules expansion
%
'$expand_term_grammar'((A-->B), C) :-
	'$translate_rule'((A-->B),C), !.
'$expand_term_grammar'(A, A).

%
% Arithmetic expansion
%
'$expand_term_arith'(G1, G2) :-
	get_value('$c_arith',true),
	'$c_arith'(G1, G2), !.
'$expand_term_arith'(G,G).


%
% Arithmetic expansion
%
'$expand_array_accesses_in_term'(Expanded0,ExpandedF) :-
	'$array_refs_compiled',
	'$c_arrays'(Expanded0,ExpandedF), !.
'$expand_array_accesses_in_term'(Expanded,Expanded).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   catch/throw implementation

% at each catch point I need to know:
% what is ball;
% where was the previous catch	
catch(G, C, A) :-
	'$catch'(C,A,_),
	yap_hacks:current_choice_point(CP0),
	'$execute'(G),
	yap_hacks:current_choice_point(CP1),
	(CP0 == CP1 -> !; true ).

% makes sure we have an environment.
'$true'.


% system_catch is like catch, but it avoids the overhead of a full
% meta-call by calling '$execute0' instead of $execute.
% This way it
% also avoids module preprocessing and goal_expansion
%
'$system_catch'(G, M, C, A) :-
	% check current trail
	'$catch'(C,A,_),
	yap_hacks:current_choice_point(CP0),
	'$execute_nonstop'(G, M),
	yap_hacks:current_choice_point(CP1),
	(CP0 == CP1 -> !; true ).

%
% throw has to be *exactly* after system catch!
%
throw(_Ball) :-
	% use existing ball
	'$get_exception'(Ball),
	!,
	'$jump_env_and_store_ball'(Ball).
throw(Ball) :-
	% get current jump point
	'$jump_env_and_store_ball'(Ball).


% just create a choice-point
'$catch'(_,_,_).
'$catch'(_,_,_) :- fail.

'$handle_throw'(_, _, _).
'$handle_throw'(C, A, _Ball) :-
	'$reset_exception'(Ball),
        % reset info
	('catch_ball'(Ball, C) ->
	    '$execute'(A)
	    ;
	    throw(Ball)
	).

'catch_ball'(Abort, _) :- Abort == '$abort', !, fail.
% system defined throws should be ignored by used, unless the
% user is hacking away.
'catch_ball'(Ball, V) :-
	var(V),
	nonvar(Ball),
	Ball = error(Type,_), % internal error ??
	functor(Type, Name, _),
	atom_codes(Name, [0'$|_]), %'0
	!, fail.
'catch_ball'(C, C).

'$run_toplevel_hooks' :-
	nb_getval('$break',0),
	recorded('$toplevel_hooks',H,_), !,
	( '$oncenotrace'(H) -> true ; true).
'$run_toplevel_hooks'.

'$enter_system_mode' :-
	nb_setval('$system_mode',on).

'$exit_system_mode' :-
	nb_setval('$system_mode',off),
	( nb_getval('$trace',on) -> '$creep' ; true).

%
% just prevent creeping from going on...
%
'$notrace'(G) :-
	'$disable_creep', !,
	(
		% creep was going on...
	 yap_hacks:current_choice_point(CP0),
	 '$execute'(G),
	 yap_hacks:current_choice_point(CP1),
	 ( CP0 == CP1 ->
	   !,
	   '$creep'
	 ;
	   (
	    '$creep'
	   ;
	    '$disable_docreep',
	    fail
	   )
	 )
	;
	 '$creep',
	 fail
	).
'$notrace'(G) :-
	'$execute'(G).

'$oncenotrace'(G) :-
	'$disable_creep', !,
	(
	 '$execute'(G)
	->
	 '$creep'
	;
	 '$creep',
	 fail
	).	
'$oncenotrace'(G) :-
	'$execute'(G), !.


'$run_at_thread_start' :-
	recorded('$thread_initialization',M:D,_),
	'$notrace'(M:D),
	fail.
'$run_at_thread_start'.


nb_getval(GlobalVariable, Val) :-
	'$nb_getval'(GlobalVariable, Val, Error),
	(var(Error)
	->
	 true
	;
	 '$getval_exception'(GlobalVariable, Val, nb_getval(GlobalVariable, Val)) ->
	 nb_getval(GlobalVariable, Val)
	;
	 '$do_error'(existence_error(variable, GlobalVariable),nb_getval(GlobalVariable, Val))
	).
		    

b_getval(GlobalVariable, Val) :-
	'$nb_getval'(GlobalVariable, Val, Error),
	(var(Error)
	->
	 true
	;
	 '$getval_exception'(GlobalVariable, Val, b_getval(GlobalVariable, Val)) ->
	 true
	;
	 '$do_error'(existence_error(variable, GlobalVariable),b_getval(GlobalVariable, Val))
	).
		    
access_file(File, Mode) :-
	swi_access_file(File, Mode).
expand_file_name(Exp, Matches) :-
	swi_expand_file_name(Exp, Matches).
file_base_name(File, Base) :-
	swi_file_base_name(File, Base).
exists_directory(Directory) :-
	swi_exists_directory(Directory).
file_directory_name(File, Dir) :-
	swi_file_directory_name(File, Dir).
file_name_extension(File, Name, Extension) :-
	swi_file_name_extension(File, Name, Extension).
is_absolute_file_name(File) :-
	swi_is_absolute_file_name(File).
prolog_to_os_filename(Prolog, OS) :-
	swi_prolog_to_os_filename(Prolog, OS).
same_file(File1, File2) :-
	swi_same_file(File1, File2).
time_file(File, Time) :-
	swi_time_file(File, Time).
working_directory(Old, New) :-
	swi_working_directory(Old, New).

cd(Dir) :- working_directory(_, Dir).

getcwd(Dir) :- working_directory(Dir, Dir).

close(Stream) :-
	swi_close(Stream).
close(Stream, Options) :-
	swi_close(Stream, Options).
open(File, Type, Stream) :-
	swi_open(File, Type, Stream).
open(File, Type, Stream, Opts) :-
	swi_open(File, Type, Stream, Opts).
open_null_stream(S) :-
	swi_open_null_stream(S).

atom_to_term(Atom, Term, Bindings) :-
	swi_atom_to_term(Atom, Term, Bindings).
term_to_atom(Term, Atom) :-
	swi_term_to_atom(Term, Atom).

with_output_to(Output, G) :-
	swi_with_output_to(Output, G).

at_end_of_stream :-
	swi_at_end_of_stream.
at_end_of_stream(Stream) :-
	swi_at_end_of_stream(Stream).
byte_count(Stream, Count) :-
%	format('~w~n',byte_count(Stream, Count)),
	swi_byte_count(Stream, Count).
character_count(Stream, Count) :-
%	format('~w~n',character_count(Stream, Count)),
	swi_character_count(Stream, Count).
current_stream(File,Mode,Stream) :-
	swi_current_stream(File,Mode,Stream).
is_stream(Stream) :-
%	format('~w~n',is_stream(Stream)),
	swi_is_stream(Stream).
line_count(Stream, Lines) :-
%	format('~w~n',line_count(Stream)),
	swi_line_count(Stream, Lines).
line_position(Stream, Position) :-
%	format('~w~n',line_position(Stream)),
	swi_line_position(Stream, Position).
set_stream(Stream, Property) :-
%	format('~w~n',set_stream(Stream,Property)),
	swi_set_stream(Stream, Property).
stream_property(Stream, Property) :-
%	format('~w~n',stream_property(Stream,Property)),
	swi_stream_property(Stream, Property).
set_stream_position(Stream, Position) :-
%	format('~w~n',stream_property(Stream,Property)),
	swi_set_stream_position(Stream, Position).

prompt1(X) :-
	swi_prompt1(X).
prompt(Old, New) :-
	swi_prompt(Old, New).

flush_output :-
	swi_flush_output.
flush_output(Stream) :-
	swi_flush_output(Stream).
ttyflush :-
	swi_ttyflush.

current_input(Stream) :-
	swi_current_input(Stream).
current_output(Stream) :-
	swi_current_output(Stream).
set_input(Stream) :-
	swi_set_input(Stream).
set_output(Stream) :-
	swi_set_output(Stream).

get(C) :-
	swi_get(C).
get(Stream, C) :-
	swi_get(Stream, C).
get0(C) :-
	swi_get0(C).
get0(Stream, C) :-
	swi_get0(Stream, C).
get_byte(C) :-
	swi_get_byte(C).
get_byte(Stream, C) :-
	swi_get_byte(Stream, C).
get_char(C) :-
	swi_get_char(C).
get_char(Stream, C) :-
	swi_get_char(Stream, C).
get_code(C) :-
	swi_get_code(C).
get_code(Stream, C) :-
	swi_get_code(Stream, C).

peek_byte(C) :-
	swi_peek_byte(C).
peek_byte(Stream, C) :-
	swi_peek_byte(Stream, C).
peek_char(C) :-
	swi_peek_char(C).
peek_char(Stream, C) :-
	swi_peek_char(Stream, C).
peek_code(C) :-
	swi_peek_code(C).
peek_code(Stream, C) :-
	swi_peek_code(Stream, C).

put(C) :-
	swi_put(C).
put(Stream, C) :-
	swi_put(Stream, C).
put_byte(C) :-
	swi_put_byte(C).
put_byte(Stream, C) :-
	swi_put_byte(Stream, C).
put_char(C) :-
	swi_put_char(C).
put_char(Stream, C) :-
	swi_put_char(Stream, C).
put_code(C) :-
	swi_put_code(C).
put_code(Stream, C) :-
	swi_put_code(Stream, C).

skip(C) :-
	swi_skip(C).
skip(Stream, C) :-
	swi_skip(Stream, C).

nl :-
	swi_nl.
nl(Stream) :-
	swi_nl(Stream).
print(T) :-
	swi_print(T).
print(Stream, T) :-
	swi_print(Stream, T).
write(T) :-
	swi_write(T).
write(Stream, T) :-
	swi_write(Stream, T).
writeq(T) :-
	swi_writeq(T).
writeq(Stream, T) :-
	swi_writeq(Stream, T).
write_canonical(T) :-
	swi_write_canonical(T).
write_canonical(Stream, T) :-
	swi_write_canonical(Stream, T).
write_term(Stream, T) :-
	swi_write_term(Stream, T).
write_term(Stream, T, Options) :-
	swi_write_term(Stream, T, Options).
tab(C) :-
	swi_tab(C).
tab(Stream, C) :-
	swi_tab(Stream, C).

append(File) :-
	swi_append(File).
see(File) :-
	swi_see(File).
seeing(File) :-
	swi_seeing(File).
seen :-
	swi_seen.
tell(File) :-
	swi_tell(File).
telling(File) :-
	swi_telling(File).
told :-
	swi_told.

format(Command, Args) :-
	swi_format(Command, Args).
format(Stream, Command, Args) :-
	swi_format(Stream, Command, Args).

is_stream(Stream) :-
	swi_is_stream(Stream).

'$raw_read'(Stream, String) :-
	'swi_$raw_read'(Stream, String).

