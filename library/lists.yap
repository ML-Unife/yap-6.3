% This file has been included as an YAP library by Vitor Santos Costa, 1999

%
% This file includes code from Bob Welham, Lawrence Byrd, and R. A. O'Keefe.
%
:- module(lists,[append/3,
	delete/3,
	is_list/1,
	last/2,
	member/2,
	memberchk/2,
	nextto/3,
	nth/3,
	nth/4,
	nth0/3,
	nth0/4,
	permutation/2,
	prefix/2,
	remove_duplicates/2,
	reverse/2,
	same_length/2,
	select/3,
	sublist/2,
	substitute/4,
	suffix/2,
	sumlist/2,
        list_concat/2
	]).


%   append(Prefix, Suffix, Combined)
%   is true when all three arguments are lists, and the members of Combined
%   are the members of Prefix followed by the members of Suffix.  It may be
%   used to form Combined from a given Prefix and Suffix, or to take a given
%   Combined apart.  E.g. we could define member/2 (from SetUtl.Pl) as
%	member(X, L) :- append(_, [X|_], L).

append([], L, L).
append([H|T], L, [H|R]) :-
	append(T, L, R).


%   delete(List, Elem, Residue)
%   is true when List is a list, in which Elem may or may not occur, and
%   Residue is a copy of List with all elements identical to Elem deleted.

delete([], _, []).
delete([Head|List], Elem, Residue) :-
	Head == Elem, !,
	delete(List, Elem, Residue).
delete([Head|List], Elem, [Head|Residue]) :-
	delete(List, Elem, Residue).


% is_list(List)
% is true when List is a proper List
%
is_list(L) :- var(L), !, fail.
is_list([]).
is_list([_|List]) :- is_list(List).

 
%   last(List, Last)
%   is true when List is a List and Last is identical to its last element.
%   This could be defined as last(L, X) :- append(_, [X], L).

last([H|List], Last) :-
	last(List, H, Last).

last([], Last, Last).
last([H|List], _, Last) :-
	last(List, H, Last).

%   member(?Element, ?Set)
%   is true when Set is a list, and Element occurs in it.  It may be used
%   to test for an element or to enumerate all the elements by backtracking.
%   Indeed, it may be used to generate the Set!

member(Element, [Element|_]).
member(Element, [_|Rest]) :-
        member(Element, Rest).


%   memberchk(+Element, +Set)
%   means the same thing, but may only be used to test whether a known
%   Element occurs in a known Set.  In return for this limited use, it
%   is more efficient when it is applicable.

memberchk(Element, [Element|_]) :- !.
memberchk(Element, [_|Rest]) :-
        memberchk(Element, Rest).

%   nextto(X, Y, List)
%   is true when X and Y appear side-by-side in List.  It could be written as
%	nextto(X, Y, List) :- append(_, [X,Y], List).
%   It may be used to enumerate successive pairs from the list.

nextto(X,Y, [X,Y|_]).
nextto(X,Y, [_|List]) :-
	nextto(X,Y, List).

%   nth0(+N, +List, ?Elem) is true when Elem is the Nth member of List,
%   counting the first as element 0.  (That is, throw away the first
%   N elements and unify Elem with the next.)  It can only be used to
%   select a particular element given the list and index.  For that
%   task it is more efficient than nmember.
%   nth(+N, +List, ?Elem) is the same as nth0, except that it counts from
%   1, that is nth(1, [H|_], H).

nth0(0, [Head|_], Head) :- !.

nth0(N, [_|Tail], Elem) :-
	nonvar(N),
	M is N-1,
	nth0(M, Tail, Elem).

nth0(N,[_|T],Item) :-		% Clause added KJ 4-5-87 to allow mode
	var(N),			% nth0(-,+,+)
	nth0(M,T,Item),
	N is M + 1.


nth(1, [Head|_], Head) :- !.

nth(N, [_|Tail], Elem) :-
	nonvar(N),
	M is N-1,			% should be succ(M, N)
	nth(M, Tail, Elem).

nth(N,[_|T],Item) :-		% Clause added KJ 4-5-87 to allow mode
	var(N),			% nth(-,+,+)
	nth(M,T,Item),
	N is M + 1.

%   nth0(+N, ?List, ?Elem, ?Rest) unifies Elem with the Nth element of List,
%   counting from 0, and Rest with the other elements.  It can be used
%   to select the Nth element of List (yielding Elem and Rest), or to 
%   insert Elem before the Nth (counting from 1) element of Rest, when
%   it yields List, e.g. nth0(2, List, c, [a,b,d,e]) unifies List with
%   [a,b,c,d,e].  nth is the same except that it counts from 1.  nth
%   can be used to insert Elem after the Nth element of Rest.

nth0(0, [Head|Tail], Head, Tail) :- !.

nth0(N, [Head|Tail], Elem, [Head|Rest]) :-
	nonvar(N),
	M is N-1,
	nth0(M, Tail, Elem, Rest).

nth0(N, [Head|Tail], Elem, [Head|Rest]) :-	% Clause added KJ 4-5-87
	var(N),					% to allow mode
	nth0(M, Tail, Elem, Rest),		% nth0(-,+,+,?).
	N is M+1.


nth(1, [Head|Tail], Head, Tail) :- !.

nth(N, [Head|Tail], Elem, [Head|Rest]) :-
	nonvar(N),
	M is N-1,
	nth(M, Tail, Elem, Rest).

nth(N, [Head|Tail], Elem, [Head|Rest]) :-	% Clause added KJ 4-5-87
	var(N),					% to allow mode
	nth(M, Tail, Elem, Rest),		% nth(-,+,+,?).
	N is M+1.


%   permutation(List, Perm)
%   is true when List and Perm are permutations of each other.  Of course,
%   if you just want to test that, the best way is to keysort/2 the two
%   lists and see if the results are the same.  Or you could use list_to_bag
%   (from BagUtl.Pl) to see if they convert to the same bag.  The point of
%   perm is to generate permutations.  The arguments may be either way round,
%   the only effect will be the order in which the permutations are tried.
%   Be careful: this is quite efficient, but the number of permutations of an
%   N-element list is N!, even for a 7-element list that is 5040.

permutation([], []).
permutation(List, [First|Perm]) :-
	select(First, List, Rest),	%  tries each List element in turn
	permutation(Rest, Perm).


% prefix(Part, Whole) iff Part is a leading substring of Whole

prefix([], _).
prefix([Elem | Rest_of_part], [Elem | Rest_of_whole]) :- 
  prefix(Rest_of_part, Rest_of_whole).

%   remove_dups(List, Pruned)
%   removes duplicated elements from List.  Beware: if the List has
%   non-ground elements, the result may surprise you.

remove_duplicates(List, Pruned) :-
	sort(List, Pruned).

%   reverse(List, Reversed)
%   is true when List and Reversed are lists with the same elements
%   but in opposite orders.  rev/2 is a synonym for reverse/2.

reverse(List, Reversed) :-
	reverse(List, [], Reversed).

reverse([], Reversed, Reversed).
reverse([Head|Tail], Sofar, Reversed) :-
	reverse(Tail, [Head|Sofar], Reversed).


%   same_length(?List1, ?List2)
%   is true when List1 and List2 are both lists and have the same number
%   of elements.  No relation between the values of their elements is
%   implied.
%   Modes same_length(-,+) and same_length(+,-) generate either list given
%   the other; mode same_length(-,-) generates two lists of the same length,
%   in which case the arguments will be bound to lists of length 0, 1, 2, ...

same_length([], []).
same_length([_|List1], [_|List2]) :-
	same_length(List1, List2).


%   select(?Element, ?Set, ?Residue)
%   is true when Set is a list, Element occurs in Set, and Residue is
%   everything in Set except Element (things stay in the same order).

select(Element, [Element|Rest], Rest).
select(Element, [Head|Tail], [Head|Rest]) :-
	select(Element, Tail, Rest).


%   sublist(Sublist, List)
%   is true when both append(_,Sublist,S) and append(S,_,List) hold.

sublist(Sublist, List) :-
	prefix(Sublist, List).
sublist(Sublist, [_|List]) :-
	sublist(Sublist, List).

%   substitute(X, XList, Y, YList)
%   is true when XList and YList only differ in that the elements X in XList
%   are replaced by elements Y in the YList.
substitute(X, XList, Y, YList) :-
	'$substitute'(XList, X, Y, YList).

'$substitute'([], _, _, []).
'$substitute'([X0|XList], X, Y, [Y|YList]) :-
	X == X0, !,
	'$substitute'(XList, X, Y, YList).
'$substitute'([X0|XList], X, Y, [X0|YList]) :-
	'$substitute'(XList, X, Y, YList).

%   suffix(Suffix, List)
%   holds when append(_,Suffix,List) holds. 
suffix(Suffix, Suffix).
suffix(Suffix, [_|List]) :-
	suffix(Suffix,List).

%   sumlist(Numbers, Total)
%   is true when Numbers is a list of integers, and Total is their sum.

sumlist(Numbers, Total) :-
	sumlist(Numbers, 0, Total).

sumlist([], Total, Total).
sumlist([Head|Tail], Sofar, Total) :-
	Next is Sofar+Head,
	sumlist(Tail, Next, Total).


%   list_concat(Lists, List)
%   is true when Lists is a list of lists, and List is the
%   concatenation of these lists.

list_concat(Lists, List) :-
	list_concat(Lists, [], List).

list_concat([], []).
list_concat([H|T], L) :-
	list_concat(H, L, Li),
	list_concat(T, Li).

list_concat([], L, L).
list_concat([H|T], [H|Lf], Li) :-
	list_concat(T, Lf, Li).




