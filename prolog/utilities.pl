:- module(utilities,[closurePartition/3, disjoint/2, write_ln/1,
 asublist/3]).
/** <module> utilities for application

@author Franco Morando 
@license GPL
*/

/*---------------------------------------------------------------*/
%! write_ln(+X) is det.
% Write X on a single line
%
% @arg X the item to be written.
write_ln(X):-write(X),!,nl.
/*---------------------------------------------------------------*/
%! asublist(?Sub:list,+L:list,+N:int) is nondet.
% upon backtrack returns all ordered sublist of L of length N.
%
% @arg Sub sublist returned.
% @arg L the complete list.
% @arg N length of the list to be returned. 
asublist(_,L,N):-length(L,N1),N>N1,!,fail.
asublist(L,L,N):-length(L,N),!.
asublist([],_,0):-!.
asublist([X|Sub],[X|L],N):-N1 is N-1,asublist(Sub,L,N1).
asublist(Sub,[_|L],N):-asublist(Sub,L,N).
/*---------------------------------------------------------------*/
%! expand(:Rel,+El,-Class:list) is det.
% One hop expand of element El.
expand(Rel,El,[El|Class]):-
	setof(X,call(Rel,El,X),Class),!.
expand(_,El,[El]).
/*---------------------------------------------------------------*/
%! expandl(:Rel,+L:list,-Class:list) is det.
% One hop expand of equivalence class  L.
expandl(Rel,L,Class):-
	maplist(expand(Rel),L,ToFlat),!,
	flatten(ToFlat,Flat),!,
	list_to_set(Flat,Class),!.
/*---------------------------------------------------------------*/
%! expandp(:Rel,+Part:list,-BigPart:list) is det.
% One hop expand of partition Part.
expandp(Rel,Part,BigPart):-
	maplist(expandl(Rel),Part,ToMerge),!,
	mergeClasses(ToMerge,BigPart),!.
/*---------------------------------------------------------------*/
%! mergeClasses(+ToMerge:list,-BigPart:list) is det.
mergeClasses([],[]).
mergeClasses([C|Others],Result):-
	mergeClasses(Others,Presult),!,
	mergeClass(C,Presult,Result),!.
mergeClass(C,[],[C]).
mergeClass(C,Others,[Merged|Disjoint]):-
	partition(disjoint(C),Others,Disjoint,ToMerge),!,
	flatten([C|ToMerge],ListMerged),!,
	list_to_set(ListMerged,Merged),!.
/*---------------------------------------------------------------*/
%! disjoint(+C1:list,+C2:list) is det.
% succeeds if C1 and C2 are disjoint.
%
% @arg C1 and C2 two lists to be considered as sets
disjoint(C1,C2):-intersection(C1,C2,[]),!.
/*---------------------------------------------------------------*/
%! closurePartition(:Rel,+L:list,-Partition:list) is det.
% returns the partition of L given by the quotient L/Rel*.
% The binary relation Rel must be
% symmetric and Rel* is the transitive closure of Rel.
%
% @arg Rel a binary relation.
% @arg L a set of elements as a list.
% @arg Partition a list of sets one for each set in the partition.
% of L defined by the transitive closure of Rel.
closurePartition(Rel,L,Partition):-
        maplist(expand(Rel),L,CoarsePartition),!,
	expandp(Rel,CoarsePartition,FinePart),!,
	maplist(length,CoarsePartition,Sign1),!,
	maplist(length,FinePart,Sign2),!,
	sort(Sign1,L1),!,
	sort(Sign2,L2),!,
	((L1==L2->Partition=FinePart);expandp(Rel,FinePart,Partition)).
