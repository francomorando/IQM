%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program  tests boundary predicate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
:-consult("iqm.pl").
main:-
  resetComplex,
  addSimplex(1,[a,b,c]),
  addSimplex(2,[b,c,d]),
  addSimplex(3,[b,a,e]),
  addSimplex(4,[a,c,f]),
  boundary(Bnd),
  print_term(Bnd,[]),
  resetComplex,
  addSimplex(5,[a,b,c,d]),
  addSimplex(6,[a,b,c,e]),
  boundary(Bnd3),
  print_term(Bnd3,[]).
