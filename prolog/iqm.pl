:- module(iqm,[resetComplex/0,
 vtToTv/1,vtIsoTv/0,addSimplex/2,addComplex/2,
 incident/2,skeleton/2,star/2,link/2,eulerX/1,
 orderOf/2,adjacent/2,adjacent/3,dm1connectedComponents/1,
 nonPseudoManifold/1,nonPseudoManifold/2,nonPseudoManifoldPair/2,
 printSimplex/1,listNonPseudoManifolds/0,listAdjacents/0,
 boundaryface/1,boundary/1,dumpTv/0,
 resetDecomposition/0,buildTotallyExploded/0,
 doGluingInstruction/2,doPseudoManifoldGluingInstruction/2,
 doVertexEquation/3,dumpDecomp/0,
 splitVertex/1,dumpSplitVertex/0,tv/2,vt/2]).
/** <module> iqm library for Initial Quasi Manifolds

 This program takes a simplicial complex defined via a TV relation and
 handles its totally decomposed version giving the possibility to
 stitch together top simplexes either using vertex equations or simplex
 gluing instructions. The totally decomposed complex is represented
 via VT relation

@author Franco Morando 
@license GPL
*/
:-consult("utilities.pl").
/*------------------------------------------------------------------*/
% Handling of the complex defined via the TV relation.
% The complex is stored as a TV relation using dynamic predicate tv/2
/*------------------------------------------------------------------*/
%! tv(?T:atom,?V:list) is nondet
%
% is a dynamic predicate used to encode the TV for complex to be decomposed. 
% If vt(foo,[a,b,c]) is stored then triangle foo with vertices a,b and c exist.
% User are advised to use carefully this predicate possibly ignoring it.
%
% @arg T an atom that is the index of some top d-simplex.
% @arg V a set of atoms  that are indexes of vertices.
:-dynamic tv/2.
/*------------------------------------------------------------------*/
%! vt(?V:atom,?T:list) is nondet
%
% is a dynamic predicate used to encode the VT for a decomposition of
% the  complex stored in the TV. 
% If tv(a,[foo,bar]) is stored then top simplexes foo and bar share
% vertex a.  To be consistent with the rest of the package no two
% entries like tv(a,[foo,bar]) and tv(a,[some,thing]) must exist.
% User are advised to use carefully this predicate possibly ignoring it.
%
% @arg V an atom that is the index of some vertex.
% @arg T a set of atoms  that are indexes of  top simplexes.
/*------------------------------------------------------------------*/
:-dynamic vt/2.
/*------------------------------------------------------------------*/
%! resetComplex is det 
%
% deletes the TV relation, complex is erased.
resetComplex:-
 retractall(tv(_,_)),!.
/*------------------------------------------------------------------*/
%! vtToTv(-Res:list) is det
% This procedure extracts in Res a TV from the VT representation of
% the complex obtained by stitching  the totally decomposed complex.
% Indeed the stitched complex is recorded in this package
% in a VT using asserts to vt(V,L). 
% This procedure leaves dynamic predicates tv/2 and vt/2 unchanged.
%
% @arg Res a list of terms of the form vt(<some vertex>,[s1,...,s1])
% representing a VT relation.
vtToTv_(_):-
 retractall(rtv(_,_)),reset_gensym,
 vt(V,Adj),gensym(V,Vs),select(T,Adj,_),
	 findrtv(T,L),assert(rtv(T,[Vs|L])),fail.
vtToTv_(Res):-!,
 findall(U,(rtv(T,VL),U=tv(T,VL)),Res).
vtToTv(Res):-!,
 vtToTv_(Res),
 retractall(rtv(_,_)).
findrtv(T,L):-retract(rtv(T,L)),!.
findrtv(_,[]):-!.
/*------------------------------------------------------------------*/
%! vtIsoTv is det
% succeeds if tv and vt are isomorphic.
% We assume that VT is created by a call to buildTotallyExploded and 
% several calls to doPseudoManifoldGluingInstruction or to doGluingInstruction.
% We recall that the effect of these two calls is the same but
% doPseudoManifoldGluingInstruction checks that stitching simplex is a pseudomanifold
% one. If not it fails.
:-dynamic isoname/2.
:-dynamic rtv/2.
isoverts(VL,RVL):-
	isorenamelist(VL,RVL1),!,
	sort(RVL1,SRVL1),!,
	sort(RVL,SRVL1),!.
isorenamelist(VL,RVL):-maplist(isorename,VL,RVL).
isorename(N,RN):-isoname(N,RN),!.
isorename(N,RN):-isoname(N,RNN),RN\=RNN,!,fail.
isorename(N,RN):-gensym(N,RN),assert(isoname(N,RN)),!.
vtIsoTv:-vtToTv_(TV1),reset_gensym,retractall(isoname(_,_)),
 findall(U,(tv(T,VL),U=tv(T,VL)),TV2),
 length(TV1,L),length(TV2,L),
 forall(tv(T,VL),(rtv(T,RVL),isoverts(VL,RVL))),
 retractall(rtv(_,_)),retractall(isoname(_,_)),!.
vtIsoTv:-retractall(rtv(_,_)),retractall(isoname(_,_)),fail.
/*------------------------------------------------------------------*/
%! addSimplex(+Gamma:atom,+SetGamma:list) is det
% adds a simplex to this complex
% Top simplexes are encoded by atoms (could be integers)
%
% @arg Gamma an atom that is used as an index for the top simplex to be added.
% @arg SetGamma a set of atoms used to encode vertices for this top simplex. 
addSimplex(Gamma,SetGamma):-
 assert(tv(Gamma,SetGamma)).
/*------------------------------------------------------------------*/
%! addComplex(+U:atom,+SimplexList:list) is det
% adds a complex made up of a list of simplexes, e.g. [[a,b,c],[b,c,d]]
% for the rectangle abdc.
% Simplex indexes are atoms U1, U2 etc. that are created randomly
% using U as prefix.
%
% @arg U an atom that is used as a prefix for top simplex indexes.
% @arg SimplexList a set of sets of vertices. 
% Vertices must be encoded by atoms. e.g.  [[a,b,c],[b,c,d]].
addComplex(U,SimplexList):-
	forall(member(Simplex,SimplexList),
		(gensym(U,Id),!,addSimplex(Id,Simplex))
		).
/*------------------------------------------------------------------*/
% Topological checks
% Once that the complex is stored in the TV relation some topological
% checks can be performed.
/*------------------------------------------------------------------*/
%! eulerX(+X:int) is det
% returns the Euler characteristics of the closed 2-complex in TV.
% Works only if in TV is stored a closed 2-complex. Otherwise 
% results are meaningless. 
%
% @arg X the Euler characteristics of surface. 
eulerX(XEuler):-
 skeleton(V,0),skeleton(E,1),skeleton(F,2),
 length(V,NV),length(E,NE),length(F,NF),
 XEuler is NV-NE+NF.
 % (orientable *-> G is (2-XEuler)/2 ;   G is 2-XEuler ).
/*------------------------------------------------------------------*/
%! link(+V:atom,-S:list) is det
% returns the set of all maximal simplices in TV in the link of V. 
%
% @arg V a vertex
% @arg S a set of lists each being a maximal simplex.
link(V,S):-
 findall(F,(tv(_,L),select(V,L,F)),S).
/*------------------------------------------------------------------*/
%! star(+V:atom,-S:list) is det 
% returns the set of all top simplices in TV in the star of V. 
%
% @arg V a vertex
% @arg S a set of atoms each being a top simplex index.
star(V,S):-
 findall(F,(tv(F,L),member(V,L)),S).
/*------------------------------------------------------------------*/
%! skeleton(-S:list,+D:int) is det
% returns the set of all D simplices in TV. 
%
% @arg S a set of sets each being a simplex.
% @arg D dimension of the simplexes to be returned D=0 for points
skeleton(Set,D):-
 L is D+1,findall(F,simplex_(F,L),S),list_to_set(S,Set).
simplex_(F,L):-tv(_,Sup),asublist(F,Sup,L).
/*------------------------------------------------------------------*/
%! incident(?Theta:atom,+S:list) is nondet
% succeeds iff  Theta is incident to simplex S.
%
% @arg Theta an atom that is the index of some top simplex.
% @arg S a set of vertices.
incident(Theta,S):- tv(Theta,Sup),subset(S,Sup).
/*------------------------------------------------------------------*/
%! orderOf(+S:list,-N:int) is det
% counts in N the number of top simplexes that are incident to
% the simplex given by the set of vertices in S.
%
% @arg S a set of vertices.
% @arg N a positive  integer that gives the number of top simplexes that are incident to S.
% 
orderOf(S,N):-findall(Theta,incident(Theta,S),Incidents),!,
		length(Incidents,N).

/*------------------------------------------------------------------*/
%! adjacent(?Theta1:atom,?Theta2:atom,?SetTheta:list) is nondet
% succeeds iff Theta1 and Theta2 are adjacent via the set of vertexes
% in SetTheta. If Theta1 is equal to Theta2 it fails.
% To succeed Theta1 and Theta2 must be two top d-simplexes and they must share
% the d-1 face in SetTheta.
%
% @arg Theta1 an atom that is the index of some top simplex.
% @arg Theta2 an atom not equal to Theta1 that is the index of some top simplex.
% @arg SetTheta a set of vertices.
adjacent(Theta1,Theta2,SetTheta):-
  tv(Theta1,SetTheta1),tv(Theta2,SetTheta2),
  Theta1\==Theta2,length(SetTheta1,D),length(SetTheta2,D),
  intersection(SetTheta1,SetTheta2,SetTheta),D1 is D-1,length(SetTheta,D1).
/*------------------------------------------------------------------*/
%! adjacent(?Theta:atom,?Theta2:atom) is nondet
% succeeds iff Theta1 and Theta2 are adjacent.
% If Theta1 is equal to Theta2 it fails.
% Theta1 and Theta2 must be two d-simplexes and they must share
% a d-1 face.
%
% @arg Theta1 an atom that is the index of some top simplex.
% @arg Theta2 an atom not equal to Theta1 that is the index of some top simplex.
adjacent(Theta1,Theta2):-adjacent(Theta1,Theta2,_).
/*------------------------------------------------------------------*/
%! dm1connectedComponents(-C:list) is det
% always succeeds and returns a set of sets. Each set contains atoms that
% are indexes of top simplices. Top simplicies in this set are uniformly
% dimensional. All top d-simplices in such sets are d-1 connected.
% 
% @arg C a set of sets of top simplexes. 
dm1connectedComponents(C):-
	findall(T,tv(T,_),TL),!,
	closurePartition(adjacent,TL,C).
/*------------------------------------------------------------------*/
%! nonPseudoManifold(?Theta1:atom) is nondet
% succeed if Theta1 is a top d-simplex d-1-adjacent to three or more 
% top d simplexes.
%
% @arg Theta1 an atom that is the index of some top simplex.
nonPseudoManifold(Theta1):-
	nonPseudoManifold(Theta1,_).
%! nonPseudoManifold(?Theta1:atom,?SetTheta:list) is nondet
% succeeds iff the top d-simplex Theta1 is d-1 adjacent via the set 
% of vertexes in SetTheta to three or more top d-simplexes. 
%
% @arg Theta1 an atom that is the index of some top d-simplex.
% @arg SetTheta a set of  d-1 vertices.
nonPseudoManifold(Theta1,SetTheta):-
  tv(Theta1,_),tv(Theta2,_),tv(Theta3,_),
  Theta2\==Theta1,Theta1\==Theta3,Theta2\==Theta3,
  adjacent(Theta1,Theta2,SetTheta),adjacent(Theta2,Theta3,SetTheta).
/*------------------------------------------------------------------*/
%! nonPseudoManifoldPair(?Theta1:atom,?Theta2:atom) is nondet
% succeed if two top d-simplexes Theta1 and Theta2 meet at a 
% non-manifold d-1-face.
%
% @arg Theta1 an atom that is the index of some top simplex.
% @arg Theta2 an atom not equal to Theta1 that is the index of some top simplex.
nonPseudoManifoldPair(Theta1,Theta2):-
  tv(Theta1,SetTheta1),tv(Theta2,SetTheta2),
  Theta2\==Theta1,
  intersection(SetTheta1,SetTheta2,SetTheta),
  adjacent(Theta1,Theta3,SetTheta),
  Theta1\==Theta3,Theta2\==Theta3.
/*------------------------------------------------------------------*/
%! printSimplex(+S:atom) is det.
% is a printing utility that prints the TV relation for the top simplex S.
%
% @arg S an atom that is the index of some top simplex.
printSimplex(S):-tv(S,L),print_term(tv(S,L),[]),nl,!.
/*------------------------------------------------------------------*/
%! listNonPseudoManifolds is det
% is a listing utility that  lists the set of vertices for all top 
% d-simplexes involved in a non-manifold d-1-adjacency for some d.
% The utility lists also all d-1 simplexes involved in a non-manifold
% adjacency of top d-simplexes for some d.
listNonPseudoManifolds:-
  nl,write("Non PseudoManifold Top d-Simplexes:"),nl,
  setof(Theta1,nonPseudoManifold(Theta1),Results),Results\=[],
  findall(SSimpl,
   (select(NMS,Results,_),tv(NMS,Simpl),sort(Simpl,SSimpl)),
  SResults),
  sort(SResults,SortedSResults),
  findall(NMS1,
  (select(NMS1,SortedSResults,_),write(NMS1),nl),_),!,
   write("Non Manifold d-1-Simplexes:"),nl,
  findall(Linking,nonPseudoManifold(_,Linking),ResultsLinks),!,
  list_to_set(ResultsLinks,NonPseudoManifoldsLinks),
  sort(NonPseudoManifoldsLinks,SortedNonPseudoManifoldsLinks),
  findall(NMS,
  (member(NMS,SortedNonPseudoManifoldsLinks),
   write(NMS),write(" of order "),orderOf(NMS,N), write(N),nl),_),!.
listNonPseudoManifolds:- write("no non-manifold"),nl,!.
/*------------------------------------------------------------------*/
%! listAdjacents is det
% is a listing utility that  lists the set of triples 
% [Theta1,Theta2,SetTheta] for all top  d-simplexes Theta1,Theta2
% involved in a d-1-adjacency  via the set of vertices in SetTheta 
% for some d.
listAdjacents:-
  setof([Theta1,Theta2,SetTheta],adjacent(Theta1,Theta2,SetTheta),Results),
  write(Results),nl,!.
listAdjacents:- write("no adjacents"),nl,!.
/*------------------------------------------------------------------*/
%! maxfaceof(?T:atom,?F:list) is nondet
% returns a top d-simplex T and, upon backtracking all its d-1 faces F
maxfaceof(T,F):-
	tv(T,V),select(_,V,F).
/*------------------------------------------------------------------*/
% succeeds if T is  a top d-simplex  and F is a d-1 face of some other
% top d-simplex.
%! covered(+T:atom,+F:list) is det
%
% @arg T an atom that is the index of some top d-simplex.
% @arg F a set of  d-1 vertices.
covered(T,F):-
	 tv(T2,V2),T\==T2,select(_,V2,F1),permutation(F,F1),!.
/*------------------------------------------------------------------*/
%! boundaryface(+F:list) is det
% succeeds if F is d-1 face of a top d-simplex  and F is on the
% boundary.
%
% @arg F a set of  d-1 vertices.
boundaryface(F):-
	tv(T,V),select(_,V,F),
	not(covered(T,F)),!.
/*------------------------------------------------------------------*/
%! boundary(-Bnd:list) is det
% returns the list of boundary simplexes.
%
% @arg Bnd a set of sets of  vertices.
boundary(Bnd):-
	findall(F,(maxfaceof(_,F),boundaryface(F)),Bnd),!.
/*------------------------------------------------------------------*/
%! dumpTv is det
% lists dynamic predicate tv(T,[v1,...]) that encodes  
% the TV for the original complex.
dumpTv:-
 nl,tv(T,L),write(['simplex ',T,'incident to vertexes ',L]),nl,fail.
dumpTv:-!.

/*------------------------------------------------------------------*/
% primitives to handle the decompositions of the complex in the TV
% by bookkeeping of
% what we called standardized vertices, e.g. $v_\Sigma$,
% encoded as vt(v,Sigma).
% This SECOND  complex is represented by means of its VT relation.
% vt(v,Sigma) iff v is in Sigma.
/*------------------------------------------------------------------*/
%! resetDecomposition is det
%
% resets the stored  decomposition  by deleting the VT.
resetDecomposition:-
 retractall(vt(_,_)),!.
/*------------------------------------------------------------------*/
%! buildTotallyExploded is det
%
% The complex given by the TV relation is turned into a totally
% exploded version where all top simplexes are  distinct components
% and vertex v is split into n-copies being n the number of top
% simplices incident to v in the TV complex.
% The result is the creation of a VT relation storing vt(foo,[bar])
% one for every vertex foo and for every top simplex bar.
% Vertices  of the form  vt(foo,[bar]) are called "vertex copies" of 
% the "splitting vertex" foo.
% As stitcing operation takes place relations like vt(foo,[bar]) and
% vt(foo,[biz]) might disappear and will be substituted by vt(foo,[bar,biz]).
buildTotallyExploded:-
 forall(tv(Gamma,SetGamma),checklist(addVertex(Gamma),SetGamma)),!.
/*------------------------------------------------------------------*/
%! addVertex(+Gamma:atom,+V:atom) is det
% adds a vertex V for the simplex Gamma in the totally exploded version.
addVertex(Gamma,V):-assert(vt(V,[Gamma])).
/*------------------------------------------------------------------*/
%! doVertexEquation(+Theta1:atom,+Theta2:atom,+V:atom) is det
%
% Vertexes =|$v_{\sigma_1}$|= and =|$v_{\sigma_2}$|= with =|$\theta_1\in\sigma_1$|= and
%  =|$\theta_2\in\sigma_2$|= are identified. Note that the resulting
%  complex, modeled by the VT relation, is a decomposition of the 
%  complex modeled by the TV relation.
%
% @arg Theta1 an atom that is the index of some top simplex.
% @arg Theta2 an atom, distinct from Theta1, that is the index of some top simplex.
% @arg V a vertex that must be a vertex of Theta1 and Theta2 in the TV relation.
doVertexEquation(Theta1,Theta2,V):-
 vt(V,Sigma1),member(Theta1,Sigma1),
 vt(V,Sigma2),member(Theta2,Sigma2),!,
 retract(vt(V,Sigma1)),retract(vt(V,Sigma2)),
 append(Sigma1,Sigma2,Sigma),!,
 assert(vt(V,Sigma)),!.
/*------------------------------------------------------------------*/
%! doGluingInstruction(+Theta1:atom,+Theta2:atom) is det
%
% merges vertexes according to simplex gluing instruction Theta1 <-> Theta2.
%
% @arg Theta1 an atom that is the index of some top simplex.
% @arg Theta2 an atom, distinct from Theta1, that is the index of some top simplex.
doGluingInstruction(Theta1,Theta2):-
  tv(Theta1,SetTheta1),tv(Theta2,SetTheta2),
  intersection(SetTheta1,SetTheta2,SetTheta),
  checklist(doVertexEquation(Theta1,Theta2),SetTheta),!.
/*------------------------------------------------------------------*/
%! doPseudoManifoldGluingInstruction(+Theta1:atom,+Theta2:atom) is det
% merges vertexes according to simplex gluing instruction Theta1 <-> Theta2.
% Before executing checks that Theta1 and Theta2 are pseudomanifold adjacent in the TV.
%
% @arg Theta1 an atom that is the index of some top simplex.
% @arg Theta2 an atom, distinct from Theta1, that is the index of some top simplex.
doPseudoManifoldGluingInstruction(Theta1,Theta2):-
  not(nonPseudoManifoldPair(Theta1,Theta2)),!,
  doGluingInstruction(Theta1,Theta2).
/*------------------------------------------------------------------*/
%! dumpDecomp is det
% dumps each simplex for  the current decomposition. The dump lists
% each vertex with the corresponding VT record.
dumpDecomp:-nl,tv(T,_),write('simplex '),write(T),write('=[ '),dump_simplex(T),write(']'),nl,fail.
dumpDecomp:-!.
/*------------------------------------------------------------------*/
%! dump_simplex is det
% is a utility that dumps the list of top simplices adjacent to V.
dump_simplex(T):-vt(V,L),member(T,L),sort(L,PL),
 write(V),write('-'),write(PL),write(' '),fail.
dump_simplex(_):-!.
/*------------------------------------------------------------------*/
%! splitVertex(?V:atom) is nondet
% is a topological check on the decomposition  that
% succeeds if the vertex V has more than one vertex copy.
splitVertex(V):-vt(V,L),vt(V,LL),LL\==L.

/*------------------------------------------------------------------*/
%! dumpSplitVertex is det
% dumps the VT for  all splitting vertexes.
dumpSplitVertex:-nl,vt(V,L),once(splitVertex(V)),print_term(vt(V,L),[]),nl,fail.
dumpSplitVertex:-!.
