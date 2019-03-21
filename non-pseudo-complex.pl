/*---------------------------------------------------------------*/
% Test 1 first builds the complex of Example 7.4.2 and checks that
% one and just one triangle has three incident tetrahedra.
% Next a copy of the totally exploded complex is generated and
% all simplex gluing instructions are applied.
% Finally we verify that the two processes lead to the same complex
/*---------------------------------------------------------------*/
:-consult("iqm.pl").
% Builds the complex of example of Example 7.4.2 
buildComplex:-
	buildFrame,
	% add the three central simplices
	addSimplex(34,[x,y,z,a]),
  	addSimplex(35,[x,y,z,b]),
  	addSimplex(36,[x,y,z,c]).
% Builds the complex of the example without 3 central simplices.
buildFrame:-buildRedFrame,buildGreenFrame,buildBlueFrame.
buildRedFrame:-
  addSimplex(1,[x,r,z,a]),
  addSimplex(2,[x,d,z,c]),
  addSimplex(3,[x,s,z,b]),
  addSimplex(16,[x,r,a,i]),
  addSimplex(17,[x,i,a,d]),
  addSimplex(18,[x,i,c,d]),
  addSimplex(19,[x,i,c,s]),
  addSimplex(20,[x,i,s,b]).
buildGreenFrame:-
  addSimplex(23,[x,y,t,a]),
  addSimplex(24,[x,y,p,c]),
  addSimplex(25,[x,y,u,b]),
  addSimplex(7,[y,j,a,t]),
  addSimplex(8,[y,j,a,p]),
  addSimplex(9,[y,j,p,c]),
  addSimplex(21,[y,j,c,u]),
  addSimplex(22,[y,j,u,b]).
buildBlueFrame:-
  addSimplex(31,[y,z,v,a]),
  addSimplex(32,[y,z,e,c]),
  addSimplex(33,[y,z,w,b]),
  addSimplex(26,[k,z,v,a]),
  addSimplex(27,[k,z,a,e]),
  addSimplex(28,[k,z,c,e]),
  addSimplex(29,[k,z,c,w]),
  addSimplex(30,[k,z,w,b]).
test1:-
  resetComplex,
  buildComplex,
  listNonPseudoManifolds,  
  resetDecomposition,
  buildTotallyExploded,
  stitchTheComplex,
  not(splitVertex(_)),
  vtIsoTv.
% if main succeds then  there is no splitting vertex i.e.
% decomposition has reached the original complex.
% Nevertheless check the two are isomorphic.
test2:-
  resetComplex,
  buildComplex,
  dm1connectedComponents(C), % the complex is 2-connected.
  print_term(C,[]),
  boundary(Bnd), % it has a boundary surface
  resetComplex,
  addComplex(bnd,Bnd),
  dm1connectedComponents(CBnd), % made up of a single surface
  print_term(CBnd,[]),
  listNonPseudoManifolds. % nine triangles on the boundary are non manifolds.
% they are triplets of triangles around edges xy xz zy.
% round xz
% [d,x,z]
% [r,x,z]
% [s,x,z]
% round xy
% [p,x,y]
% [t,x,y]
% [u,x,y]
% round yz
% [e,y,z]
% [v,y,z]
% [w,y,z]
% hence all non manifold edges are of order three.
% [x,y] of order 3
% [x,z] of order 3
% [y,z] of order 3

test3:-
  resetComplex,
  buildFrame,
  dm1connectedComponents(C), % the complex has three  2-connected
  % components. Each is made up of eight tetrahedra.
  print_term(C,[]),
  boundary(Bnd), % it has a  boundary surface that has a 1-connected
  % boundary component of 54 triangles.
  resetComplex,
  addComplex(bnd,Bnd),
  dumpTv,
  dm1connectedComponents(CBnd), % made up of a single surface
  print_term(CBnd,[]),
  listNonPseudoManifolds.
% 36 triangles on the boundary are non manifolds.
% these are organized in groups of four around xg yg zg for g=a,c,b
% This seems to sum up to 36 but each group of four shares two triangles
% with others. For instance the group round xa is xat xay xad xaz 
% but xay is 
% shared with group round ya and xaz is shared with group round za.
% Thus instead of 36 we have 27.
% round ax
% [a,x,y] shared
% [a,d,x]
% [a,t,x]
% [a,x,z] shared
% round ay
% [a,p,y]
% [a,v,y]
% [a,y,z] shared
% round az
% [a,e,z]
% [a,r,z]
% round bx
% [b,i,x]
% [b,u,x]
% [b,x,y]
% [b,x,z]
% round by
% [b,j,y]
% [b,w,y]
% [b,y,z]
% round bz
% [b,k,z]
% [b,s,z]
% round cx
% [c,p,x]
% [c,s,x]
% [c,x,y]
% [c,x,z]
% round cy
% [c,e,y]
% [c,u,y]
% [c,y,z]
% round cz
% [c,d,z]
% [c,w,z]
% The reamaining nine are in group of three around xy, xz and zy
% round xy
% [p,x,y]
% [t,x,y]
% [u,x,y]
% round xz
% [d,x,z]
% [r,x,z]
% [s,x,z]
% round yz
% [e,y,z]
% [v,y,z]
% [w,y,z]
% [x,a] of order 4
% [x,b] of order 4
% [x,c] of order 4
% [y,a] of order 4
% [y,b] of order 4
% [y,c] of order 4
% [z,a] of order 4
% [z,b] of order 4
% [z,c] of order 4
% The nine edges center of these groups i.e.  xg yg zg for g=a,c,b are no
% the only non manifold edges also xz, yz, zy are adjacent to six 
% triangles each but they are in the above mentionde group above. 
% [x,y] of order 6
% [x,z] of order 6
% [y,z] of order 6
test4:-
  resetComplex,
  buildComplex,
  resetDecomposition,
  buildTotallyExploded,
  stitchTheFrame, % just build the three cones.
  dumpSplitVertex, 
  listNonPseudoManifolds,
  vtToTv(Frame), % compute the  TV relation for this frame.
  resetComplex, % consider this and forget the other
  maplist(assert,Frame), % this install the new TV in Frame 
  listNonPseudoManifolds, % no non-pseudomanifolds
  dm1connectedComponents(C), % complex has 6 2-connected components
  % Three are isolated tetrahedra.
  print_term(C,[]),
  boundary(Bnd), % it has a boundary surface that has 3 1-connected
  % components each of 16 triangles plus three tetrahedra.
  resetComplex,
  addComplex(bnd,Bnd), 
  % now the TV holds the boundary
  dm1connectedComponents(CBnd), % boundary is made up of 6 surfaces
  print_term(CBnd,[]),
  listNonPseudoManifolds. % no pseudomanifolds in the boundary
% one can see that the three cones in figure form three distinct complexes. Indeed we have three vertex 
% copies for a,b,c and within each complex x,y,z are split in three.
% Be careful z is split in three in the complex in the red frame 
% and not split in the blue frame
% x is split in three in the complex in the green frame 
% and not split in the red frame
% y is split in three in the complex in the blue frame 
% and not split in the green frame
% Thus we have four splits for x,y,z and three for a,b,c for a total of 21.
stitchTheComplex:-
  stitchTheRedFrame,
  doPseudoManifoldGluingInstruction(1,34),
  doPseudoManifoldGluingInstruction(23,34),
  stitchTheGreenFrame,
  doPseudoManifoldGluingInstruction(31,34),
  stitchTheBlueFrame,
  doPseudoManifoldGluingInstruction(3,35),
  doPseudoManifoldGluingInstruction(25,35),
  doPseudoManifoldGluingInstruction(33,35),
  doPseudoManifoldGluingInstruction(2,36),
  doPseudoManifoldGluingInstruction(24,36),
  doPseudoManifoldGluingInstruction(32,36).
stitchTheFrame:-
     stitchTheRedFrame,stitchTheGreenFrame,stitchTheBlueFrame.
stitchTheRedFrame:-
  doPseudoManifoldGluingInstruction(1,16),
  doPseudoManifoldGluingInstruction(16,17),
  doPseudoManifoldGluingInstruction(17,18),
  doPseudoManifoldGluingInstruction(18,2),
  doPseudoManifoldGluingInstruction(18,19),
  doPseudoManifoldGluingInstruction(19,20),
  doPseudoManifoldGluingInstruction(20,3).
% completed the cone of x
stitchTheGreenFrame:-
  doPseudoManifoldGluingInstruction(23,7),
  doPseudoManifoldGluingInstruction(7,8),
  doPseudoManifoldGluingInstruction(8,9),
  doPseudoManifoldGluingInstruction(9,24),
  doPseudoManifoldGluingInstruction(9,21),
  doPseudoManifoldGluingInstruction(21,22),
  doPseudoManifoldGluingInstruction(25,22).
% completed the cone of y
stitchTheBlueFrame:-
  doPseudoManifoldGluingInstruction(31,26),
  doPseudoManifoldGluingInstruction(26,27),
  doPseudoManifoldGluingInstruction(27,28),
  doPseudoManifoldGluingInstruction(28,32),
  doPseudoManifoldGluingInstruction(28,29),
  doPseudoManifoldGluingInstruction(29,30),
  doPseudoManifoldGluingInstruction(33,30).
% completed the cone of z
% Shellable construction of this complex adds a couple of spanning % tetrahedra in the end. A shellable complex is homotopy is
% equivalent to a wedge sum of spheres, one for each spanning
% simplex and of corresponding dimension.
% They are 35 and 36. 
stitchOneSpanning:-
  stitchTheRedFrame,
  doPseudoManifoldGluingInstruction(1,34),
  doPseudoManifoldGluingInstruction(23,34),
  stitchTheGreenFrame,
  doPseudoManifoldGluingInstruction(31,34),
  stitchTheBlueFrame,
  doPseudoManifoldGluingInstruction(3,35),
  doPseudoManifoldGluingInstruction(25,35),
  doPseudoManifoldGluingInstruction(33,35).
% and the corresponding test is
test5:-
  resetComplex,
  buildComplex,
  retract(tv(36,_)), % forget last tetrahedron to add,
  resetDecomposition,
  buildTotallyExploded,
  stitchTheRedFrame,
  doPseudoManifoldGluingInstruction(1,34),
  doPseudoManifoldGluingInstruction(23,34),
  stitchTheGreenFrame,
  doPseudoManifoldGluingInstruction(31,34),
  stitchTheBlueFrame,
  doPseudoManifoldGluingInstruction(3,35),
  doPseudoManifoldGluingInstruction(25,35),
  doPseudoManifoldGluingInstruction(33,35),
  %stitchOneSpanning,
  %build 3 cones and add one spanning tetrahedron 
  dumpSplitVertex, 
  vtToTv(Sphere), % compute the  TV relation for this frame,
  resetComplex, % consider this and forget the other
  maplist(assert,Sphere),% this install the new TV in Frame 
  listNonPseudoManifolds,% must be a sphere no non-pseudomanifolds
  dm1connectedComponents(C),% complex has a 2-connected components
  print_term(C,[]),
  boundary(Bnd),% its boundary has a 1-connected surfaces
  print_term(Bnd,[]),
  resetComplex,
  addComplex(bnd,Bnd), 
  dm1connectedComponents(CBnd),% boundary is made up of a surfaces
  print_term(CBnd,[]),
  listNonPseudoManifolds, % no pseudomanifolds in the boundary
  boundary(BBnd), % the boundary of the boundary is empty
  print_term(BBnd,[]), % there are no boundary simplices.
  eulerX(XChar),nl,write(XChar),nl.
