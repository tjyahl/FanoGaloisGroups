---------------------
---------------------
-- FanoProblems.m2 --
---------------------
---------------------
--
-- Thomas Yahl
-- January 25, 2022
--

newPackage(
    "FanoProblems",
    Version=>"0.1.0",
    Authors=>{{
	Name=>"Thomas Yahl",
	Email=>"Thomasjyahl@tamu.edu",
	Homepage=>"https://math.tamu.edu/~thomasjyahl"
	}},
    Headline=>"Methods to numerically compute monodromy groups of finite Fano problems.",
    PackageImports=>{"MonodromySolver","NumericalAlgebraicGeometry"},
    PackageExports=>{"NumericalAlgebraicGeometry"},
    DebuggingMode=>true
    )

export{
    --------------------
    --Monodromy object--
    --------------------
    "monodromy",
    "family",
    "basePoint",
    "baseSolutions",
    "computedPermutations",
    
    ----------------
    --Fano Methods--
    ----------------
    "FanoMonodromy",
    "FanoEquations",
    "FanoNumSolns",
    
    ---------------------
    --Monodromy methods--
    ---------------------
    "monodromyLoop",
    "changeBasePoint",
    "saveMonodromy",
    "loadMonodromy",
    "printCycleTally",
    "printForGAP",
    
    -----------
    --options--
    -----------
    "BaseElements",
    "BaseRing",
    "NumPaths",
    "RefineSolns",
    "RefineStartSolns",
    "RefineEndSolns",
    "LoopScale"
    }


------------------------------------------
------------ Monodromy Object ------------
------------------------------------------

Monodromy = new Type of HashTable
Monodromy.synonym = "monodromy group"
Monodromy.GlobalAssignHook = globalAssignFunction
Monodromy.GlobalReleaseHook = globalReleaseFunction

net (Monodromy) := String => M->(
    "Monodromy of family "|toString(M#family)
    )

--Main constructor of Monodromy objects.
monodromy = method(Options=>{BaseElements=>{}})
monodromy (PolySystem,Point,List) := Monodromy => o->(F,basePt,baseSolns)->(
    monod := new Monodromy from {
	family => F,
	basePoint => basePt,
	baseSolutions => baseSolns,
	computedPermutations => new MutableHashTable from o.BaseElements,
	cache => new CacheTable
	};
    monod    
    )

--?
monodromy (List,List,List) := Monodromy => o->(F,basePt,baseSolns)->(
    monodromy(polySystem F, point {basePt}, apply(baseSolns,p-> point {p}))
    )


--Checks validity of Monodromy object.
----change this.
isWellDefined (Monodromy) := Boolean => M->(
    F := M#family;
    
    R := ring F;
    numVars := numgens R;
    
    C := coefficientRing R;
    numParams := numgens C;
    
    --Checks for types, square system, number of coordinates
    if not all(M#baseSolutions,p->class p === Point) then (
	print("--isWellDefined: Expected solutions of class Point");
	return false
	);
    if not (class M#basePoint === Point) then (
	print("--isWellDefined: Expected base point to be of class Point");
	return false
	);
    if not (F#NumberOfPolys === F#NumberOfVariables) then (
	print("--isWellDefined: Expected square system");
	return false
	);
    if not ((#(coordinates M#basePoint) === numParams) and all(M#baseSolutions,p->#(coordinates p) === numVars)) then (
	print("--isWellDefined: Dimension mismatch in base point or solutions");
	return false
	);
    
    --Check all elements of group are permutations of solutions
    if not all(keys M#computedPermutations,P->(set P === set toList(0..#(M#baseSolutions)-1)) and (#P === #(M#baseSolutions))) then (
	print("--isWellDefined: Not all group elements are permutations");
	return false
	);
    
    true
    )


--------------------------------------
------------ Fano Methods ------------
--------------------------------------

--Constructs a Monodromy object corresponding to the problem of r-planes in P^n 
----on the intersection of polynomials of degrees d = (d_1,..,d_s).
FanoMonodromy = method(Options=>{Verbose=>false})
FanoMonodromy (ZZ,ZZ,Sequence) := Monodromy => o->(r,n,degs)->(
    if not ((n-r)*(r+1) === sum(degs,d->binomial(d+r,r))) then error "--FanoMonodromy: Expected family of non-empty finite Fano problems.";
    
    eqns := FanoEquations(r,n,degs);
    numSolns := FanoNumSolns(r,n,degs);
    
    --Solve via MonodromySolver
    (basePt,baseSolns) := solveFamily(eqns,TargetSolutionCount=>numSolns,Verbose=>o.Verbose,NumberOfNodes=>3);
    baseSolns = points baseSolns;
    if (#baseSolns =!= numSolns) then (
	print("--FanoMonodromy: Not all solutions found, please try again");
	return(0)
	);
    
    monodromy(eqns,basePt,baseSolns)
    )


--Generates the equations of a Fano problem
FanoEquations = method(Options=>{BaseRing=>CC})
FanoEquations (ZZ,ZZ,Sequence) := PolySystem => o->(r,n,degs)->(
    if not ((n-r)*(r+1) === sum(degs,d->binomial(d+r,r))) then error "--FanoMonodromy: Expected family of non-empty finite Fano problems.";
    
    s := #degs;
    
    --Coefficient ring for hyperplane equations
    w := symbol w;
    params := splice apply(s,i->w_(i,1)..w_(i,binomial(degs#i+n,n)));
    C := (o.BaseRing)[params];
    
    --Ring for hyperplane equations
    y := symbol y;
    S := C[y_1..y_n];
    
    --Coefficient ring for hyperplane equations restricted to linear space (defined by the x_(i,j))
    x := symbol x;
    GrVars := toList(x_(1,1)..x_(n-r,r+1));
    R := C[GrVars];
    
    --Ring for hyperplane equations restricted to linear space
    t := symbol t;
    T := R[t_1..t_r];
    
    --Hyperplane equations
    F := apply(s,i->basis(0,degs#i,S)*(transpose matrix {toList(w_(i,1)..w_(i,binomial(degs#i+n,n)))}));
    
    --Parameterization of linear space determined by the x_(i,j)
    L := apply(n-r,j->sum(toList(x_(j+1,1)..x_(j+1,r)),toList(t_1..t_r),(z,v)->z*v)+x_(j+1,r+1))|toList(t_1..t_r);
    restrict := map(T,S,L);
    
    FonL := apply(F,f->restrict f);
    coeffs := polySystem flatten apply(FonL,f->flatten entries sub(last coefficients f,R));
    
    rewriteEqns(coeffs,BaseRing=>o.BaseRing)
    )


--Computes the number of solutions to a generic finite Fano problem of 
----the given data.
FanoNumSolns = method()
FanoNumSolns (ZZ,ZZ,Sequence) := ZZ => (r,n,degs)->(
    x := symbol x;
    R := ZZ[x_0..x_r];
    
    Q := product(degs,d->(
    	partitionList := apply(subsets(d+r,r),S->(
       		S = {-1}|S|{d+r};
       		apply(r+1,i->S#(i+1)-S#i-1)
       		));
       	product(partitionList,s->sum(s,gens R,(a,x)->a*x))
       	));
    V := product(r+1,i-> product((i+1)..r,j->x_i-x_j) );
    mon := product(r+1,i->x_i^(n-i));
    
    coefficient(mon,Q*V)
    )


-----------------------------------------
------------ Utility Methods ------------
-----------------------------------------

--Specializes a family of systems to a base point(s) in the parameter space.
specializeSys = method()
specializeSys (PolySystem,Point) := PolySystem => (F,pt)->(
    R := ring F;
    n := numgens R;
    z := symbol z;
    S := CC[z_1..z_n];
    phi := map(S,R,toList(z_1..z_n)|(coordinates pt));
    
    polySystem apply(equations F,f->phi(f))
    )

specializeSys (PolySystem,List) := List => (F,pts)->(
    R := ring F;
    n := numgens R;
    z := symbol z;
    S := CC[z_1..z_n];

    for pt in pts list (
	phi := map(S,R,toList(z_1..z_n)|(coordinates pt));
    	polySystem apply(equations F,f->phi(f))
	)
    )


--rewrites equations in k[parms][variables] without doublely-indexed variables.
rewriteEqns = method(Options=>{BaseRing=>CC})
rewriteEqns (PolySystem) := PolySystem => o->F->(
    R := ring F;
    n := numgens R;
    m := numgens coefficientRing R;
    
    x := symbol x;
    a := symbol a;
    
    newR := (o.BaseRing)[a_1..a_m][x_1..x_n];
    phi := map(newR,R,toList(x_1..x_n)|toList(a_1..a_m));

    polySystem apply(equations F,f->phi(f))
    )


--Determines if a set of Points corresponds to a permutation of the solutions
----over the base of a given Monodromy object.
----can change max to 2-norm or other norm. (max scales well with dimension)
permutation = method(Options=>{Tolerance=>1e-6})
permutation (Monodromy,List) := List => o->(M,S)->(
    numSolns := #(M#baseSolutions);
    P := apply(S,p->position(toList(0..numSolns-1),i->norm(2,coordinates M#baseSolutions#i - coordinates p) < o.Tolerance));
    
    if (set P === set toList(0..numSolns-1)) then (
	P
	) else (
	error "--permutation: No permutation found between given list of Points and base solutions."
	)
    )


--Gives a permutation of {0,..,n} given in list form in cycle notation
printPermutation = method()
printPermutation (List) := String => P->(
    cycles := {};
    L := toList(0..#P-1);
    while (#L>0) do (
	i := first L;
	if (P#i =!= i) then (
	    c := {i+1};
	    while (P#i =!= first L) do (
	    	L = delete(P#i,L);
	    	c = append(c,P#i+1);
	    	i = P#i;
	    	);
	    cycles = append(cycles,c)
	    );
	L = drop(L,1)
	);
    
    if (#cycles > 0) then (
	concatenate apply(cycles,l->toString toSequence l)
	) else (
	"()"
	)
    )


--Refines the solutions of a given monodromy object. Uses
----the same options as 'refine' given in "NumericalAlgebraicGeometry".
refine (Monodromy) := Monodromy => o->M->(
    F' := specializeSys(M#family,M#basePoint);
    newSolns := refine(F',M#baseSolutions,o);
    grp := pairs M#computedPermutations;
    
    monodromy(M#family,M#basePoint,newSolns,BaseElements=>grp)
    )


-------------------------------------------
------------ Monodromy Methods ------------
-------------------------------------------

--Tracks solutions from base point to given points and back to base point.
----The resulting permutation of solutions is recorded.
monodromyLoop = method(Options=>{RefineSolns=>false,Tolerance=>1e-4,Verbosity=>0,LoopScale=>.5,stepIncreaseFactor=>3,tStep=>1e-8,maxCorrSteps=>2,CorrectorTolerance=>1e-6,tStepMin=>1e-30})
monodromyLoop (Monodromy,List) := Monodromy => o-> (M,pts)->(
    F := M#family;
    R := ring F;
    basePt := M#basePoint;
    baseSolns := M#baseSolutions;
	    
    pts = {basePt}|pts|{basePt};

    trackedSolns := baseSolns;
    
    --Specialized systems
    specSystems := specializeSys(F,pts);
    
    --Tracking options consolidated for NAG
    trackingOpts := new OptionTable from {
	stepIncreaseFactor=>o.stepIncreaseFactor,
	tStep=>o.tStep,
	maxCorrSteps=>o.maxCorrSteps,
	CorrectorTolerance=>o.CorrectorTolerance,
	tStepMin=>o.tStepMin
	};
    
    --Tracking base solutions along loop defined by given points
    for i in 0..(#pts-2) do (
	trackedSolns = track(specSystems#i,specSystems#(i+1),trackedSolns,trackingOpts);
	if (o.Verbosity>2) then print("Finished tracking part "|toString(i+1)|" of "|toString(#pts-1));
	if any(trackedSolns,p->status p != Regular) then error "--monodromyLoop: Tracking failure"
	);
    
    --Refine solutions 
    if (o.RefineSolns) then (
	trackedSolns = refine(specSystems#0,trackedSolns,Bits=>30)
	);
    
    --Checking that tracked solutions correspond to permutation of base solutions
    try(perm := permutation(M,trackedSolns,Tolerance=>o.Tolerance)) else (
	if (o.Verbosity>0) then print("--monodromyLoop: No permutation from tracked solutions");
	if (o.Verbosity>1) then error "--monodromyLoop: Break in function";
	return M
	);
    if (o.Verbosity>0) then print printPermutation perm;
    
    --Adjoin permutation to the list of permutations
    if M#computedPermutations#?perm then (
	M#computedPermutations#perm = M#computedPermutations#perm + 1
	) else (
	M#computedPermutations#perm = 1
	);
    
    M
    )

--Tracks solutions from base point along a random triangular path. The 
----resulting permutation of solutions is recorded.
monodromyLoop (Monodromy) := Monodromy => o-> M->(
    n := numgens coefficientRing ring (M#family);
    pts := apply(2,i->point {coordinates M#basePoint + o.LoopScale*apply(n,j->(random(CC)-random(CC)) )});
    
    monodromyLoop(M,pts,o)
    )

--Tracks solutions from base point along a given number of random triangular
----paths. The resulting permutations of solutions are recorded.
monodromyLoop (Monodromy,ZZ) := Monodromy => o->(M,n)->(
    for i in 1..n do (
	try(monodromyLoop(M,o)) else if (o.Verbosity>0) then print("Path tracking failed");
	if (o.Verbosity>0) then print("Loop "|toString(i)|" complete")
	);
    
    M
    )


--Tracks solutions from the base point to a new chosen point. A new
----Monodromy object is returned with the chosen base point.
----NumPaths tracks along multiple paths and accumulates solutions
changeBasePoint = method(Options=>{RefineStartSolns=>false,RefineEndSolns=>false,Tolerance=>1e-4,NumPaths=>10,stepIncreaseFactor=>1.2,tStep=>1e-4,maxCorrSteps=>40,CorrectorTolerance=>1e-8,Verbose=>false})
changeBasePoint (Monodromy,Point) := Monodromy => o->(M,bp)->(
    F := M#family;
    basePt := M#basePoint;
    baseSolns := M#baseSolutions;
    grp := pairs M#computedPermutations;
    
    --Specialized systems
    specSystems := specializeSys(F,{basePt,bp});
    
    --Tracking options consolidated for NAG
    trackingOpts := new OptionTable from {
	gamma=>exp(2*pi*ii*random(RR)),
	stepIncreaseFactor=>o.stepIncreaseFactor,
	tStep=>o.tStep,
	maxCorrSteps=>o.maxCorrSteps,
	CorrectorTolerance=>o.CorrectorTolerance
	};
    
    --Refine solutions 
    if (o.RefineStartSolns) then (
	baseSolns = refine(specSystems#0,baseSolns)
	);
    
    --Tracking solutions
    if o.Verbose then print("Tracking 1st time");
    trackedSolns := track(specSystems#0,specSystems#1,baseSolns,trackingOpts);
    trackedSolns = select(trackedSolns,p->status p == Regular);

    --Refine solutions 
    if (o.RefineEndSolns) then (
	trackedSolns = refine(specSystems#1,trackedSolns,Bits=>30)
	);

    newBaseSolns := {};
    for soln in trackedSolns do (
	if all(newBaseSolns,oldSoln->norm(2,coordinates oldSoln - coordinates soln) > o.Tolerance) then (
	    if o.Verbose then print("solution #"|toString(#newBaseSolns)|" found");
	    if o.Verbose then print(toString(coordinates soln));
	    if o.Verbose then print("");
	    newBaseSolns = newBaseSolns|{soln}
	    )
	);
    if o.Verbose then print("Counted "|toString(#newBaseSolns)|" solutions out of "|toString(#M#baseSolutions));
    COUNT := 1;
    while (#newBaseSolns < #baseSolns) and (COUNT < o.NumPaths-1) do (
	if o.Verbose then print("Tracking "|toString(COUNT+1)|"-th time");
	trackedSolns = track(specSystems#0,specSystems#1,baseSolns,trackingOpts++{gamma=>random(CC)});
	trackedSolns = select(trackedSolns,p->status p == Regular);
	
	--Refine solutions 
    	if (o.RefineEndSolns) then (
	    trackedSolns = refine(specSystems#1,trackedSolns)
	    );
	
	for soln in trackedSolns do (
	    if all(newBaseSolns,oldSoln->norm(2,coordinates oldSoln - coordinates soln) > o.Tolerance) then (
		newBaseSolns = newBaseSolns|{soln};
		if o.Verbose then print("solution #"|toString(#newBaseSolns)|" found");
	    	if o.Verbose then print(toString(coordinates soln));
	    	if o.Verbose then print("")
		)
	    );
	if o.Verbose then print("Counted "|toString(#newBaseSolns)|" solutions out of "|toString(#M#baseSolutions));
	COUNT = COUNT + 1
	);
    
    if (#newBaseSolns === #baseSolns) then (
	if o.Verbose then print("All solutions to new base system were found");
	monodromy(F,bp,newBaseSolns,BaseElements=>grp)
	) else (
	error "--changeBasePoint: Couldn't find all solutions over new base."
	)
    )


--Writes a Monodromy object to a file.
saveMonodromy = method()
saveMonodromy (Monodromy,String) := Nothing => (M,s)->(
    F := M#family;
    R := ring F;
    C := coefficientRing R;
    s << toString gens R << endl << toString gens C << endl << toExternalString (M#family) << endl << toExternalString (M#basePoint) << endl << toExternalString (M#baseSolutions) << endl << toString (pairs M#computedPermutations) << endl << toString (apply(keys M#cache,h->h=>M#cache#h)) << close
    )


--Recovers a Monodromy object from a file.
loadMonodromy = method()
loadMonodromy (String) := Monodromy => s->(
    L := lines get s;
    polyVars := value L#0;
    parms := value L#1;
    R := CC[parms][polyVars];
    F := value L#2;
    basePt := value L#3;
    baseSolns := value L#4;
    grp := value L#5;
    monodromyCache := value L#6;
    
    monodromy(F,basePt,baseSolns,BaseElements=>grp)
    )


--Prints tracked permutations in cycle notation.
printCycleTally = method()
printCycleTally (Monodromy) := Nothing => M->(
    T := M#computedPermutations;
    scan(keys T,sigma-> print(printPermutation(sigma) | ": " | toString(T#sigma)))
    )


--Outputs a string that represents the group generated by tracked
----permutations. Output is GAP readable.
printForGAP = method()
printForGAP (Monodromy) := String => M->(
    T := keys M#computedPermutations;
    if (#T > 0) then (
	S := apply(T,sigma->printPermutation(sigma));
    	"Group("|fold(S,(s,t)->s|","|t)|")"
	) else (
	"Group()"
	)
    )


---------------------------------------
------------ Documentation ------------
---------------------------------------









end

restart
loadPackage("FanoMonodromy")

