----------------------
----------------------
-- FanoMonodromy.m2 --
----------------------
----------------------
--
-- Thomas Yahl
-- January 25, 2022
--

newPackage(
    "FanoMonodromy",
    Version=>"0.1.0",
    Authors=>{{
	Name=>"Thomas Yahl",
	Email=>"Thomasjyahl@tamu.edu",
	Homepage=>"https://math.tamu.edu/~thomasjyahl"
	}},
    Headline=>"Methods to numerically compute monodromy groups of finite Fano problems.",
    PackageImports=>{"MonodromySolver","NumericalAlgebraicGeometry","NumericalCertification"}, --maybe not numericalcert?
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
    "alphaConstants", --change to some form of soft verification
    "group", --change name
    
    ----------------
    --Constructors--
    ----------------
    "FanoMonodromy",
    
    ----------------
    --Core methods--
    ----------------
    "monodromyLoop",
    "changeBasePoint",
    "saveMonodromy",
    "loadMonodromy",
    "printCycleTally",
    "printForGAP",
    "testAlternatingMonodromy"
    
    -----------
    --options--
    -----------
    "BaseElements",
    "MonteCarlo", --change to accumulateSolns?
    "RefineSolns"
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
monodromy (PolySystem,Point,List,List) := Monodromy => o->(F,basePt,baseSolns,cons)->(
    monod := new Monodromy from {
	family => F,
	basePoint => basePt,
	baseSolutions => baseSolns,
	alphaConstants => cons,
	group => new MutableHashTable from o.BaseElements,
	cache => new CacheTable
	};
    
    monod    
    )

monodromy (PolySystem,Point,List) := Monodromy => o->(F,basePt,baseSolns)->(
    F' := specializeSys(F,basePt);
    cons := apply(baseSolns,p->computeConstants(F',p));
    
    monod := monodromy(F,basePt,baseSolns,cons,o);
    monod
    )

monodromy (List,List,List) := Monodromy => o->(F,basePt,baseSolns)->(
    M := monodromy(polySystem F, point {basePt}, apply(baseSolns,p-> point {p}));
    M
    )


--Checks validity of Monodromy object. Solutions are valid if they are certified.
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
    if not all(keys M#group,P->(set P === set toList(0..#(M#baseSolutions)-1)) and (#P === #(M#baseSolutions))) then (
	print("--isWellDefined: Not all group elements are permutations");
	return false
	);
    
    --Check all solutions are certified and distinct
    if any(M#alphaConstants,l-> l#0 > .03) then (
	if any(M#alphaConstants,l-> l#0 > (13-3*sqrt(17))/4) then (
	    print("--isWellDefined: Not all solutions certified");
	    return false
	    ) else (
	    print("--isWellDefined: Solutions not refined for tracking");
	    return false
	    )
	);
    baseSolns := M#baseSolutions;
    n := #baseSolns;
    S := subsets(n,2);
    if not all(S,s->norm(2,point{coordinates baseSolns#(s#0) - coordinates baseSolns#(s#1)})>2*(M#alphaConstants#(s#0)#1 + M#alphaConstants#(s#1)#1)) then (
	print("--isWellDefined: Not all solutions certified as distinct");
	return false
	);
    
    true
    )


--------------------------------------
------------ Constructor ------------
--------------------------------------

--Constructs a Monodromy object corresponding to the problem of r-planes in P^n 
----on the intersection of polynomials of degrees d = (d_1,..,d_s).
FanoMonodromy = method(Options=>{Verbose=>false})
FanoMonodromy (ZZ,ZZ,Sequence) := Monodromy => o->(r,n,degs)->(
    if not ((n-r)*(r+1) === sum(degs,d->binomial(d+r,r))) then error "--FanoMonodromy: Expected family of finite Fano problems.";
    
    s := #degs;
    numSolns := FanoNumSolns(r,n,degs);
    
    --Coefficient ring for hyperplane equations
    w := symbol w;
    params := splice apply(s,i->w_(i,1)..w_(i,binomial(degs#i+n,n)));
    C := CC[params];
    
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
    F := apply(s,i->basis(0,degs#i,module S)*(transpose matrix {toList(w_(i,1)..w_(i,binomial(degs#i+n,n)))}));
    
    --Parameterization of linear space determined by the x_(i,j)
    L := apply(n-r,j->sum(toList(x_(j+1,1)..x_(j+1,r)),toList(t_1..t_r),(z,v)->z*v)+x_(j+1,r+1))|toList(t_1..t_r);
    restrict := map(T,S,L);
    
    FonL := apply(F,f->restrict f);
    coeffs := polySystem flatten apply(FonL,f->flatten entries sub(last coefficients f,R));
    
    --Solve via MonodromySolver
    (basePt,baseSolns) := solveFamily(coeffs,TargetSolutionCount=>numSolns,Verbose=>o.Verbose,NumberOfNodes=>3);
    baseSolns = points baseSolns;
    if (#baseSolns =!= numSolns) then (
	print("--FanoMonodromy: Not all solutions found");
	return(0)
	);
    
    M := monodromy(rewriteEqns(coeffs),basePt,baseSolns);
    M
    )


--------------------------------------------------------
------------ Utility Methods (Not exported) ------------
--------------------------------------------------------

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
    
    numSols := coefficient(mon,Q*V);
    numSols
    )


--Specializes a family of systems to a base point(s) in the parameter space.
specializeSys = method()
specializeSys (PolySystem,Point) := PolySystem => (F,pt)->(
    R := ring F;
    n := numgens R;
    z := symbol z;
    S := CC[z_1..z_n];
    phi := map(S,R,toList(z_1..z_n)|(coordinates pt));
    F' := polySystem apply(equations F,f->phi(f));
    F'
    )

specializeSys (PolySystem,List) := List => (F,pts)->(
    R := ring F;
    n := numgens R;
    z := symbol z;
    S := CC[z_1..z_n];
    for pt in pts list (
	phi := map(S,R,toList(z_1..z_n)|(coordinates pt));
    	F' := polySystem apply(equations F,f->phi(f));
    	F'
	)
    )


--rewrites equations in k[parms][variables] without doublely-indexed variables.
rewriteEqns = method()
rewriteEqns (PolySystem) := PolySystem => F->(
    R := ring F;
    n := numgens R;
    m := numgens coefficientRing R;
    
    x := symbol x;
    a := symbol a;
    
    newR := CC[a_1..a_m][x_1..x_n];
    phi := map(newR,R,toList(x_1..x_n)|toList(a_1..a_m));
    polySystem apply(equations F,f->phi(f))
    )


--Determines if a set of Points corresponds to a permutation of the solutions
----over the base of a given Monodromy object (via alpha theory).
permutation = method()
permutation (Monodromy,List) := List => (M,S)->(
    numSolns := #(M#baseSolutions);
    P := apply(S,p->position(toList(0..numSolns-1),i->norm(2,point {coordinates M#baseSolutions#i - coordinates p}) < 1/(20*M#alphaConstants#i#2)));
    if (set P === set toList(0..numSolns-1)) then (
	P
	) else (
	error "--permutation: No permutation between given list of Points and base solutions."
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


--Counts the number of cycles of a permutation written in list form
numCycles = method()
numCycles (List) := ZZ => P->(
    cycles := 0;
    L := toList(0 .. #P-1);
    while (#L>0) do (
	i := first L;
	while (P#i != first L) do (
	    L = delete(P#i,L);
            i = P#i
	    );
	L = drop(L,1);
	cycles = cycles + 1
	);
    cycles
    )


--Computes the cycle type of a permutation given in list form.
cycleType = method()
cycleType (List) := String => P->(
    cType := {};
    L := toList(0 .. #P-1);
    while (#L>0) do (
	len := 1;
	i := first L;
	while (P#i != first L) do (
	    L = delete(P#i,L);
            i = P#i;
	    len = len + 1
	    );
	L = drop(L,1);
	cType = append(cType,len)
	);
    toString toSequence sort cType
    )


--------------------------------------
------------ Core Methods ------------
--------------------------------------

--Refines the solutions of a given monodromy object. Uses
----the same options as 'refine' given in "NumericalAlgebraicGeometry".
refine (Monodromy) := Monodromy => o->M->(
    F' := specializeSys(M#family,M#basePoint);
    newSolns := refine(F',M#baseSolutions,o);
    grp := pairs M#group;
    N := monodromy(M#family,M#basePoint,newSolns,BaseElements=>grp);
    N
    )


--Tracks solutions from base point to given points and back to base point.
----The resulting permutation of solutions is recorded.
monodromyLoop = method(Options=>{RefineSolns=>false,Verbosity=>0,stepIncreaseFactor=>2,tStep=>.01,maxCorrSteps=>1,CorrectorTolerance=>1e-6,tStepMin=>1e-8})
monodromyLoop (Monodromy,List) := Monodromy => o-> (M,pts)->(
    F := M#family;
    R := ring F;
    basePt := M#basePoint;
    baseSolns := M#baseSolutions;
	    
    pts = {basePt}|pts|{basePt};

    trackedSolns := baseSolns;
    
    --Specialized systems
    specSystems := specializeSys(F,pts);
    
    --Tracking base solutions along path defined by given points
    trackingOpts := new OptionTable from {
	stepIncreaseFactor=>o.stepIncreaseFactor,
	tStep=>o.tStep,
	maxCorrSteps=>o.maxCorrSteps,
	CorrectorTolerance=>o.CorrectorTolerance,
	tStepMin=>o.tStepMin
	};
    for i in 0..(#pts-2) do (
	trackedSolns = track(specSystems#i,specSystems#(i+1),trackedSolns,trackingOpts);
	if any(trackedSolns,p->status p != Regular) then error "--monodromyLoop: Tracking failure"
	);
    
    --Refine solutions 
    if (o.RefineSolns) then (
	trackedSolns = refine(specSystems#0,trackedSolns,Bits=>65)
	);
    
    --Checking that tracked solutions correspond to permutation of base solutions
    try(perm := permutation(M,trackedSolns)) else (
	if (o.Verbosity>0) then print("--monodromyLoop: No permutation from tracked solutions");
	if (o.Verbosity>1) then error "--monodromyLoop: Break in function";
	return M
	);
    if (o.Verbosity>0) then print printPermutation perm;
    
    --Adjoin permutation to the list of permutations
    if M#group#?perm then (
	M#group#perm = M#group#perm + 1
	) else (
	M#group#perm = 1
	);
    M
    )

--Tracks solutions from base point along a random triangular path. The 
----resulting permutation of solutions is recorded.
monodromyLoop (Monodromy) := Monodromy => o-> M->(
    n := numgens coefficientRing ring (M#family);
    r := norm(2,M#basePoint);
    pts := apply(2,i->point {apply(n,j->5*(random(CC)-random(CC)) )});
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
--
--Monte carlo solutions instead of tracking along one path
changeBasePoint = method(Options=>{MonteCarlo=>10,stepIncreaseFactor=>3,tStep=>.001,maxCorrSteps=>3,CorrectorTolerance=>1e-6})
changeBasePoint (Monodromy,Point) := Monodromy => o->(M,bp)->(
    F := M#family;
    basePt := M#basePoint;
    baseSolns := M#baseSolutions;
    grp := pairs M#group;
    
    trackingOpts := new OptionTable from {
	gamma=>random(RR)+ii*random(RR),
	stepIncreaseFactor=>o.stepIncreaseFactor,
	tStep=>o.tStep,
	maxCorrSteps=>o.maxCorrSteps,
	CorrectorTolerance=>o.CorrectorTolerance
	};
    specSystems := specializeSys(F,{basePt,bp});
    
    trackedSolns := track(specSystems#0,specSystems#1,baseSolns,trackingOpts);
    newBaseSolns := select(trackedSolns,p->status p == Regular);
    COUNT := 1;
    while (#newBaseSolns != #baseSolns) and (COUNT < o.MonteCarlo) do (
	trackedSolns = track(specSystems#0,specSystems#1,baseSolns,trackingOpts++{gamma=>random(CC)});
	trackedSolns = select(trackedSolns,p->status p == Regular);
	additionalSolns := select(trackedSolns,p->not any(newBaseSolns,q->p==q));
	newBaseSolns = newBaseSolns|additionalSolns;
	COUNT = COUNT + 1
	);
    
    if (#newBaseSolns === #baseSolns) then (
	monodromy(F,bp,newBaseSolns,BaseElements=>grp)
	) else (
	error "--changeBasePoint: Couldn't find all solutions over new base."
	)
    )


--Writes a Monodromy object to a file.
----Fix to save alphaConstants
saveMonodromy = method()
saveMonodromy (Monodromy,String) := Nothing => (M,s)->(
    F := M#family;
    R := ring F;
    C := coefficientRing R;
    s << toString gens R << endl << toString gens C << endl << toExternalString (M#family) << endl << toExternalString (M#basePoint) << endl << toExternalString (M#baseSolutions) << endl << toString (M#alphaConstants) << endl << toString (pairs M#group) << endl << toString (apply(keys M#cache,h->h=>M#cache#h)) << close
    )


--Recovers a Monodromy object from a file.
--Fix to load alphaConstants
loadMonodromy = method()
loadMonodromy (String) := Monodromy => s->(
    L := lines get s;
    polyVars := value L#0;
    parms := value L#1;
    R := CC[parms][polyVars];
    F := value L#2;
    basePt := value L#3;
    baseSolns := value L#4;
    alpha := value L#5;
    grp := value L#6;
    monodromyCache := value L#7;
    
    M := monodromy(F,basePt,baseSolns,alpha,BaseElements=>grp);
    M
    )


--Prints tracked permutations in cycle notation.
printCycleTally = method()
printCycleTally (Monodromy) := Nothing => M->(
    T := M#group;
    scan(keys T,sigma-> print(printPermutation(sigma) | ": " | toString(T#sigma)))
    )


--Outputs a string that represents the group generated by tracked
----permutations. Output can be read into GAP.
printForGAP = method()
printForGAP (Monodromy) := String => M->(
    T := keys M#group;
    if (#T > 0) then (
	S := apply(T,sigma->printPermutation(sigma));
    	"Group("|fold(S,(s,t)->s|","|t)|")"
	) else (
	"Group()"
	)
    )


--Determines the size of the group generated by tracked permutations.
size (Monodromy) := M->(
    if (#(M#group) == 0) then return 0;
    inFile := temporaryFileName();
    outFile := temporaryFileName();
    grp := printForGAP M;
    inFile << "G := "|grp|";; Size(G);" << close; 
    run("(cat "|inFile|" | gap -q) 1>"|outFile);
    n := get outFile;
    removeFile inFile;
    removeFile outFile;
    value n
    )


--Determines whether the group generated by tracked permutations is
----a subgroup of the alternating group.
testAlternatingMonodromy = method()
testAlternatingMonodromy (Monodromy) := Boolean => M->(
    T := keys M#group;
    if (#T > 0) then (
    	m := length (first T);
    	all(T,L-> (m - numCycles L) % 2 === 0)
	) else (
	false
	)
    )


---------------------------------------
------------ Documentation ------------
---------------------------------------









end

restart
loadPackage("FanoMonodromy")

