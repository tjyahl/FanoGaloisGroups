---------------------------------------------------------------------------------
---- CODE TO GENERATE (GAUSSIAN) RATIONAL SYSTEMS WITH A SINGLE DOUBLE POINT ----
---------------------------------------------------------------------------------

restart
loadPackage("FanoProblems",FileName=>"../Software/FanoProblems.m2",Reload=>true)
loadPackage("alphaCertified",FileName=>"../Software/alphaCertified.m2",Reload=>true)

--Choice of Fano problem
FanoData = (1,6,(2,2,3))


--generate defining equations
GR = QQ[I]/(I^2+1);
F = FanoEquations(FanoData,BaseRing=>GR)

R = ring F
C = coefficientRing R
n = numgens R
m = numgens C


--Choice of singular solution
singSoln = apply(n,i->(random(-3,3)*random(100)/100+I*random(-3,3)*random(100)/100));


--linear constraints on coefficients so that singSoln is a solution of F
----phi substitutes singSoln for variables
phi = map(C,R,singSoln);
K = sub(transpose last coefficients(transpose phi(F.PolyMap),Monomials=>gens C),GR);


--Dumb trick, choice of hyperplane containing derivatives so that singSoln is singular solution
----choose kervec creatively?
J = phi jacobian F;
kervec = matrix table(1,n,(i,j)->I^(random(4))*(random(100)/10 + I*random(100)/10))
K' = sub(transpose last coefficients(kervec*J,Monomials=>gens C),GR);


--Choice of coefficients. Column span of kergens is the set of usable coefficients.
----use least squares to determine coefficients closest to a reasonable vector of coefficients, v.
----choose vector of coefficients so that solutions have moderate size/condition numbers?
kergens = gens kernel (K||K');
v = transpose matrix {apply(m,i->random(CC))};

CCmap = map(CC,GR,{ii+0.0});
xLS = solve(CCmap(kergens),v,ClosestFit=>true,MaximalRank=>true);

xLSGR = matrix table(numgens target xLS,1,(i,j)->lift(realPart xLS_(i,j),QQ) + I*lift(imaginaryPart xLS_(i,j),QQ));
singBase = flatten entries transpose (kergens*xLSGR);


--Check chosen system is singular.
----for base system substitutes singSoln for variables and singBase for parameters
----shows kervec is in the kernel of the jacobian evaluated at singSoln
singTest = map(GR,R,singSoln|singBase);
singTest(F.PolyMap)
kervec*(singTest(jacobian F))


--Complexify singular base and solution for numerical methods
CCsingBase = singBase/CCmap
CCsingSoln = singSoln/CCmap

max apply(CCsingBase,abs)


--Start system data
FanoProblem = fold((a,b)->a|"_"|b,(splice FanoData)/toString)
problemFile = "../StartSystems/"|FanoProblem|".txt"
try (M = loadMonodromy(problemFile)) else (M = FanoMonodromy(FanoData,Verbose=>true))


--Systems needed for numerical methods
----start data from StartSystems folder
----NAG doesn't track with higher precision, use Bertini and a higher precision ring as needed
----NAG will refine to arbitrary precision (hence higher precision for trackingRing)
p = 200
GRRing = GR[y_1..y_n]
trackingRing = CC_p[z_1..z_n]

GRsingSystem = polySystem (map(GRRing,R,(gens GRRing)|singBase))(F.PolyMap);
CCsingSystem = polySystem (map(trackingRing,GRRing,(gens trackingRing)|{toCC(p,ii+0.0)}))(GRsingSystem.PolyMap);

startSystem = polySystem (map(trackingRing,R,(gens trackingRing)|apply(coordinates M#basePoint,z->toCC(p,z))|{toCC(p,ii+0.0)}))(F.PolyMap);
startSolns = M#baseSolutions;


--Tracking solutions from startSystem to obtain solutions to our singular system
----there are many options that can be changed here if needed
----moving the start system via N adds a bit of randomness into the solving, not necessary
refinementFails = {}
certifiedSolns = {}
count = 1
while (#certifiedSolns<#M#baseSolutions-2) do (
    print("starting count: "|toString(count));
    --track solutions from startSystem to CCSystem
    trackedSolns = track(startSystem,CCsingSystem,startSolns,gamma=>random(CC),stepIncreaseFactor=>1.5,tStep=>1e-10,maxCorrSteps=>20,CorrectorTolerance=>1e-6,tStepMin=>1e-30);
    print("tracking completed");
    
    --checking which, if any, are new solutions
    regularSolns = select(trackedSolns,s->status s == Regular);
    regularSolns = select(regularSolns,z->(
	    GRz = apply(coordinates z,a->{lift(realPart a,QQ),lift(imaginaryPart a,QQ)});
	    all(certifiedSolns,s->sum apply(#GRz,i->sum(GRz#i-s#i,j->j^2)) > 1e-3) and (norm(2,coordinates z - CCsingSoln) > 1e-3)
	    ));
    if (#regularSolns == 0) then (print("no potential new solutions"); print(""); count = count+1; continue);
    print(toString(#regularSolns)|" potential new solutions");
    
    --checking if there are any duplicate solutions from tracking
    temp = {};
    for z in regularSolns do (
	if all(temp,s->norm(2,coordinates s - coordinates z) > 1e-3) then temp = temp|{z};
	);
    regularSolns = temp;
    print(toString(#regularSolns)|" distinct potential new solutions found");
    
    --certify new solutions by a refinement filter
        for z in regularSolns do (
	print("certifying "|toString(position(regularSolns,s->s==z)+1)|"-th solution out of "|toString(#regularSolns));
	print(coordinates z);
      	refineCount = 0;
	while (refineCount < 30) do (
	    z = first refine(CCsingSystem,{z},Bits=>50+5*refineCount);
	    if (status z == RefinementFailure) then (print("refinement failed"); print("#L = "|toString(#L)); refinementFails = refinementFails|{z}; print(""); break);
	    refineCount = refineCount+1;
	    print("refined "|toString(refineCount)|" times");
	    try(certifiedSolns = certifiedSolns|certifySolutions(GRsingSystem,{toGaussianRational(z,ExtraPrecision=>5*(refineCount-1))})) else continue;
	    print("new solution certified");
	    print("#certifiedSolns = "|toString(#certifiedSolns));
	    break
	    );
	);
    
    count = count + 1;
    print(" ")
    )


--save all and check certifiedSolns are all certified distinct solutions (with singSoln)
----systemData.txt has Fano data, coefficientRing, singBase, AND singSoln
----nonsingSolns.txt contains all NONSINGULAR solns
systemDataFile = openOut(fanoProb|"/systemData.txt")
nonsingSolnsFile = openOut(fanoProb|"/nonsingSolns.txt")

systemDataFile << "(1,5,(2,4))" << endl << "QQ[I]/(I^2+1)" << endl << toString(singBase) << endl << toString(singSoln) << close;
nonsingSolnsFile << toString(certifiedSolns) << close;
