------------------------------------------------------------------------------
---- CODE TO COMPUTE AND WRITE DATA FOR SMALL LOOP AROUND SINGULAR SYSTEM ----
------------------------------------------------------------------------------

restart
loadPackage("FanoProblems",FileName=>"../Software/FanoProblems.m2",Reload=>true)

--Choice of Fano problem
FanoData = (1,10,(2,2,2,2,2,2))

--Reformating FanoData for saving info and loading pre-existing data if possible
FanoProblem = fold((a,b)->a|"_"|b,(splice FanoData)/toString)
problemFile = "../StartSystems/"|FanoProblem|".txt"
try (M = loadMonodromy(problemFile)) else (M = FanoMonodromy(FanoData,Verbose=>true))

--Generate system
F = FanoEquations(FanoData)
R = ring F;
C = coefficientRing R;
n = numgens R;
p = numgens C;

--Choice of singular solution
----K is the linear constraints on the coefficients for singSoln to be a solution of F
singSoln = point random(CC^1,CC^n);
phi = map(C,R,coordinates singSoln);
K = sub(transpose last coefficients(transpose phi(F.PolyMap),Monomials=>gens C),CC);

--Choice of linear dependence of derivatives at singSoln
----K' is the linear constraints on the coefficients for singSoln to be a singular solution (with specified dependence on derivatives)
kervec = random(CC^1,CC^n);
K' = sub(transpose last coefficients(kervec*(phi jacobian F),Monomials=>gens C),CC);

--Least squares method to choose coefficients closest to a specified vector of coefficients v
----for some reason "gens ker (K||K')" seems unstable? current code taken from M2 Core "kernel" method without "mingens".
kergens = syz gb(K||K',Syzygies=>true);
--v = random(CC^p,CC^1);
v = transpose matrix{coordinates M#basePoint};
xLS = solve(kergens,v,ClosestFit=>true,MaximalRank=>true);
singBase = flatten entries transpose (kergens*xLS);

max apply(singBase,abs);
norm(2,singBase - coordinates M#basePoint)



--scan for best
bestDist = norm(2,singBase - coordinates M#basePoint)
bestSingSoln = singSoln
bestKervec = kervec
bestSingBase = singBase
for i from 1 to 10000 do (
    print("checking "|toString(i)|"-th singular system");
    
    singSoln = point random(CC^1,CC^n);
    phi = map(C,R,coordinates singSoln);
    K = sub(transpose last coefficients(transpose phi(F.PolyMap),Monomials=>gens C),CC);
    
    kervec = random(CC^1,CC^n);
    K' = sub(transpose last coefficients(kervec*(phi jacobian F),Monomials=>gens C),CC);
    
    kergens = syz gb(K||K',Syzygies=>true);
    v = transpose matrix{coordinates M#basePoint};
    xLS = solve(kergens,v,ClosestFit=>true,MaximalRank=>true);
    singBase = flatten entries transpose (kergens*xLS);
    
    dist = norm(2,singBase - coordinates M#basePoint);
    print("dist = "|toString(dist));
    if (dist < bestDist) then (
	bestDist = dist;
	bestSingBase = singBase;
	bestSingSoln = singSoln;
	bestKervec = kervec;
	print("new best found")
	)
    )
singSoln = bestSingSoln;
singBase = bestSingBase;
kervec = bestKervec;


--Check chosen system is singular.
----for base system substitutes singSoln for variables and singBase for parameters
----shows kervec is in the kernel of the jacobian evaluated at singSoln and the jacobian is singular
singTest = map(CC,R,(coordinates singSoln)|singBase);
singTest(F.PolyMap);
kervec*(singTest(jacobian F));

--The affine line (1-t)(M#basePoint)+t(singBase) has the singular base at t=1. The points newBase, pt1, and pt2 form a small triangle around singBase.
----might need to play with eps. 
--eps = .5;
--print("eps = "|toString(eps));
--t0 = 1 - eps;
--newBase = point {(1-t0)*(coordinates M#basePoint) + t0*singBase};

--Change base point to newBase
----Option NumPaths does a Monte Carlo method to obtain solution over the new base by tracking along at most the given number of paths.
--print("changing base");
--tol = 1e-5;
--try(N = changeBasePoint(M,newBase,NumPaths=>150,Tolerance=>tol,Verbose=>true)) else (print("changing base failed"); continue);

--Track solutions along the small triangular loop described by newBase, pt1, and pt2. 
----Verbosity=>2 prints all errors and the obtained permutation
--print("tracking along loop");
--try(N = monodromyLoop(N,{pt1,pt2},Verbosity=>2)) else (print("monodromyLoop failed"); continue);
N = M
eps = .4
tol = 1e-5
while(eps > .001) do (
    eps = eps/2;
    print("eps = "|toString(eps));
    t0 = 1 - eps; 
    newBase = point {(1-t0)*(coordinates M#basePoint) + t0*singBase};
    
    print("changing base");
    try(N = changeBasePoint(N,newBase,NumPaths=>150,Tolerance=>tol,Verbose=>true)) else (print("changing base failed"); break);
    
    if (eps < .01) then (
	t1 = 1 - eps*exp(2*pi*ii/3);
    	t2 = 1 - eps*exp(4*pi*ii/3);
	pt1 = point {(1-t1)*(coordinates M#basePoint) + t1*singBase};
    	pt2 = point {(1-t2)*(coordinates M#basePoint) + t2*singBase};
	
	print("tracking along loop");
    	try(N = monodromyLoop(N,{pt1,pt2},Verbosity=>2)) else (print("monodromyLoop failed"); break);
	)
    );

--saveMonodromy(N,FanoProblem|"/monod.txt")
--myFile = FanoProblem|"/pts.txt"
--myFile << toExternalString({pt1,pt2}) << close


end





