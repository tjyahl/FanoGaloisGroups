--Code to certify Fano problems with a single double point
----(1,7,(2,2,2,2)): 512 solutions
----(1,6,(2,2,3)): 720 solutions
----(2,8,(2,2,2)): 1024 solutions
----(1,5,(3,3)): 1053 solutions
----(1,4,(5)): 2875 solutions
----
---------------------------------------------------------------------
---------------------------------------------------------------------

restart
loadPackage("FanoProblems",FileName=>"../Software/FanoProblems.m2",Reload=>true)

M = loadMonodromy("../Data/1_7_2_2_2_2.txt")
M = loadMonodromy("../Data/2_8_2_2_2.txt")

F = transpose (M#family).PolyMap
n = numgens ring F
p = numgens coefficientRing ring F
s = toString F;
GR = QQ[I]/(I^2+1);
C = GR[a_1..a_p];
R = C[x_1..x_n];
--r = "(GR"|substring(6,toString ring F);
--R = value r;
F = value s;

randmat = (n,m,h,scale)->matrix table(n,m,(i,j)->I^(random(4))*(random(h)/(h*scale)+I*random(h)/(h*scale)))

--Choice of point to be the singular solution
singSoln = apply(n,i->I^(random(4))*(random(10,100)*random(100)/100+I*random(10,100)*random(100)/100));
phi = map(C,R,singSoln);

K = sub(transpose last coefficients phi(F),C)

--Dumb trick
J = phi jacobian F;
D0 = sub(transpose last coefficients(randmat(1,n,100,1)*J,Monomials=>gens C),C);
--D0 = sub(transpose last coefficients((randmat(1,6,10,1)|matrix{{0,0,0,0}})*J,Monomials=>gens C),C);
--D0 = sub(transpose last coefficients(randmat(1,3,20,10)*J^{0,4,8},Monomials=>gens C),C);
--D0 = sub(transpose last coefficients(transpose (J*randmat(n,1,20,10)),Monomials=>gens C),C);

--least squares?
CCmap = map(CC,GR,{ii+0.0})
A = sub(gens kernel (K||D0),GR);
xLS = solve(CCmap(A),transpose matrix {toList(p:(1 + 0*ii))},ClosestFit=>true,MaximalRank=>true);
B = flatten entries transpose ((sub(gens kernel (K||D0),GR))*transpose matrix{apply(flatten entries xLS,z->lift(realPart z,QQ) + I*lift(imaginaryPart z,QQ))});

--check sing system is sing.
singBase = B;
f = map(GR,R,singSoln|singBase);
f(F)
det f(jacobian F)

--complexify singular soln
CCsingBase = singBase/CCmap
CCsingSoln = singSoln/CCmap


--scaling 
----first scale constant to zero
F' = transpose (M#family).PolyMap
R' = ring F'
S = CC[y_1..y_n]
rho = map(S,R',(gens S)|CCsingBase)
maxc = apply(entries transpose last coefficients rho(F'),l->max apply(l,z->abs sub(z,CC)))
--G = apply(10,i->(flatten entries rho(F'))#i/maxc#i)
G = flatten entries rho(F')
test = map(CC,S,CCsingSoln)
G/test

rho' = map(S,R',(gens S)|(coordinates M#basePoint))

--solve
--L = {}
--N = M
count = 1
while (#L<718) and (count < 10) do (
    print("starting count: "|toString(count));
    try(N = changeBasePoint(N,point {coordinates N#basePoint + apply(p,i->random(CC))})) else print("issues changing base point");
    rho' = map(S,R',(gens S)|(coordinates N#basePoint));
    solns = track(flatten entries rho'(F'),G,(N#baseSolutions)/coordinates,gamma=>random(CC),stepIncreaseFactor=>10,tStep=>1e-10,maxCorrSteps=>30,CorrectorTolerance=>1e-6,EndZoneFactor=>.05,tStepMin=>1e-30);
    print("tracking completed");
    
    --unique apply(solns,status)
    regularsolns = select(solns,s->status s == Regular);
    print(toString(#regularsolns)|" regular solutions");
    
    refinedsolns = refine(G,solns,ErrorTolerance=>1e-12);
    refinedsolns = select(refinedsolns,s->status s != RefinementFailure);
    print(toString(#refinedsolns)|" refined solutions");
    
    T = GR[z_1..z_10];
    delta = map(T,R,(gens T)|singBase);
    try(distinct = AlphaCertified(polySystem transpose delta(F),refinedsolns)) else (print("no distinct certified solutions"); continue);
    print(toString(#distinct)|" distinct solutions");
    
    print(toString(number(distinct,z->not any(L,s->norm(2,apply(#z,i->sum(z#i-s#i,j->j^2))) < 1e-6)))|" new certified solutions");
    L = L|select(distinct,z->not any(L,s->norm(2,apply(#z,i->sum(z#i-s#i,j->j^2))) < 1e-6));
    print("#L = "|toString(#L));
    --L is ongoing list
    count = count + 1;
    print(" ")
    )




count = 1;
while (#L<718) and (count < 10) do (
    print("starting count: "|toString(count));
    CCL = apply(L,s->apply(s,l->(first l)+ii*(last l)));
    p1 = map(S,R',(gens S)|(coordinates N#basePoint + apply(140,i->.5*random(CC))));
    p2 = map(S,R',(gens S)|(coordinates N#basePoint + apply(140,i->.5*random(CC))));
    
    trackedsolns = track(G,flatten entries p1(F'),CCL,gamma=>random(CC),stepIncreaseFactor=>5,tStep=>1e-10,maxCorrSteps=>30,CorrectorTolerance=>1e-6,EndZoneFactor=>.05,tStepMin=>1e-80);
    trackedsolns = select(trackedsolns,s->status s == Regular);
    print("first tracking done");
    
    trackedsolns = track(flatten entries p1(F'),flatten entries p2(F'),trackedsolns,gamma=>random(CC),stepIncreaseFactor=>5,tStep=>1e-10,maxCorrSteps=>30,CorrectorTolerance=>1e-6,EndZoneFactor=>.05,tStepMin=>1e-80);
    trackedsolns = select(trackedsolns,s->status s == Regular);
    print("second tracking done");
    
    trackedsolns = track(flatten entries p2(F'),G,trackedsolns,gamma=>random(CC),stepIncreaseFactor=>5,tStep=>1e-10,maxCorrSteps=>30,CorrectorTolerance=>1e-6,EndZoneFactor=>.05,tStepMin=>1e-80);
    trackedsolns = select(trackedsolns,s->status s == Regular);
    print("third tracking done");
    print("#trackedsolns = "|toString(#trackedsolns));
    
    refinedsolns = refine(G,trackedsolns,ErrorTolerance=>1e-35);
    refinedsolns = select(refinedsolns,s->status s != RefinementFailure);
    print(toString(#refinedsolns)|" refined solutions");
    
    T = GR[b_1..b_10];
    delta = map(T,R,(gens T)|singBase);
    try(distinct = AlphaCertified(polySystem transpose delta(F),refinedsolns)) else (print("no distinct certified solutions"); continue);
    print(toString(#distinct)|" distinct solutions");
    
    print(toString(number(distinct,z->not any(L,s->norm(2,apply(#z,i->sum(z#i-s#i,j->j^2))) < 1e-6)))|" new certified solutions");
    L = L|select(distinct,z->not any(L,s->norm(2,apply(#z,i->sum(z#i-s#i,j->j^2))) < 1e-6));
    print("#L = "|toString(#L));
    --L is ongoing list
    count = count + 1;
    print(" ")
    )
