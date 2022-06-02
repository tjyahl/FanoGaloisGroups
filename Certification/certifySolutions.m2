--------------------------------------------------------------------
---- CODE TO CERTIFY DISTINCT SOLUTIONS FOR GIVEN FANO PROBLEMS ----
--------------------------------------------------------------------

restart
loadPackage("FanoProblems",FileName=>"../Software/FanoProblems.m2",Reload=>true)
loadPackage("alphaCertified",FileName=>"../Software/alphaCertified.m2",Reload=>true)

--Add or remove problems as deemed necessary.
--1_7_2_2_2_2, 1_6_2_2_3, 1_5_3_3 ~3hrs
--1_5_2_4 ~6hrs
--2_8_2_2_2 ~28hrs
for FanoProblem in {"1_7_2_2_2_2","1_6_2_2_3","2_8_2_2_2","1_5_3_3","1_5_2_4"} do (
    print("Starting Fano problem "|FanoProblem);
    systemData = (lines get (FanoProblem|"/systemData.txt"))/value;
    nonsingSolns = value get (FanoProblem|"/nonsingSolns.txt");
    
    K = systemData#1;
    eqns = FanoEquations(systemData#0,BaseRing=>K);
    R = ring eqns;
    C = K[z_1..z_(numgens R)];
    
    basePt = systemData#2;
    specializeSystem = map(C,R,(gens C)|basePt);
    F = polySystem apply(equations eqns,f->specializeSystem(f));
    
    singSoln = apply(systemData#3,c->{coefficient(1_K,sub(c,K)),coefficient(K_0,sub(c,K))});
    solns = {singSoln}|nonsingSolns;
    
    --SHOW SINGULAR SOLUTION IS ISOLATED
    ----MAKE M2 COMPUTE QUOTIENTS IN GR??
    ----USE CONTAINMENT CODE FOR BELOW
    phi = map(K,C,apply(singSoln,l->l#0 + I*l#1));
    jac = phi jacobian F;
    print("Dimension of ambient space: "|toString(numgens R));
    print("Rank of Jacobian DF at singular solution x: "|toString(rank jac));
    print("Choosing tangent vector v at singular solution x");
    tangVec = (gens ker jac)_{0};
    w = fold(apply(equations F,f->(transpose tangVec)*(phi jacobian transpose jacobian f)*tangVec),(a,b)->a||b);
    print("Does D^2F(x)(v,v) % Im DF(x) == 0: "|toString(zero(w % image phi jacobian F)));
    
    
    print("Certifying solutions");
    --CHECK IF THIS NEEDS TO BE FIXED
    elapsedTime (S = certifySolutions(F,solns));
    print(#S|" certified solutions of singular system");
    print((FanoNumSolns systemData#0)|" expected solutions")
    )



