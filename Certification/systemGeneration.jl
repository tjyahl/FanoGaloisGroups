#######################
##SYSTEMGENERATION.JL##
#######################
#
# Thomas Yahl
# Thomasjyahl@tamu.edu
# Created 5/27/22
#
# Code for creating Fano problems with a simple double solution
#

include("../Software/Fano.jl")

#######################################
#Choosing Fano problem

r = 1
n = 8
degs = [2,2,2,4]
(F,y,p) = FanoEquations(r,n,degs);
sys = System(F,variables=y,parameters=p)
numSolns = FanoNumSolns(r,n,degs)


#######################################
#create start system (if it doesn't exist)
FanoFile = string(r)*"_"*string(n)*"_"*foldl((a,b)->a*"_"*b,map(string,degs))
if !isfile("../StartSystems/juliaData/"*FanoFile*"_coeffs.txt") || !isfile("../StartSystems/juliaData/"*FanoFile*"_solns.txt")
    M = monodromy_solve(sys,target_solutions_count=numSolns)
    startCoeffs = M.parameters;
    startSolns = [s.solution for s in M.results];
    saveStartData(r,n,degs,startCoeffs,startSolns)
else 
    (startCoeffs,startSolns) = loadStartData(r,n,degs)
end


#######################################
#find close singular system, parameter homotopy to solve
singSoln = vcat([1],zeros(Complex{Rational{BigInt}},(r+1)*(n-r)-1))
singSoln = rand(-2:2,(r+1)*(n-r))+zeros(Complex{Rational{BigInt}},(r+1)*(n-r))
tangVec = rand(-5:5,(r+1)*(n-r))+zeros(Complex{Rational{BigInt}},(r+1)*(n-r))
(systemCoeffs,singSoln,tangVec) = FanoSimpleSingularSystem(r,n,degs,singSoln,tangVec,startCoeffs)
S = HomotopyContinuation.solve(sys,startSolns,start_parameters=startCoeffs,target_parameters=systemCoeffs)

#find remaining solutions
#use monodromy_solve to compute solutions over base
#need smooth solutions to start monodromy_solve, obtain via parameter homotopy or smart choice of coefficients
#best to use when parameter homotopy from method 1 fails and already have many solutions
nonsingSolns = [s.solution for s in S.path_results];
M = monodromy_solve(sys,nonsingSolns,systemCoeffs)


#######################################
#certify nonsingular solutions and singular solution disjoint from boxes
#@time C = certify(sys,S,systemCoeffs)
@time C = certify(sys,M)
nonsingCertified = filter(c->c.certified,C.certificates);
all(c->!(singSoln in c.I),nonsingCertified)


#######################################
#save data
nonsingSolns = [c.solution for c in nonsingCertified];
saveData(r,n,degs,systemCoeffs,singSoln,tangVec,nonsingSolns)


#######################################
#check data loads
(systemCoeffs,singSoln,tangVec,nonsingularSolns) = loadData(r,n,degs)
specSys = subs.(F,p=>systemCoeffs)
C = certify(specSys,nonsingSolns)
all(c->!(singSoln in c.I),C.certificates)



