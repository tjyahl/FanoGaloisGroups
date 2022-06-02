###################
##FANOPROBLEMS.JL##
###################
#
# Thomas Yahl
# Thomasjyahl@tamu.edu
# Created 5/27/22
#
# Code for exploring Fano problems and their Galois groups
#

using IterTools
using MultivariatePolynomials
using HomotopyContinuation
using LinearAlgebra
using Nemo


#####################
##UTILITY FUNCTIONS##
#####################

function subsetToPartition(S,n)
    [if i==1
        S[1]-1
        elseif i==length(S)+1
        n+length(S)-last(S)
        else
        S[i]-S[i-1]-1
        end
    for i in 1:length(S)+1]
end



#############
##FUNCTIONS##
#############

#Checks whether the Fano problem is finite
function isFiniteFano(r::Int,n::Int,degs::Vector{Int})
    (n-r)*(r+1) == sum(binomial.(degs.+r,r))
end



#Computes the defining equations of the desired linear spaces in a coordinate chart of the appropriate Grassmanian
function FanoEquations(r::Int,n::Int,degs::Vector{Int})
    if !isFiniteFano(r,n,degs)
        error("Not a finite Fano problem")
    end
    s = length(degs)
    
    #equations for complete intersection of polynomials with given degrees
    numParameters = sum(binomial.(degs.+n,n))
    @polyvar x[1:n+1] t[1:r+1] y[1:(n-r)*(r+1)] p[1:numParameters]
    monoms = [transpose(HomotopyContinuation.monomials(x,d)) for d in degs]
    splits = vcat([0],accumulate(+,binomial.(degs.+n,n)))
    eqParams = [p[splits[i]+1:splits[i+1]] for i=1:s]
    F = monoms.*eqParams
    
    #linear space parameterization and substitution
    ℓ = vcat([sum(y[(r+1)*(i-1)+1:(r+1)*i].*t) for i=1:n-r],t)
    F̂ = subs(F,x=>ℓ)
    
    #pull t coefficients (painful in Julia?)
    eqns = [];
    for i=1:s
    	monoms = HomotopyContinuation.monomials(t,degs[i])
    	for m in monoms
	    append!(eqns,[coefficient(F̂[i],m,t)])
	end
    end    
    eqns = subs.(eqns,t=>ones(Int,r+1))

    (eqns,y,p)
end



#Computes the number of solutions to a finite Fano problem by the formula of Debarre and Manivel
function FanoNumSolns(r::Int,n::Int,degs::Vector{Int})
    @polyvar x[1:r+1]
    degreePartitions = foldl(vcat,[map(S->subsetToPartition(S,d),subsets(1:r+d,r)) for d in degs])
    Q = prod(P->sum(P.*x),degreePartitions)
    V = prod(i->prod(j->x[i]-x[j],i+1:r+1),1:r)
    m = prod(i->x[i]^(n+1-i),1:r+1)
    coefficient(Q*V,m)
end



#Computes a start system based off a given (exact) start solution
function FanoStartSystem(r::Int,n::Int,degs::Vector{Int},startSoln::Vector{Complex{Rational{BigInt}}})
    #create equations and linear constraints for startSoln to be a solution
    (F,y,p) = FanoEquations(r,n,degs)
    evalEqns = subs.(F,y=>startSoln)
    
    #created for exact computation of matrix kernel
    (R,x) = PolynomialRing(QQ,"x")
    (QQi,ii) = NumberField(x^2+1,"I")

    #create matrix from linear constraints, convert to Nemo matrix for exact nullspace
    A = transpose(hcat(map(f->coefficient.(f,p),evalEqns)...))
    A = map(z->real(z)+ii*imag(z),A)
    (numrows,numcols) = size(A)
    M = MatrixSpace(QQi,numrows,numcols)

    #compute exact nullspace, transform back to Complex{Rational{BigInt}} matrix
    K = nullspace(M(A))[2]
    K = map(z->Rational(coeff(z,0))+im*Rational(coeff(z,1)),Matrix(K))

    #find system with coefficients closest to all ones vector
    leastSquareSoln = K \ ones(sum(binomial.(degs.+n,n)))
    leastSquareSoln = map(z->Rational(real(z))+im*Rational(imag(z)),leastSquareSoln)
    coeffs = K*leastSquareSoln
    
    #substitute coefficients into system
    (F,y,p,coeffs,startSoln)
end

function FanoStartSystem(r::Int,n::Int,degs::Vector{Int})
    startSoln = [BigInt(rand(-10:10))//rand(1:10) + rand(-10:10)//rand(1:10)*im for i in 1:(n-r)*(r+1)]
    FanoStartSystem(r,n,degs,startSoln)
end



#Computes a complex rational system with a complex rational simple double root (with high probability)
function FanoSimpleSingularSystem(r::Int,n::Int,degs::Vector{Int},singSoln::Vector{Complex{Rational{BigInt}}},tangVec::Vector{Complex{Rational{BigInt}}},startSoln::Vector{Complex{Rational{BigInt}}})
    #create equations
    (F,y,p) = FanoEquations(r,n,degs)

    #linear constraints for singSoln to be a solution
    evalEqns = subs.(F,y=>singSoln)

    #linear constraints for singSoln to be a simple multiple root with tangent direction tangVec
    J = transpose(differentiate.(transpose(F),y))
    jacEqns = subs.(J,y=>singSoln)*tangVec

    #linear constraints for startSoln to be a solution (mostly for starting monodromy)
    evalStartEqns = subs.(F,y=>startSoln)

    #created for exact computation of matrix kernel
    (R,x) = PolynomialRing(QQ,"x")
    (QQi,ii) = NumberField(x^2+1,"I")

    #create matrix from linear constraints evalEqns and jacEqns, convert to Nemo matrix for exact nullspace
    A = transpose(hcat(map(f->coefficient.(f,p),vcat(evalEqns,jacEqns,evalStartEqns))...))
    A = map(z->real(z)+ii*imag(z),A)
    (numrows,numcols) = size(A)
    M = MatrixSpace(QQi,numrows,numcols)

    #compute exact nullspace, transform back to Complex{Rational{BigInt}} matrix
    K = nullspace(M(A))[2]
    K = map(z->Rational(coeff(z,0))+im*Rational(coeff(z,1)),Matrix(K))

    #find system with coefficients closest to all ones vector
    leastSquareSoln = K \ ones(sum(binomial.(degs.+n,n)))
    leastSquareSoln = map(z->Rational(real(z))+im*Rational(imag(z)),leastSquareSoln)
    coeffs = K*leastSquareSoln
    
    #substitute coefficients into system
    G = subs.(F,p=>coeffs)
    (F,y,p,coeffs,singSoln,tangVec,startSoln)
    
end

function FanoSimpleSingularSystem(r::Int,n::Int,degs::Vector{Int})
    singSoln = [BigInt(rand(-10:10))//rand(1:10) + rand(-10:10)//rand(1:10)*im for i in 1:(n-r)*(r+1)]
    
    tangVec = [BigInt(rand(-10:10))//rand(1:10) + rand(-10:10)//rand(1:10)*im for i in 1:(n-r)*(r+1)]
    startSoln = [BigInt(rand(-10:10))//rand(1:10) + rand(-10:10)//rand(1:10)*im for i in 1:(n-r)*(r+1)]
    FanoSimpleSingularSystem(r,n,degs,singSoln,tangVec,startSoln)
end





#Check that singSoln is singular solution and startSoln is solution as well. These should all be zero.
r = 1
n = 6
degs = [2,2,3]
(F,y,p,coeffs,singSoln,tangVec,startSoln) = FanoSimpleSingularSystem(r,n,degs);
G = subs.(F,p=>coeffs);
subs.(G,y=>startSoln)
subs.(G,y=>singSoln)
transpose(subs.(differentiate.(transpose(G),y),y=>singSoln))*tangVec

(F,y,p,coeffs,startSoln) = FanoStartSystem(r,n,degs);

#restrict to line in parameter space
@polyvar t
ℓ = coeffs.*(1-t)+[BigInt(rand(-5:5))//rand(1:5) + rand(-50:50)//rand(1:5)*im for i in 1:sum(binomial.(degs.+n,n))].*t;

F̂ = subs.(F,p=>ℓ);
monodromy_solve(F̂,[startSoln],[0],parameters=[t])













#Just don't use?
#Loads data from Certification folders
function valueQQ(s::AbstractString)
    s = replace(s,r"\(|\)|\{|\}| "=>"")
    L = split(s,"/")
    if length(L)==1
	parse(BigInt,L[1])
    else
	parse(BigInt,L[1])//parse(BigInt,L[2])
    end
end

function valueQQI(s::AbstractString)
    s = replace(s,r"\(|\)|\{|\}| "=>"")
    L = filter(!isempty,split(s,"I"))
    sum(a->(contains(a,'*') && !isempty(a)) ? valueQQ(replace(a,r"\*"=>""))*im : valueQQ(a),L) + 0//1*im
end

function loadFanoData(r::Int,n::Int,degs::Vector{Int})
    FanoProblem = "$(r)_$(n)_"*foldl((a,b)->a*"_"*b,map(string,degs))*"/"
    dataFile = open(FanoProblem*"systemData.txt")
    systemData = readlines(dataFile)
    close(dataFile)

    coeffs = valueQQI.(split(systemData[3],","))
    startSoln = valueQQI.(split(systemData[4],","))

    (F,y,p) = FanoEquations(r,n,degs)
    F̂ = subs.(F,p=>coeffs,y=>startSoln)

end






















