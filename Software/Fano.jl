###########
##FANO.JL##
###########
#
# Thomas Yahl
# Thomasjyahl@tamu.edu
# Created 5/27/22
#
# Code for exploring Fano problems and their Galois groups
#

using IterTools
using DelimitedFiles
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


#Computes the number of solutions to a finite Fano problem by the formula of Debarre and Manivel
function FanoNumSolns(r::Int,n::Int,degs::Vector{Int})
    if !isFiniteFano(r,n,degs)
        error("Not a finite Fano problem")
    end
    
    @polyvar x[1:r+1]
    degreePartitions = foldl(vcat,[map(S->subsetToPartition(S,d),subsets(1:r+d,r)) for d in degs])
    Q = prod(P->sum(P.*x),degreePartitions)
    V = prod(i->prod(j->x[i]-x[j],i+1:r+1),1:r)
    m = prod(i->x[i]^(n+1-i),1:r+1)
    coefficient(Q*V,m)
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

    #System(eqns,variables=y,parameters=p)
    (eqns,y,p)
end


#Computes a complex rational system with a complex rational simple double root (with high probability)
function FanoSimpleSingularSystem(r::Int,n::Int,degs::Vector{Int},singSoln::Vector{Complex{Rational{BigInt}}},tangVec::Vector{Complex{Rational{BigInt}}},targetCoeffs::Vector{ComplexF64})
    #create equations
    (F,y,p) = FanoEquations(r,n,degs)

    #linear constraints for singSoln to be a solution
    evalEqns = subs.(F,y=>singSoln)

    #linear constraints for singSoln to be a simple multiple root with tangent direction tangVec
    J = transpose(differentiate.(transpose(F),y))
    jacEqns = subs.(J,y=>singSoln)*tangVec

    #created for exact computation of matrix kernel
    (R,x) = PolynomialRing(QQ,"x")
    (QQi,ii) = NumberField(x^2+1,"I")

    #create matrix from linear constraints evalEqns and jacEqns, convert to Nemo matrix for exact nullspace
    A = transpose(hcat(map(f->coefficient.(f,p),vcat(evalEqns,jacEqns))...))
    A = map(z->real(z)+ii*imag(z),A)
    (numrows,numcols) = size(A)
    M = MatrixSpace(QQi,numrows,numcols)

    #compute exact nullspace, transform back to Complex{Rational{BigInt}} matrix
    K = nullspace(M(A))[2]
    K = map(z->Rational(coeff(z,0))+im*Rational(coeff(z,1)),Matrix(K))

    #find system with coefficients closest to targetCoeffs vector, round leastSquareSoln to small rational numbers
    leastSquareSoln = K \ targetCoeffs
    leastSquareSoln = map(z->BigInt(floor(10^4*real(z)))//10^4+im*BigInt(floor(10^4*imag(z)))//10^4,leastSquareSoln)
    systemCoeffs = K*leastSquareSoln
    
    #return coefficients, singular solution, and tangent vector
    (systemCoeffs,singSoln,tangVec)
end


function saveStartData(r::Int,n::Int,degs::Vector{Int},systemCoeffs::Vector{ComplexF64},nonsingSolns::Vector{Vector{ComplexF64}})
    #create output files
    FanoFile = string(r)*"_"*string(n)*"_"*foldl((a,b)->a*"_"*b,map(string,degs))
    sysCoeffsFile = open("../StartSystems/juliaData/"*FanoFile*"_coeffs.txt","w")
    nonsingSolnsFile = open("../StartSystems/juliaData/"*FanoFile*"_solns.txt","w")

    #write data to files
    writedlm(sysCoeffsFile,systemCoeffs)
    writedlm(nonsingSolnsFile,nonsingSolns)

    #close files
    close(sysCoeffsFile)
    close(nonsingSolnsFile)
end


function loadStartData(r::Int,n::Int,degs::Vector{Int})
    #find data files
    FanoFile = string(r)*"_"*string(n)*"_"*foldl((a,b)->a*"_"*b,map(string,degs))
    sysCoeffsFile = open("../StartSystems/juliaData/"*FanoFile*"_coeffs.txt","r")
    nonsingSolnsFile = open("../StartSystems/juliaData/"*FanoFile*"_solns.txt","r")
    
    #read from data files
    systemCoeffsData = readdlm(sysCoeffsFile,'\t',ComplexF64,'\n')
    nonsingSolnsData = readdlm(nonsingSolnsFile,'\t',ComplexF64,'\n')

    #convert data from string
    systemCoeffs = systemCoeffsData[:,1]
    numNonsingSolns = FanoNumSolns(r,n,degs)
    nonsingSolns = [nonsingSolnsData[i,:] for i in 1:numNonsingSolns]

    #output
    (systemCoeffs,nonsingSolns)
end



function saveData(r::Int,n::Int,degs::Vector{Int},systemCoeffs::Vector{Complex{Rational{BigInt}}},singSoln::Vector{Complex{Rational{BigInt}}},tangVec::Vector{Complex{Rational{BigInt}}},nonsingSolns::Vector{Vector{ComplexF64}})
    #create output files
    FanoFolder = string(r)*"_"*string(n)*"_"*foldl((a,b)->a*"_"*b,map(string,degs))*"/"
    sysCoeffsFile = open(FanoFolder*"juliaData/systemCoefficients.txt","w")
    singSolnFile = open(FanoFolder*"juliaData/singSoln.txt","w")
    tangVecFile = open(FanoFolder*"juliaData/tangVec.txt","w")
    nonsingSolnsFile = open(FanoFolder*"juliaData/nonsingSolns.txt","w")

    #convert data to write
    systemCoeffsToWrite = vcat(map(z->[numerator(real(z)), denominator(real(z)), numerator(imag(z)), denominator(imag(z))],systemCoeffs)...)
    singSolnToWrite = vcat(map(z->[numerator(real(z)), denominator(real(z)), numerator(imag(z)), denominator(imag(z))],singSoln)...)
    tangVecToWrite = vcat(map(z->[numerator(real(z)), denominator(real(z)), numerator(imag(z)), denominator(imag(z))],tangVec)...)
    
    #write data to files
    writedlm(sysCoeffsFile,systemCoeffsToWrite)
    writedlm(singSolnFile,singSolnToWrite)
    writedlm(tangVecFile,tangVecToWrite)
    writedlm(nonsingSolnsFile,nonsingSolns)

    #close files
    close(sysCoeffsFile)
    close(singSolnFile)
    close(tangVecFile)
    close(nonsingSolnsFile)
end

function loadData(r::Int,n::Int,degs::Vector{Int})
    #find data files
    FanoFolder = string(r)*"_"*string(n)*"_"*foldl((a,b)->a*"_"*b,map(string,degs))*"/"
    sysCoeffsFile = open(FanoFolder*"juliaData/systemCoefficients.txt","r")
    singSolnFile = open(FanoFolder*"juliaData/singSoln.txt","r")
    tangVecFile = open(FanoFolder*"juliaData/tangVec.txt","r")
    nonsingSolnsFile = open(FanoFolder*"juliaData/nonsingSolns.txt","r")
    
    #read from data files
    systemCoeffsData = readdlm(sysCoeffsFile,Int64)
    singSolnData = readdlm(singSolnFile,Int64)
    tangVecData = readdlm(tangVecFile,Int64)
    nonsingSolnsData = readdlm(nonsingSolnsFile,'\t',ComplexF64,'\n')

    #convert data from string
    numParams = sum(binomial.(degs.+n,n))
    numVars = (r+1)*(n-r)
    numNonsingSolns = FanoNumSolns(r,n,degs)-2
    
    systemCoeffs = [systemCoeffsData[4*i+1]//systemCoeffsData[4*i+2] + im*systemCoeffsData[4*i+3]//systemCoeffsData[4*i+4] for i in 0:numParams-1]
    singSoln = [singSolnData[4*i+1]//singSolnData[4*i+2] + im*singSolnData[4*i+3]//singSolnData[4*i+4] for i in 0:numVars-1]
    tangVec = [tangVecData[4*i+1]//tangVecData[4*i+2] + im*tangVecData[4*i+3]//tangVecData[4*i+4] for i in 0:numVars-1]
    nonsingSolns = [nonsingSolnsData[i,:] for i in 1:numNonsingSolns]

    #output
    (systemCoeffs,singSoln,tangVec,nonsingSolns)
end
