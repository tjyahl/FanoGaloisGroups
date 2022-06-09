#######################
##CERTIFYSOLUTIONS.JL##
#######################
#
# Thomas Yahl
# Thomasjyahl@tamu.edu
# Created 6/7/22
#
# Code for certifying data given in specified data folders
#

include("../Software/Fano.jl")

FanoProblems = [(1,7,[2,2,2,2]),(1,6,[2,2,3]),(2,8,[2,2,2]),(1,5,[3,3]),(1,5,[2,4]),(1,10,[2,2,2,2,2,2]),(1,9,[2,2,2,2,3]),(2,10,[2,2,2,2]),(1,8,[2,2,3,3])]

for FanoData in FanoProblems
    ####################
    ##LOAD SYSTEM DATA##
    ####################
    #r is the dimension of the sought linear spaces
    #n is the dimension of the ambient space
    #degs is the list of degrees of generic homogeneous polynomials
    (r,n,degs) = FanoData
    FanoDegree = FanoNumSolns(r,n,degs)
    FanoFolder = string(r)*"_"*string(n)*"_"*foldl((a,b)->a*"_"*b,map(string,degs))*"/"
    print("Checking Fano problem "*string(FanoData)*"\n")
    
    #Load system data and create system
    #F is the family of systems computing r-planes on the intersection of polynomials f₁,..,fₛ of degrees d₁,..,dₛ.
    #The polynomials f₁,..,fₛ are written in graded lex monomial order with coefficients given by p₁,..,pₖ. The p₁,..,pₖ are then the parameters of the family F.
    #systemCoeffs is the choice of parameters p₁,..,pₖ for a specific Fano problem, giving a system of polynomial equations F̄.
    #singSoln is the singular solution to F̄.
    #tangVec is a tangent vector at the singular solution to F̄.
    #nonsingSolns is the list of smooth solutions to F̄.
    (F,y,p) = FanoEquations(r,n,degs)
    (systemCoeffs,singSoln,tangVec,nonsingSolns) = loadData(r,n,degs)
    F̄ = subs.(F,p=>systemCoeffs)
    print("Data for singular system F̄ loaded\n")
    
    #############################
    ##VERIFY SIMPLE DOUBLE ROOT##
    #############################
    #Check that singSoln is a solution of F̄
    #Note: the convert statement is needed for type checking
    evaluatedSystem = subs.(F̄,y=>singSoln)
    if (evaluatedSystem == convert(typeof(evaluatedSystem),zeros((r+1)*(n-r))))
        print("F̄(singSoln) = 0. Given exact point is exact zero of system\n")
    else
        print("Error: Given exact point is NOT exact zero of system\n")
    end
    
    #Check that tangVec ∈  ker DF̄(singSoln) so that singSoln is a singular solution of F̄
    #J is the jacobian of the system
    J = transpose(differentiate.(transpose(F̄),y))
    JEvaluated = convert(Matrix{Complex{Rational{Int64}}},subs.(J,y=>singSoln))
    if (JEvaluated*tangVec == convert(typeof(J*tangVec),zeros((r+1)*(n-r))))
        print("DF̄(singSoln)*tangVec = 0. Given exact point is singular solution of system\n")
    else
        print("Error: Issue showing exact point is singular solution of system\n")
    end
    
    #Check that singSoln is a simple double root in the sense of Shub
    #This amounts to checking:
    #1) rank DF̄(singSoln) = (r+1)*(n-r)-1 or equivalently, dim ker DF̄(singSoln) = 1
    #2) D²F̄(singSoln)(tangVec,tangVec) ∉  Im DF̄(singSoln)
    #Note: julia's linear algebra package does not do exact computations, Nemo is necessary for this purpose. The number field QQi is the rational numbers adjoin the imaginary unit i, written as ii. Nemo works over this number ring.
    (R,x) = PolynomialRing(QQ,"x")
    (QQi,ii) = NumberField(x^2+1,"I")
    M = MatrixSpace(QQi,(r+1)*(n-r),(r+1)*(n-r))
    
    #Checking 1 above
    JNemo = M(map(z->real(z)+ii*imag(z),JEvaluated))
    (nullity,v) = nullspace(JNemo)
    if (nullity == 1) 
        print("dim ker DF̄(singSoln) = 1. Tangent space at exact singular solution has dimension 1\n")
    else
        print("Error: Tangent space does not have dimension 1\n")
    end

    #Checking 2 above
    H = [subs.(differentiate(f,y,2),y=>singSoln) for f in F̄]
    HContracted = [transpose(tangVec)*hes*tangVec for hes in H]
    
    N = MatrixSpace(QQi,(r+1)*(n-r),1)
    HContractedNemo = N(map(z->real(convert(Complex{Rational{Int64}},z))+ii*imag(convert(Complex{Rational{Int64}},z)),HContracted))

    if (rank(hcat(JNemo,HContractedNemo)) == (r+1)*(n-r))
        print("D²F̄(singSoln)(tangVec,tangVec) ∉  Im DF̄(singSoln) so singSoln is a simple double root\n")
    else
        print("Error: Second condition failed\n")
    end

    ###############################
    ##CERTIFY REMAINING SOLUTIONS##
    ###############################
    #Certify remaining solutions
    print("Degree of Fano problem: "*string(FanoDegree)*"\n")
    @time C = certify(F̄,nonsingSolns)
    certifiedSolns = filter(c->c.certified,C.certificates)
    print("Number of certified solutions: "*string(length(certifiedSolns))*"\n")
    if all(c->!(singSoln in c.I),certifiedSolns)
        print("singSoln is not contained in any certification interval\n\n")
    else
        print("singSoln is contained in some certification interval\n\n")
    end 
end
