-----------------------
-----------------------
-- alphaCertified.m2 --
-----------------------
-----------------------

newPackage(
    "alphaCertified",
    Version=>"0.1.0",
    Authors=>{{
	Name=>"Thomas Yahl",
	Email=>"Thomasjyahl@tamu.edu",
	Homepage=>"https://math.tamu.edu/~thomasjyahl"
	}},
    Headline=>"Methods to utilize alphaCertified to certify solutions to Gaussian rational systems with Gaussian rational approximate solutions",
    PackageImports=>{"NumericalAlgebraicGeometry"},
    PackageExports=>{"NumericalAlgebraicGeometry"},
    DebuggingMode=>true
    )

export{
    -------------------------
    --Certification methods--
    -------------------------
    "certifySolutions",
    
    -------------------
    --Utility methods--
    -------------------
    "toGaussianRational",
    
    -----------
    --Options--
    -----------
    "alphaCertifiedExe",
    "alphaCertifiedConfig",
    "ExtraPrecision"
    }


-----------------------------------------
------------ Utility Methods ------------
-----------------------------------------

--Method for approximating complex solutions by Gaussian rational solutions for alphaCertified.
toGaussianRational = method(Options=>{ExtraPrecision=>0})
toGaussianRational (CC) := List => o->z->(
    {lift(realPart (10^(o.ExtraPrecision)*z),QQ)/10^(o.ExtraPrecision),lift(imaginaryPart (10^(o.ExtraPrecision)*z),QQ)/10^(o.ExtraPrecision)}
    )

toGaussianRational (Point) := List => o->p->(
    apply(coordinates p,z->toGaussianRational(z,o))
    )

toGaussianRational (List) := List => o->L->(
    apply(L,p->toGaussianRational(p,o))
    )


--Method for creating input files for AlphaCertified
----Intended for Gaussian rationals systems and solutions
writeAlphaCertifiedFile = method()
writeAlphaCertifiedFile (PolySystem) := String => F->(
    temp := temporaryFileName();
    file := openOutAppend(temp);
    C := coefficientRing ring F;
    n := F#NumberOfPolys;
    file << toString(n)|" "|toString(n) << endl << endl;
    
    eqns := equations F;
    for f in eqns do (
	exps := exponents f;
	coeffs := apply(flatten entries sub(last coefficients f,C),c-> flatten entries last coefficients(c,Monomials=>flatten entries basis C));
	numTerms := #exps;
	
	file << toString(numTerms) << endl;
	for i from 0 to numTerms - 1 do (
	    file << fold((a,b)->a|toString(b)|" ","",exps#i|coeffs#i) << endl
	    )
	);
    file << close;
    
    temp
    )

writeAlphaCertifiedFile (List) := Nothing => solns->(
    temp := temporaryFileName();
    file := openOutAppend(temp);
    numSolns := #solns;
    n := #(first solns);
    
    file << toString(numSolns) << endl;
    for s in solns do (
	file << endl;
	for z in s do (
	    file << fold((a,b)->a|toString(b)|" ","",z) << endl
	    )
	);
    file << close;
    
    temp
    )


-----------------------------------------------
------------ Certification Methods ------------
-----------------------------------------------

--Outputs the certified distinct approximate solutions from AlphaCertified with input from F and solns.
----Intended for Gaussian rational systems and solutions
certifySolutions = method(Options=>{alphaCertifiedExe=>"../Software/alphaCertified", alphaCertifiedConfig=>"../Software/settings"})
certifySolutions (PolySystem, List) := List => o->(F,solns)->(
    f1 := writeAlphaCertifiedFile(F);
    f2 := writeAlphaCertifiedFile(solns);
    run(o.alphaCertifiedExe|" "|f1|" "|f2|" "|o.alphaCertifiedConfig|" > ACOUTPUT");
    s := lines get "distinctSolns";
    run("rm ACOUTPUT");
    run("rm approxSolns");
    run("rm constantValues");
    run("rm distinctSolns");
    run("rm isApproxSoln");
    run("rm isDistinctSoln");
    run("rm redundantSolns");
    run("rm refinedPoints");
    run("rm summary");
    run("rm unknownPoints");
    
    if (#s == 1) then return({});
    numSolns := value s#0;
    n := sub((#s-2)/numSolns - 1,ZZ);
    s = apply(delete("",drop(s,1)),s->apply(separate(" ",s),value));
    
    distinctSolns := table(numSolns,n,(i,j)->s#(n*i+j));
    distinctSolns
    )


---------------------------------------
------------ Documentation ------------
---------------------------------------

beginDocumentation()

undocumented {[toGaussianRational,ExtraPrecision],[certifySolutions,alphaCertifiedExe],[certifySolutions,alphaCertifiedConfig]}

document {
    Key => alphaCertified,
    Headline => "Methods for certifying distinct Gaussian rational solutions to Gaussian rational systems",
    "Provides methods for approximating complex systems and solutions by Gaussian rational systems and solutions, as well as certifying solutions to these approximations."
    }

document {
    Key => toGaussianRational,
    Headline => "Approximate complex data by Gaussian rational data",
    Usage => "toGaussianRational(z)
            toGaussianRational(p)
	    toGaussianRational(L)",
    Inputs => {
        CC => "z" => {"a complex number"},
        Point => "p" => {"a complex point"},
        List => "L" => {"a list of complex points"}
        },
    Outputs => { List => {"an approximation of the complex number by a Gaussian rational data"},
	List => {"an approximation of the complex point by a Gaussian rational data"},
        List => {"an approximation of the complex points by Gaussian rational data"}
	},
    PARA {""}
    }

document {
    Key => certifySolutions,
    Headline => "Certify distinct Gaussian rational solutions to Gaussian rational systems",
    Usage => "certifySolutions(F,solns)",
    Inputs => {
	PolySystem => "F" => {"a Gaussian rational polynomial system"},
        List => "solns" => {"a list of approximate Gaussian rational solutions"}
        },
    Outputs => { List => "a list of input solutions certified by alphaCertified to be distinct approximate solutions"},
    PARA {""}
    }


end

restart
loadPackage("alphaCertified")


