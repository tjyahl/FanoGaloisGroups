------------------------------------------------------------------------
---- FILE FOR TRACKING SMALL LOOP AROUND SINGULAR SYSTEMS FROM DATA ----
------------------------------------------------------------------------

restart
loadPackage("FanoProblems",FileName=>"../Software/FanoProblems.m2")

FanoProblems = sort select(toList (set readDirectory "." - set {".",".."}),isDirectory)

for FanoFolder in FanoProblems do (
    print("Loading monodromy from folder "|FanoFolder);
    if not isSubset({"monod.txt","pts.txt"},readDirectory FanoFolder) then (print("No data for specified problem"); continue);
    try(M = loadMonodromy(FanoFolder|"/monod.txt")) else (print("Error loading data for Fano problem") continue);
    try(pts = value get (FanoFolder|"/pts.txt")) else (print("Error loading data for Fano problem") continue);
    if isWellDefined M then (
	print("Computing monodromy loop for folder "|FanoFolder);
	s = elapsedTime monodromyLoop(M,pts,Verbosity=>2);
	) else (
	print("Monodromy object not well-defined")
	)
    )


