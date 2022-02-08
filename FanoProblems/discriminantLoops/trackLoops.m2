--file for tracking discriminant loops from each directory
--restart
loadPackage("Monodromy",FileName=>"../../../Monodromy.m2")

FanoProbs = select(readDirectory "./",s->(s != ".") and (s != "..") and (s != "trackLoops.m2") and (s != "trackLoops.m2~") and (s != "1_5_2_4"))

for FanoFolder in FanoProbs do (
    print("Loading monodromy from folder "|FanoFolder);
    M = loadMonodromy(FanoFolder|"/monod.txt");
    pts = value get (FanoFolder|"/pts.txt");
    if isWellDefined M then (
	print("Computing monodromy loop for folder "|FanoFolder);
	try(monodromyLoop(M,pts,Verbosity=>2)) else print("--Failed")
	) else (
	print("Monodromy object not well-defined")
	)
    )
