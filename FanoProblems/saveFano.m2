restart
loadPackage("Monodromy",FileName=>"../../Monodromy.m2")
M = FanoMonodromy(1,4,1:5,Verbose=>true)
isWellDefined M
max apply(M#alphaConstants,last)
sum(M#alphaConstants,last)/(#M#baseSolutions)
max apply(M#alphaConstants,l->l#1)
sum(M#alphaConstants,l->l#1)/(#M#baseSolutions)
max apply(M#alphaConstants,first)
sum(M#alphaConstants,first)/(#M#baseSolutions)
monodromyLoop(M,10,Verbosity=>1)
printCycleTally M
saveMonodromy(M,"1_4_5.txt")

--refinement
M = refine(M,Bits=>60)
--for checking
N = loadMonodromy("1_4_5.txt")
