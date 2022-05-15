FanoGaloisGroups
================

Repository for software and computations related to Galois groups of finite Fano problems
-----------------------------------------------------------------------------------------

Repository created and maintained by:

> Thomas Yahl
> 
> Thomasjyahl@tamu.edu
> 
> Texas A&M University

The problem of enumerating linear spaces of a fixed dimension on a variety is known as a Fano problem. Those Fano problems with finitely many solutions have an associated Galois group that acts on the set of solutions. For a class of Fano problems, Hashimoto and Kadets determined the Galois group completely and showed that for all other Fano problems the Galois group contains the alternating group on its solutions. For Fano problems of moderate size with as yet undetermined Galois group, computational methods has been used to prove that the Galois group is the full symmetric group. For larger examples numerical methods provide evidence the Galois group is the full symmetric group. This repository contains code written in Macaulay2 (M2) and data for replicating these computations.

The repository is split into several folders:

1) The folder 'MEGA_Abstract' contains an article describing finite Fano problems, their Galois groups, and the computational methods used in other folders.
2) The folder 'Software' contains softwares necessary for running computations from other folders. The contents are as follows:
   - an executable version of the software 'alphaCertified'.
   - a M2 package 'alphaCertified' for transcribing systems and solutions to alphaCertified readable files.
   - a M2 package 'FanoProblems' for numerically tracking solutions of a parameterized system of equations around a designated path in the parameter space.
3) The folder 'Certification' contains several files:
   - several folders containing data for different finite Fano problems.
   - a M2 file 'systemGeneration.m2' for generating contents of the data contained in the mentioned folders.
   - a M2 file 'certifySolutions.m2' for running contents of the data folders through alphaCertified.
4) The folder 'Monodromy' contains several files:
   - several folders containing data for different finite Fano problems.
   - a M2 file 'discriminant_loops.m2' for generating contents of the data contained in previous folders.
   - a M2 file 'trackLoops.m2' for running contents of the data folders through alphaCertified.
5) The folder 'StartSystems' contains data used for more quickly computing contents in previous folders.

