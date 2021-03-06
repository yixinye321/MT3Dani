# MT3Dani
1. This code is an finite element modeling program for solving the 3D MT problem in anisotropic media using vector-scalar potentials with unstructured mesh.
 
2. The program is written in language Fortran 90, and performed in Linux with Intel Fortran compiler.

3. Three files are needed to execute the program, one is the "RUNFILE", which inputs the parameters for MT responses computations, another one is the ".node" file, which consists of coordinates of all vertices in the mesh, and the third one is the element file “.ele”  that contains all the elements in the mesh. The ".node" and “.ele” are the mesh files generated by the tetrahedral mesh generator TetGen.

4. The ".mt" file is output results of the MT responses calculated using the program.

5. The input file ".poly" represents a piecewise linear complex (PLC) as well as some additional information for generating mesh with TetGen. It consists of four parts, which are a list of points, a list of facets, a list of (volume) hole points, and a list of region attributes, respectively. The first three parts are mandatory, but the fourth part is optional. 

6. The ".vol" file is the volume constraint file for local mesh refinement in the nex step, it was created using the matlab program "lvr.m" according to the mesh refinement requirement.
