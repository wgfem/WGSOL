# WGSOL
WG MatLab functions for PDE solving

WGSOL a collection of MATLAB functions which implement the weak Galerkin (WG) finite element method in a simplified formulation (known as SWG – Simplified Weak Galerkin) for numerical solving of PDEs in two dimensions. The current version contains SWG solvers for the Poisson equation, the convection-diffusion-reaction equation, and the Stokes equation.

Steps to follow for trying out the code:

To use the sofware functions:

1. Download the released version, unzip the file, and save the source code into your favorite folder, say ./WGSOL

2. Set path configuration: go to the underlying directory (e.g., ./WGSOL) 
   and run the file PathConfigure.m.

3. In each subdirectory, .\Stokes, .\Poisson, .\Convection-diffusion,
    run the "Main_*_swg.m" file to test each solver. Read the header of the 
    main file "Main*_swg.m" for setting new configurations (such as PDE info, domain, 
    and partition options). In the main file, one can also set different 
    values for the stablizer parameter of the WG and element types.

Details regarding SWG can be found in the following articles:

1)  Simplified Weak Galerkin and New Finite Difference Schemes for the Stokes Equation, 
    2018, arXiv:1803.00120

2)  A Simplified Weak Galerkin Finite Element Method: Algorithm and Error Estimates,
    2018, arXiv:1808.08667
  


By Yujie Liu and Junping Wang
September 2018
