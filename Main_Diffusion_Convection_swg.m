% PURPOSE:
% Solve the Convection Diffusion equation by Lowest Order Weak Galerkin FEM 
%
% Convection Diffusion equation in 2 Dimension case:
%
%       -div(\alpha \nabla u) + \beta \cdot \nabla u + c u = f, 
%        u = g on \partial \Omega
%
% Run this function to get the numerical solutions. But the following 
% configurations are needed:
%
%    Part I: Enter the PDE information. This can be done in two ways 
%
%            I.1: Edit the file Diffusion_convection_MyPDE.m, and select 
%                 this file later in Main_Stokes_swg.m (this file).
%            I.2: Select one of the test examples later in this file.
%
%    Part II: Enter the domain and finite element partition options later
%             in Main_Diffusion_convection_swg.m (this file).
%
%  Moreover, the computation on Lshape domain is possible, only have to
%  replace line 121-122 by line 116-119
%
%  Copyright (C)  Yujie LIU.  Junping WANG. See COPYRIGHT.txt for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pde Model Problem structure 
%
%% Select your pde Option 1: try the code with own PDE and domain data  
%
% pde = Diffusion_convection_MyPDE;  
% domain = [0.0 2.0 0.0 2.0];   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PDE Option 2: Test example 1 as specified in Diffusion_convection_poly_test1.m
%     
% pde = Diffusion_convection_poly_test1;          
% domain = [0.0 1.0 0.0 1.0]; 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PDE Option 3: Test example 2 as specified in  Diffusion_convection_poly_test2.m
%
% pde = Diffusion_convection_poly_test2;          
% domain = [0.0 1.0 0.0 1.0]; 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%PDE Option 4: Test example 3 as specified in  Diffusion_convection_poly_test3.m
%
% pde = Diffusion_convection_sinsindata_test3;    
% domain = [0.0 1.0 0.0 1.0]; 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PDE Option 5: Test example 4 as specified in  Diffusion_convection_poly_test4.m
%
 pde = Diffusion_convection_sinsindata_test4;    
 domain = [0.0 2.0 0.0 2.0];
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Enter domain partition information
%
% Enter the level of grid refinement: Set maxIt =3 or larger values if you know 
% the exact sol and would like to see the errors and convergence.
maxIt = 4;
%
% Computation parameter: meshsize
h=[1/2 1/4 1/8 1/16 1/32 1/64 1/128];
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Mesh options
%
%   option.type  ='triangle';  % triangular elements
%
    option.type  ='quadrangle'; % rectangular elements
%
%  option.type  ='hexagon';
%      Must select suboption for hexamesh. 
%        option.hexamesh ='regular';
%        option.hexamesh ='irregular_one';
%        option.hexamesh ='irregular_two'; 
        option.hexamesh ='irregular_three';

% option.type   ='octagon';
%       Must select suboption for octamesh
%           option.octamesh ='regular';
           option.octamesh ='irregular_one';
%           option.octamesh ='irregular_two';
%           option.octamesh ='irregular_three';
%
% DONE PDE and Domain Configratuin. You can run this funtion NOW.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N     = zeros(maxIt,1);
H     = zeros(maxIt,1);
ERR_triangle   = zeros(maxIt,2);
ERR_hexagon    = zeros(maxIt,2);
ERR_quadrangle = zeros(maxIt,2);
ERR_octagon    = zeros(maxIt,2);
%
% Option solver
option.solver = 'direct';
option.exact =pde.KnownSol;
%
% Stabliser parameters
% this is the stabilization parameter in SWG, any value>\epsilon >0 
% is admissible
option.kappa=4.0; 
%  kappa * h^{alpha} is the factor of the stabilizer, usually we use
%  alpha= -1.0.
option.alpha=-1 ; 
%
%% Computation parameters
for k = 1:maxIt
    H(k)= h(k);
    
    meshtype =option.type;
    switch meshtype
        
    % triangle partition       
    case 'triangle'
        
%     Lshape  domain computation     
%     [node,elem] = squaremesh([-1,1,-1,1],h(k));
%     [node,elem] = delmesh(node,elem,'x>0 & y<0');
%     bdFlag = setboundary(node,elem,'Dirichlet');
     
     [node,elem] = squaremesh(domain, h(k));
     bdFlag = setboundary(node,elem,'Dirichlet');
     [u,AD,b,err_tria,info] = ...
     Diffusion_Convection_SWG_triangle(node,elem,pde,bdFlag,option);   
     ERR_triangle(k,1)=err_tria.L2;
     ERR_triangle(k,2)=err_tria.H1;     
     %
     figure(1);
     showmesh(node,elem);
     hold off;
     %
     fprintf('The computation for meshsize h = 1/%d is finished. \n', 1/h(k)); 
     %
     if (option.exact==1)
     if (k==maxIt)
         if (maxIt > 2)
     figure('Name',' L^2 Error '); 
     showrate(H, ERR_triangle(:,1),1,'-*','L2 Error ');
     figure('Name',' H^1 Error '); 
     showrate(H, ERR_triangle(:,2),1,'-*','H^1 Error ');
     hold off;
         end
     %
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     if ~exist('Saved_Results','dir')
         mkdir Saved_Results;    
     end
     fileID = fopen('.\Saved_Results\Triangle_Convection_Diffusion.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s \n','Mid_Edge_x','Mid_Edge_y','Num u     ', 'Num_u - Exact_u');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Mid_Edge(:,1) Mid_Edge(:,2) u info.error]');
     fclose(fileID);
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Triangle_Convection_Diffusion.txt in the folder Saved_Results.\n', 1/h(maxIt))
     %
     rateL2 =zeros(maxIt,1);
     rateH1 =zeros(maxIt,1);
     for kk = 2:maxIt
     rateL2(kk)     = (log (ERR_triangle(kk,1))-log (ERR_triangle(kk-1,1)))/ (log (h(kk))-log (h(kk-1))) ;   
     rateH1(kk)     = (log (ERR_triangle(kk,2))-log (ERR_triangle(kk-1,2)))/ (log (h(kk))-log (h(kk-1))) ;  
     end
     A = [ERR_triangle(:,1) rateL2 ERR_triangle(:,2) rateH1]';
     fileID = fopen('.\Saved_Results\Triangle_Convection_Diffusion_error.txt','wt');
     fprintf(fileID,'%12s %6s %12s %6s \n','errL2   ','rateL2','errH1   ','rateH1');
     fprintf(fileID,'%12.2e %6.2f %12.2e %6.2f\n',A);
     fclose(fileID);
     fprintf('The error and convergence rate are saved in file Triangle_Convection_Diffusion_error.txt in folder Saved_Results.\n');
     end
     elseif (option.exact==0)
     if (k==maxIt)   
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     if ~exist('Saved_Results','dir')
         mkdir Saved_Results;    
     end
     fileID = fopen('.\Saved_Results\Triangle_Convection_Diffusion.txt','wt');
     fprintf(fileID,'%12s %12s %12s  \n','Mid_Edge_x','Mid_Edge_y','Num u     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f \n',[Mid_Edge(:,1) Mid_Edge(:,2) u ]');
     fclose(fileID);
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Triangle_Convection_Diffusion.txt in the folder Saved_Results.\n', 1/h(maxIt)) 
     end
     end
     
    % hexagon partition
    case 'hexagon'
    
     [node,elem] = squaremesh(domain, h(k));
     bdFlag = setboundary(node,elem,'Dirichlet');
     [Nodeh,Elemh,u,AD,b,err_hexa,info] = ...
     Diffusion_Convection_SWG_hexagon(node,elem,pde,bdFlag,option);
     ERR_hexagon(k,1)=err_hexa.L2;
     ERR_hexagon(k,2)=err_hexa.H1; 
     %
     fprintf('The computation for meshsize h = 1/%d is finished. \n', 1/h(k)); 
     %
     if (option.exact==1)
     if (k==maxIt)
         if (maxIt > 2)
     figure('Name',' L^2 Error '); 
     showrate(H, ERR_hexagon(:,1),1,'-*','L^2 Error ');
     figure('Name',' H^1 Error '); 
     showrate(H, ERR_hexagon(:,2),1,'-*','H^1 Error ');
     hold off;
         end
          %
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     if ~exist('Saved_Results','dir')
         mkdir Saved_Results;    
     end
     fileID = fopen('.\Saved_Results\Hexagon_Convection_Diffusion.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s \n','Mid_Edge_x','Mid_Edge_y','Num u     ', 'Num_u - Exact_u');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Mid_Edge(:,1) Mid_Edge(:,2) u info.error]');
     fclose(fileID);
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Hexagon_Convection_Diffusion.txt in folder Saved_Results.\n', 1/h(maxIt))
     %
     rateL2 =zeros(maxIt,1);
     rateH1 =zeros(maxIt,1);
     for kk = 2:maxIt
     rateL2(kk)     = (log (ERR_hexagon(kk,1))-log (ERR_hexagon(kk-1,1)))/ (log (h(kk))-log (h(kk-1))) ;   
     rateH1(kk)     = (log (ERR_hexagon(kk,2))-log (ERR_hexagon(kk-1,2)))/ (log (h(kk))-log (h(kk-1))) ;  
     end
     A = [ERR_hexagon(:,1) rateL2 ERR_hexagon(:,2) rateH1]';
     fileID = fopen('.\Saved_Results\Hexagon_Convection_Diffusion_error.txt','wt');
     fprintf(fileID,'%12s %6s %12s %6s \n','errL2   ','rateL2','errH1   ','rateH1');
     fprintf(fileID,'%12.2e %6.2f %12.2e %6.2f\n',A);
     fclose(fileID);
     fprintf('The error and convergence rate are saved in file Hexagon_Convection_Diffusion_error.txt in folder Saved_Results.\n');
     end
     elseif (option.exact==0)
     if (k==maxIt)  
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     if ~exist('Saved_Results','dir')
         mkdir Saved_Results;    
     end
     fileID = fopen('.\Saved_Results\Hexagon_Convection_Diffusion.txt','wt');
     fprintf(fileID,'%12s %12s %12s \n','Mid_Edge_x','Mid_Edge_y','Num u     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f  \n',[Mid_Edge(:,1) Mid_Edge(:,2) u ]');
     fclose(fileID);
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Hexagon_Convection_Diffusion.txt in folder Saved_Results.\n', 1/h(maxIt))  
     end
     end
     
    % quadrangle partition
    case 'quadrangle'
 
    [nodeq,elemq] = squarequadmesh(domain,h(k));
    bdFlag = setboundary(nodeq,elemq,'Dirichlet');
    
   % Lshape domain computation  
   % [nodeq,elemq] = squarequadmesh([-1,1,-1,1],h(k));
   % [nodeq,elemq] = delmesh(nodeq,elemq,'x>0 & y<0');
   % bdFlag = setboundary(nodeq,elemq,'Dirichlet');
     
    [u,AD,b,err_quad,info] = ...
    Diffusion_Convection_SWG_quadrangle(nodeq,elemq,pde,bdFlag,option);
    ERR_quadrangle(k,1)=err_quad.L2;
    ERR_quadrangle(k,2)=err_quad.H1;     
     figure(2);
     showmesh(nodeq,elemq);
     hold off;
     %
     fprintf('The computation for meshsize h = 1/%d is finished. \n', 1/h(k)); 
     %
     if (option.exact==1)
     if (k==maxIt)
         if (maxIt > 2)
     figure('Name',' L^2 Error '); 
     showrate(H, ERR_quadrangle(:,1),1,'-*','L^2 Error ');
     figure('Name',' H^1 Error '); 
     showrate(H, ERR_quadrangle(:,2),1,'-*','H^1 Error ');
     hold off;
         end
     %
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     if ~exist('Saved_Results','dir')
         mkdir Saved_Results;    
     end
     fileID = fopen('.\Saved_Results\Quadrangle_Convection_Diffusion.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s \n','Mid_Edge_x','Mid_Edge_y','Num u     ', 'Num_u - Exact_u');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Mid_Edge(:,1) Mid_Edge(:,2) u info.error]');
     fclose(fileID);
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Quadrangle_Convection_Diffusion.txt in folder Saved_Results.\n', 1/h(maxIt))
     %
     rateL2 =zeros(maxIt,1);
     rateH1 =zeros(maxIt,1);
     for kk = 2:maxIt
     rateL2(kk)     = (log (ERR_quadrangle(kk,1))-log (ERR_quadrangle(kk-1,1)))/ (log (h(kk))-log (h(kk-1))) ;   
     rateH1(kk)     = (log (ERR_quadrangle(kk,2))-log (ERR_quadrangle(kk-1,2)))/ (log (h(kk))-log (h(kk-1))) ;  
     end
     A = [ERR_quadrangle(:,1) rateL2 ERR_quadrangle(:,2) rateH1]';
     fileID = fopen('.\Saved_Results\Quadrangle_Convection_Diffusion_error.txt','wt');
     fprintf(fileID,'%12s %6s %12s %6s \n','errL2   ','rateL2','errH1   ','rateH1');
     fprintf(fileID,'%12.2e %6.2f %12.2e %6.2f\n',A);
     fclose(fileID);
     fprintf('The error and convergence rate are saved in file Quadrangle_Convection_Diffusion_error.txt in folder Saved_Results.\n');
     end
     elseif (option.exact==0)
     if (k==maxIt)  
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     if ~exist('Saved_Results','dir')
         mkdir Saved_Results;    
     end
     fileID = fopen('.\Saved_Results\Quadrangle_Convection_Diffusion.txt','wt');
     fprintf(fileID,'%12s %12s %12s \n','Mid_Edge_x','Mid_Edge_y','Num u     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f \n',[Mid_Edge(:,1) Mid_Edge(:,2) u]');
     fclose(fileID);
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Quadrangle_Convection_Diffusion.txt in folder Saved_Results.\n', 1/h(maxIt))
     end
     end
    % octagon partition
    case 'octagon'
     [nodeq,elemq] = squarequadmesh(domain,h(k));
     bdFlag = setboundary(nodeq,elemq,'Dirichlet');
     [u,AD,b,err_octa,info] = ...
     Diffusion_Convection_SWG_octagon(nodeq,elemq,pde,bdFlag,option);
     ERR_octagon(k,1)=err_octa.L2;
     ERR_octagon(k,2)=err_octa.H1;     
     %
     fprintf('The computation for meshsize h = 1/%d is finished. \n', 1/h(k)); 
     %
     if (option.exact==1)
     if (k==maxIt)
         if (maxIt > 2)
     figure('Name',' L^2 Error '); 
     showrate(H, ERR_octagon(:,1),1,'-*','L^2 Error ');
     figure('Name',' H^1 Error '); 
     showrate(H, ERR_octagon(:,2),1,'-*','H^1 Error ');
     hold off;
         end
     %
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     if ~exist('Saved_Results','dir')
         mkdir Saved_Results;    
     end
     fileID = fopen('.\Saved_Results\Octagon_Convection_Diffusion.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s \n','Mid_Edge_x','Mid_Edge_y','Num u     ', 'Num_u - Exact_u');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Mid_Edge(:,1) Mid_Edge(:,2) u info.error]');
     fclose(fileID);
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Octagon_Convection_Diffusion.txt in folder Saved_Results.\n', 1/h(maxIt))
     %
     rateL2 =zeros(maxIt,1);
     rateH1 =zeros(maxIt,1);
     for kk = 2:maxIt
     rateL2(kk)     = (log (ERR_octagon(kk,1))-log (ERR_octagon(kk-1,1)))/ (log (h(kk))-log (h(kk-1))) ;   
     rateH1(kk)     = (log (ERR_octagon(kk,2))-log (ERR_octagon(kk-1,2)))/ (log (h(kk))-log (h(kk-1))) ;  
     end
     A = [ERR_octagon(:,1) rateL2 ERR_octagon(:,2) rateH1]';
     if ~exist('Saved_Results','dir')
         mkdir Saved_Results;    
     end
     fileID = fopen('.\Saved_Results\Octagon_Convection_Diffusion_error.txt','wt');
     fprintf(fileID,'%12s %6s %12s %6s \n','errL2   ','rateL2','errH1   ','rateH1');
     fprintf(fileID,'%12.2e %6.2f %12.2e %6.2f\n',A);
     fclose(fileID);
     fprintf('The error and convergence rate are saved in file Octagon_Convection_Diffusion_error.txt in folder Saved_Results.\n');
     end
     elseif (option.exact==0)
     if (k==maxIt)  
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     if ~exist('Saved_Results','dir')
         mkdir Saved_Results;    
     end
     fileID = fopen('.\Saved_Results\Octagon_Convection_Diffusion.txt','wt');
     fprintf(fileID,'%12s %12s %12s \n','Mid_Edge_x','Mid_Edge_y','Num u     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f \n',[Mid_Edge(:,1) Mid_Edge(:,2) u ]');
     fclose(fileID);
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Octagon_Convection_Diffusion.txt in folder Saved_Results.\n', 1/h(maxIt))
     end
     end
    end
end
       

    
