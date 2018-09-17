% PURPOSE:
% Solve the Stokes Equation by Lowest Order Weak Galerkin FEM 
% 
% Stokes equation in 2-Dimension: 
%     - div (nu \nabla (u,v)) + \nabla p = (f1, f2),  in \Omega 
%                 div (u,v) = f3,     in \Omega
%                       u = g on \partial \Omega
%
% Run this function to get the numerical solutions. But the following 
% configurations are needed:
%
%    Part I: Enter the PDE information. This can be done in two ways 
%
%            I.1: Edit the file Stokes_MyPDE.m, and select this file later 
%                 in Main_Stokes_swg.m (this file)
%            I.2: Select one of the test examples later in this file.
%
%    Part II: Enter the domain and finite element partition options later
%             in Main_Stokes.swg.m (this file).
%
%
%  Copyright (C) Yujie LIU. Junping WANG. See COPYRIGHT.txt for details.
%  2018.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Select your pde Option 1: try the code with own PDE and domain data  

% pde = Stokes_MyPDE;  
% domain = [0.0 1.5 0.0 1.5];   

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PDE Option 2: Test example 1 as specified in Stokes_sincos.m
%
 pde = Stokes_sincos;  
 domain = [0.0 2 0.0 2];  
%
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PDE Option 3: Test example 2 as specified in Stokes_test.m
%
% pde = Stokes_test;  
% domain = [0.0 1 0.0 1];
%
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PDE Option 4: Test example 3 as specified in Stokes_polynom.m
%
% pde = Stokes_polynom;  
% domain = [0.0 1 0.0 1];
%
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PDE Option 5: Test example 4 as specified in Stokes_viscosity_function.m
%
% pde = Stokes_viscosity_function;  
% domain = [0.0 1 0.0 1]; 
%
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PDE Option 6: Test example 5 as specified in Stokes_lid_driven.m
%               (exact solution is not known)
%            
% pde = Stokes_lid_driven_cavity;
% domain = [0.0 1.0 0.0 1.0];
%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Enter domain partition information
%
% Enter the level of grid refinement: Set maxIt = 3 or larger values if you know 
% the exact sol and would like to see the errors and convergence.
maxIt = 3;
%
% Computation parameter: meshsize
h=[1/4 1/8 1/16 1/32];
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Mesh options
%
% option.type  ='triangle';  % triangular elements
%
 option.type  ='quadrangle'; % rectangular elements
%
% option.type  ='hexagon';
%      Must select one of the suboptions for hexamesh. 
%        option.hexamesh ='regular';
        option.hexamesh ='irregular_one';
%        option.hexamesh ='irregular_two'; 
%        option.hexamesh ='irregular_three';

%  option.type   ='octagon';
%       Must select one of the suboptions for octamesh
%           option.octamesh ='regular';
           option.octamesh ='irregular_one';
%           option.octamesh ='irregular_two';
%           option.octamesh ='irregular_three';
%
% DONE PDE and Domain Configratuin. You can run this funtion NOW.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
N     = zeros(maxIt,1);
H     = zeros(maxIt,1);
ERR_triangle   = zeros(maxIt,3);
ERR_hexagon    = zeros(maxIt,3);
ERR_quadrangle = zeros(maxIt,3);
ERR_octagon    = zeros(maxIt,3);
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
%% Computation
for k = 1:maxIt
    H(k)= h(k);
    
    meshtype =option.type;
    switch meshtype
        
    % triangle partition        
    case 'triangle'
     [node,elem] = squaremesh(domain, h(k));
     bdFlag = setboundary(node,elem,'Dirichlet');
     [u,AD,b,err_tria,info] = Stokes_SWG_triangle(node,elem,domain,pde,bdFlag,option);   
     ERR_triangle(k,1)=err_tria.L2;
     ERR_triangle(k,2)=err_tria.H1;     
     ERR_triangle(k,3)=err_tria.p;

     figure(1);
     showmesh(node,elem);
     hold off;
     
     fprintf('The computation for meshsize h = 1/%d is finished. \n', 1/h(k));
     
     if (option.exact==1)
     if (k==maxIt)
         if (maxIt > 2)
     figure('Name',' L^2 for velocity (u,v)_h '); 
     showrate(H, ERR_triangle(:,1),1,'-*','L^2 for velocity (u,v)_h ');
     figure('Name',' H^1 for velocity (u,v)_h '); 
     showrate(H, ERR_triangle(:,2),1,'-*','H^1 for velocity (u,v)_h ');
     figure('Name',' L^2 for pressure p_h '); 
     showrate(H, ERR_triangle(:,3),1,'-*','L^2 for pressure p_h ');
     hold off;
         end
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     NE         = info.NE;
     fileID = fopen('.\Save_Results\Triangle_Stokes.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s  %12s %12s\n','Mid_Edge_x','Mid_Edge_y','Num u     ', 'Num_u - Exact_u','Num v     ', 'Num_v - Exact_v');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n',[Mid_Edge(:,1) Mid_Edge(:,2) u(1:NE) info.erroru u(NE+1:2*NE) info.errorv]' );
     fprintf(fileID,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
     fprintf(fileID,'%12s %12s %12s %12s \n','Elem_center_x','Elem_center_y','Num p     ', 'Num_p - Exact_p');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Elem_center(:,1) Elem_center(:,2) u(2*NE+1:end) info.errorp]' );
     fclose(fileID);
     %
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Triangle_Stokes.txt in folder Save_Results.\n', 1/h(maxIt))
     %
     rateL2 =zeros(maxIt,1);
     rateH1 =zeros(maxIt,1);
     rateL2p =zeros(maxIt,1);
     for kk = 2:maxIt
     rateL2(kk)     = (log (ERR_triangle(kk,1))-log (ERR_triangle(kk-1,1)))/ (log (h(kk))-log (h(kk-1))) ;   
     rateH1(kk)     = (log (ERR_triangle(kk,2))-log (ERR_triangle(kk-1,2)))/ (log (h(kk))-log (h(kk-1))) ;  
     rateL2p(kk)    = (log (ERR_triangle(kk,3))-log (ERR_triangle(kk-1,3)))/ (log (h(kk))-log (h(kk-1))) ; 
     end
     A = [ERR_triangle(:,1) rateL2 ERR_triangle(:,2) rateH1 ERR_triangle(:,3) rateL2p]';
     fileID = fopen('.\Save_Results\Triangle_Stokes_error.txt','wt');
     fprintf(fileID,'%12s %6s %12s %6s %12s %6s \n','errL2   ','rateL2','errH1   ','rateH1','errL2 p ','rateL2 p');
     fprintf(fileID,'%12.2e %6.2f %12.2e %6.2f %12.2e %6.2f\n',A);
     fclose(fileID);
     fprintf('The error and convergence rate are saved in file Triangle_Stokes_error.txt in folder Save_Results.\n');
     end
     elseif (option.exact==0)
     if (k==maxIt)
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     NE         = info.NE;
     fileID = fopen('.\Save_Results\Triangle_Stokes.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s \n','Mid_Edge_x','Mid_Edge_y','Num u     ','Num v     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Mid_Edge(:,1) Mid_Edge(:,2) u(1:NE) u(NE+1:2*NE)]' );
     fprintf(fileID,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
     fprintf(fileID,'%12s %12s %12s \n','Elem_center_x','Elem_center_y','Num p     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f \n',[Elem_center(:,1) Elem_center(:,2) u(2*NE+1:end) ]' );
     fclose(fileID);
     %
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Triangle_Stokes.txt in folder Save_Results.\n', 1/h(maxIt))
     % 
     end
     end;
    %hexagon partition     
    case 'hexagon'
     %domain = [0.0 1 0.0 1]; 
     [node,elem] = squaremesh(domain, h(k));
     bdFlag = setboundary(node,elem,'Dirichlet');
     [Nodeh,Elemh,u,AD,b,err_hexa,info] = Stokes_SWG_hexagon(node,elem,domain,pde,bdFlag,option);
     ERR_hexagon(k,1)=err_hexa.L2;
     ERR_hexagon(k,2)=err_hexa.H1;     
     ERR_hexagon(k,3)=err_hexa.p;
     fprintf('The computation for meshsize h = 1/%d is finished. \n', 1/h(k));
     %
     if (option.exact==1)
     if (k==maxIt)
         if (maxIt > 2)
     figure('Name',' L^2 for velocity (u,v)_h '); 
     showrate(H, ERR_hexagon(:,1),1,'-*','L^2 for velocity (u,v)_h ');
     figure('Name',' H^1 for velocity (u,v)_h '); 
     showrate(H, ERR_hexagon(:,2),1,'-*','H^1 for velocity (u,v)_h ');
     figure('Name',' L^2 for pressure p_h '); 
     showrate(H, ERR_hexagon(:,3),1,'-*','L^2 for pressure p_h ');
     hold off;
         end
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     NE         = info.NE;
     fileID = fopen('.\Save_Results\Hexagon_Stokes.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s  %12s %12s\n','Mid_Edge_x','Mid_Edge_y','Num u     ', 'Num_u - Exact_u','Num v     ', 'Num_v - Exact_v');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n',[Mid_Edge(:,1) Mid_Edge(:,2) u(1:2*NE) info.erroru u(2*NE+1:4*NE) info.errorv]' );
     fprintf(fileID,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
     fprintf(fileID,'%12s %12s %12s %12s \n','Elem_center_x','Elem_center_y','Num p     ', 'Num_p - Exact_p');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Elem_center(:,1) Elem_center(:,2) u(4*NE+1:end) info.errorp]' );
     fclose(fileID);
     %
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Hexagon_Stokes.txt in folder Save_Results.\n', 1/h(maxIt))
     %
     rateL2 =zeros(maxIt,1);
     rateH1 =zeros(maxIt,1);
     rateL2p =zeros(maxIt,1);
     for kk = 2:maxIt
     rateL2(kk)     = (log (ERR_hexagon(kk,1))-log (ERR_hexagon(kk-1,1)))/ (log (h(kk))-log (h(kk-1))) ;   
     rateH1(kk)     = (log (ERR_hexagon(kk,2))-log (ERR_hexagon(kk-1,2)))/ (log (h(kk))-log (h(kk-1))) ;  
     rateL2p(kk)    = (log (ERR_hexagon(kk,3))-log (ERR_hexagon(kk-1,3)))/ (log (h(kk))-log (h(kk-1))) ; 
     end
     A = [ERR_hexagon(:,1) rateL2 ERR_hexagon(:,2) rateH1 ERR_hexagon(:,3) rateL2p]';
     fileID = fopen('.\Save_Results\Hexagon_Stokes_error.txt','wt');
     fprintf(fileID,'%12s %6s %12s %6s %12s %6s \n','errL2   ','rateL2','errH1   ','rateH1','errL2 p ','rateL2 p');
     fprintf(fileID,'%12.2e %6.2f %12.2e %6.2f %12.2e %6.2f\n',A);
     fclose(fileID);
     fprintf('The error and convergence rate are saved in file Hexagon_Stokes_error.txt in folder Save_Results.\n');
     end
     elseif (option.exact==0)
     if (k==maxIt)
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     NE         = info.NE;
     fileID = fopen('.\Save_Results\Hexagon_Stokes.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s \n','Mid_Edge_x','Mid_Edge_y','Num u     ','Num v     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Mid_Edge(:,1) Mid_Edge(:,2) u(1:2*NE) u(2*NE+1:4*NE)]' );
     fprintf(fileID,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
     fprintf(fileID,'%12s %12s %12s \n','Elem_center_x','Elem_center_y','Num p     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f \n',[Elem_center(:,1) Elem_center(:,2) u(4*NE+1:end) ]' );
     fclose(fileID);
     %
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Hexagon_Stokes.txt in folder Save_Results.\n', 1/h(maxIt))
     % 
     end
     end;
    % quadrangle partition     
    case 'quadrangle'
    %domain = [0,pi,0,pi];
    %domain = [0,1,0,1];
    [nodeq,elemq] = squarequadmesh(domain,h(k));
    bdFlag = setboundary(nodeq,elemq,'Dirichlet');
    [u,AD,b,err_quad,info] = Stokes_SWG_quadrangle(nodeq,elemq,domain,pde,bdFlag,option);
     ERR_quadrangle(k,1)=err_quad.L2;
     ERR_quadrangle(k,2)=err_quad.H1;     
     ERR_quadrangle(k,3)=err_quad.p;
     
     figure(2);
     showmesh(nodeq,elemq);
     hold off;
     
     fprintf('The computation for meshsize h = 1/%d is finished. \n', 1/h(k));
     
     if (option.exact==1)
     if (k==maxIt)
         if (maxIt > 2)
     figure('Name',' L^2 for velocity (u,v)_h '); 
     showrate(H, ERR_quadrangle(:,1),1,'-*','L^2 for velocity (u,v)_h ');
     figure('Name',' H^1 for velocity (u,v)_h '); 
     showrate(H, ERR_quadrangle(:,2),1,'-*','H^1 for velocity (u,v)_h ');
     figure('Name',' L^2 for pressure p_h '); 
     showrate(H, ERR_quadrangle(:,3),1,'-*','L^2 for pressure p_h ');
     hold off;
         end
     %
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     NE         = info.NE;
     fileID = fopen('.\Save_Results\Quadrangle_Stokes.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s  %12s %12s\n','Mid_Edge_x','Mid_Edge_y','Num u     ', 'Num_u - Exact_u','Num v     ', 'Num_v - Exact_v');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n',[Mid_Edge(:,1) Mid_Edge(:,2) u(1:NE) info.erroru u(NE+1:2*NE) info.errorv]' );
     fprintf(fileID,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
     fprintf(fileID,'%12s %12s %12s %12s \n','Elem_center_x','Elem_center_y','Num p     ', 'Num_p - Exact_p');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Elem_center(:,1) Elem_center(:,2) u(2*NE+1:end) info.errorp]' );
     fclose(fileID);
     %
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Quadrangle_Stokes.txt in folder Save_Results.\n', 1/h(maxIt))
     %
     rateL2 =zeros(maxIt,1);
     rateH1 =zeros(maxIt,1);
     rateL2p =zeros(maxIt,1);
     for kk = 2:maxIt
     rateL2(kk)     = (log (ERR_quadrangle(kk,1))-log (ERR_quadrangle(kk-1,1)))/ (log (h(kk))-log (h(kk-1))) ;   
     rateH1(kk)     = (log (ERR_quadrangle(kk,2))-log (ERR_quadrangle(kk-1,2)))/ (log (h(kk))-log (h(kk-1))) ;  
     rateL2p(kk)    = (log (ERR_quadrangle(kk,3))-log (ERR_quadrangle(kk-1,3)))/ (log (h(kk))-log (h(kk-1))) ; 
     end
     A = [ERR_quadrangle(:,1) rateL2 ERR_quadrangle(:,2) rateH1 ERR_quadrangle(:,3) rateL2p]';
     fileID = fopen('.\Save_Results\Quadrangle_Stokes_error.txt','wt');
     fprintf(fileID,'%12s %6s %12s %6s %12s %6s \n','errL2   ','rateL2','errH1   ','rateH1','errL2 p ','rateL2 p');
     fprintf(fileID,'%12.2e %6.2f %12.2e %6.2f %12.2e %6.2f\n',A);
     fclose(fileID);
     fprintf('The error and convergence rate are saved in file Quadrangle_Stokes_error.txt in folder Save_Results.\n');
     end
     elseif (option.exact==0)
     if (k==maxIt)
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     NE         = info.NE;
     fileID = fopen('.\Save_Results\Quadrangle_Stokes.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s \n','Mid_Edge_x','Mid_Edge_y','Num u     ','Num v     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Mid_Edge(:,1) Mid_Edge(:,2) u(1:NE) u(NE+1:2*NE)]' );
     fprintf(fileID,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
     fprintf(fileID,'%12s %12s %12s \n','Elem_center_x','Elem_center_y','Num p     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f \n',[Elem_center(:,1) Elem_center(:,2) u(2*NE+1:end) ]' );
     fclose(fileID);
     %
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Quadrangle_Stokes.txt in folder Save_Results.\n', 1/h(maxIt))
     % 
     end
     end;
    % octagon partition     
    case 'octagon'
     %domain = [0,1,0,1];
     [nodeq,elemq] = squarequadmesh(domain,h(k));
     bdFlag = setboundary(nodeq,elemq,'Dirichlet');
     [u,AD,b,err_octa,info] = Stokes_SWG_octagon(nodeq,elemq,domain,pde,bdFlag,option);
     ERR_octagon(k,1)=err_octa.L2;
     ERR_octagon(k,2)=err_octa.H1;     
     ERR_octagon(k,3)=err_octa.p;
     fprintf('The computation for meshsize h = 1/%d is finished. \n', 1/h(k));
     
     if (option.exact==1)
     if (k==maxIt)
         if (maxIt > 2)
     figure('Name',' L^2 for velocity (u,v)_h '); 
     showrate(H, ERR_octagon(:,1),1,'-*','L^2 for velocity of (u,v)_h ');
     figure('Name',' H^1 for velocity (u,v)_h '); 
     showrate(H, ERR_octagon(:,2),1,'-*','H^1 for velocity (u,v)_h ');
     figure('Name',' L^2 for pressure p_h '); 
     showrate(H, ERR_octagon(:,3),1,'-*','L^2 for pressure p_h ');
     hold off;
         end
     %
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     NE         = info.NE;
     fileID = fopen('.\Save_Results\Octagon_Stokes.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s  %12s %12s\n','Mid_Edge_x','Mid_Edge_y','Num u     ', 'Num_u - Exact_u','Num v     ', 'Num_v - Exact_v');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n',[Mid_Edge(:,1) Mid_Edge(:,2) u(1:2*NE) info.erroru u(2*NE+1:4*NE) info.errorv]' );
     fprintf(fileID,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
     fprintf(fileID,'%12s %12s %12s %12s \n','Elem_center_x','Elem_center_y','Num p     ', 'Num_p - Exact_p');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Elem_center(:,1) Elem_center(:,2) u(4*NE+1:end) info.errorp]' );
     fclose(fileID);
     %
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Octagon_Stokes.txt in folder Save_Results.\n', 1/h(maxIt))
     %
     rateL2 =zeros(maxIt,1);
     rateH1 =zeros(maxIt,1);
     rateL2p =zeros(maxIt,1);
     for kk = 2:maxIt
     rateL2(kk)     = (log (ERR_octagon(kk,1))-log (ERR_octagon(kk-1,1)))/ (log (h(kk))-log (h(kk-1))) ;   
     rateH1(kk)     = (log (ERR_octagon(kk,2))-log (ERR_octagon(kk-1,2)))/ (log (h(kk))-log (h(kk-1))) ;  
     rateL2p(kk)    = (log (ERR_octagon(kk,3))-log (ERR_octagon(kk-1,3)))/ (log (h(kk))-log (h(kk-1))) ; 
     end
     A = [ERR_octagon(:,1) rateL2 ERR_octagon(:,2) rateH1 ERR_octagon(:,3) rateL2p]';
     fileID = fopen('.\Save_Results\Octagon_Stokes_error.txt','wt');
     fprintf(fileID,'%12s %6s %12s %6s %12s %6s \n','errL2   ','rateL2','errH1   ','rateH1','errL2 p ','rateL2 p');
     fprintf(fileID,'%12.2e %6.2f %12.2e %6.2f %12.2e %6.2f\n',A);
     fclose(fileID);
     fprintf('The error and convergence rate are saved in file Octagon_Stokes_error.txt in folder Save_Results.\n');
     end  
     elseif (option.exact==0)
     if (k==maxIt)
     Mid_Edge   = info.Mid_Edge;
     Elem_center= info.Elem_center;
     NE         = info.NE;
     fileID = fopen('.\Save_Results\Octagon_Stokes.txt','wt');
     fprintf(fileID,'%12s %12s %12s %12s \n','Mid_Edge_x','Mid_Edge_y','Num u     ','Num v     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f \n',[Mid_Edge(:,1) Mid_Edge(:,2) u(1:2*NE) u(2*NE+1:4*NE)]' );
     fprintf(fileID,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
     fprintf(fileID,'%12s %12s %12s \n','Elem_center_x','Elem_center_y','Num p     ');
     fprintf(fileID,'%12.4f %12.4f %12.4f \n',[Elem_center(:,1) Elem_center(:,2) u(4*NE+1:end) ]' );
     fclose(fileID);
     %
     fprintf('The numerical results for meshsize h = 1/%d is saved in file Octagon_Stokes.txt in folder Save_Results.\n', 1/h(maxIt))
     % 
     end
     end;
    end
end
       

    
