function pde = Stokes_lid_driven_cavity
%% This is the lid driven cavity problem for the stokes equation,
%  there is no exact solution for this problem, the boundary condtion for 
%  this problem is:
%      u1=1, 
%      u2=0, on y=1, x in (0,1)
% and 
%      u1=0,
%      u2=0, otherwise,
%
% the right hand term function is f1=0, f2=0.
% to simplify the computation other terms is set to zero.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%'f',        @f : load data (right hand side function)
%'exactu',   @exactu : exact solution of u
%'g_D',      @g_D :  Dirichlet boundary condition
%'Du',       @Du : Derivative of the exact solution u
%'d1',       @Diff1 : This is the viscosity of the fluid.
%'Df',       @Df: Derivative of the right hand side function
%'pp',       @pp: Exact soltion of p
%
%  Copyright (C)  Yujie LIU. Junping WANG. See COPYRIGHT.txt for details.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%


    % load data (right hand side function)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  f(p)
    % Get the point
    x = p(:,1); y = p(:,2);
    % define the first component of f in the momentum equation
    rhs(:,1) =0.*x;
    % define the second component of f in the momentum equation
    rhs(:,2) =0.*y;
    % define the RHS function in the mass conservation equation. Should be
    % zero for incompressible fluid. But the code takes any function. Make
    % sure that Dirichlet boundary data is compatible with this function
    rhs(:,3) =0.*x;
    end
    % Derivative of the right hand side function
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1,1) =  0.*x;
    rhs(:,1,2) =  0.*x;
    rhs(:,2,1) =  0.*y;
    rhs(:,2,2) =  0.*y;
    rhs(:,3,1) =  0.*x;
    rhs(:,3,2) =  0.*x;
    end
    % exact solution of u
    KnownSol=0; % set value = 1 if knowing the exact solution
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u(:,1) = 0.*x;
    u(:,2) = 0.*y;
    end
    % Dirichlet boundary condition
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function u =  g_D(p)
    x = p(:,1); y = p(:,2);
    testy1=(y>=0.9999);
    testxp=(x>=0.0001);
    testxm=(x<=0.9999);
    testbdvalue=(testy1+testxp+testxm)/3;
    testbd=(testbdvalue>=1 );
    u( testbd,1) =   1.+ 0.*x( testbd);
    u(~testbd,1) =       0.*x(~testbd);
    u(:,2) =  0.*x;      
    end
    % Derivative of the exact solution: Could be anything here because no
    % exact solution is known.
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1,1) =  0.*y;
    uprime(:,1,2) =  0.*x;
    uprime(:,2,1) =  0.*y;
    uprime(:,2,2) =  0.*x;
    end
   % Diffusion coefficient tensor in case of scalar. This is the viscosity
   % of the fluid in the Stokes solver.
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function coef =  Diff1(p)
    x = p(:,1); y = p(:,2);
    coef = 0*x+1;
    end
   % Exact soltion of p
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  pp(p)
    x = p(:,1); y = p(:,2);
    rhs =0.*x;
    end
%
pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du, ...
            'viscosity',@Diff1,'Df',@Df,'pp',@pp, 'KnownSol', KnownSol );
end