function pde = Stokes_test
%% This is a model problem pde data for the stokes equation:
%
%    u1  = @(x,y) sin(x).*cos(y)+x.*y;
%    u2  = @(x,y)-cos(x).*sin(y)-x.*y;
%    gu1x= @(x,y) cos(x).*cos(y)+y;
%    gu1y= @(x,y)-sin(x).*sin(y)+x;
%    gu2x= @(x,y) sin(x).*sin(y)-y;
%    gu2y= @(x,y)-cos(x).*cos(y)-x;;
%    p   = @(x,y) x.*x-y.*y;
%    f1  = @(x,y) 2*sin(x).*cos(y)+2*x;
%    f2  = @(x,y)-2*cos(x).*sin(y)-2*y;
%    f3  = @(x,y) y-x;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%'f',        @f : load data (right hand side function)
%'exactu',   @exactu : exact solution of u
%'g_D',      @g_D :  Dirichlet boundary condition
%'Du',       @Du : Derivative of the exact solution u
%'viscosity',@Diff1 : Diffusion coefficient tensor in case of scalar
%'Df',       @Df: Derivative of the right hand side function
%'pp',       @pp: Exact soltion of p
%
%  Copyright (C)  Yujie LIU. Junping WANG. See COPYRIGHT.txt for details.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%KnownSolution=0;


% Enter viscosity 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function viscosity =  Diff1(p)
    x = p(:,1); y = p(:,2);
    viscosity = 0*x+1;
    end
% load data (right hand side function)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  f(p)
    % Get the point
    x = p(:,1); y = p(:,2);
    % define the first component of f in the momentum equation
    rhs(:,1) =2*sin(x).*cos(y)+2*x;
    % define the second component of f in the momentum equation
    rhs(:,2) =-2*cos(x).*sin(y)-2*y;
    % define the RHS function in the mass conservation equation. Should be
    % zero for incompressible fluid. But the code takes any function. Make
    % sure that Dirichlet boundary data is compatible with this function
    rhs(:,3) = y-x;
    end
    % Derivative of the right hand side function
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1,1) =  2*cos(x).*cos(y)+2;
    rhs(:,1,2) = -2*sin(x).*sin(y);
    rhs(:,2,1) =  2*sin(x).*sin(y);
    rhs(:,2,2) = -2*cos(x).*cos(y)-2;
    rhs(:,3,1) = -1+0.*x;
    rhs(:,3,2) =  1+0.*x;
    end
% Dirichlet boundary condition
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function u =  g_D(p)
    u =  exactu(p);
    end
%
% Enter information on exact solution of u: 
%
% set KnownSol = 1 if you know the exact solution and would like errors be
% computed.
% set KnownSol = 0 if you do not know the exact solution 
   KnownSol = 1;
%
% Enter exact solution: DO NOT TOUCH THE FOLLOWING FUNCTION IF YOU HAVE NO
% EXACT SOLUTION
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u(:,1) = sin(x).*cos(y)+x.*y;
    u(:,2) =-cos(x).*sin(y)-x.*y;
    end
 % Derivative of the exact solution. Needed for computing the gradient
 % error. But this could be anything if no exact solution is known.
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1,1) =  cos(x).*cos(y)+y;
    uprime(:,1,2) = -sin(x).*sin(y)+x;
    uprime(:,2,1) =  sin(x).*sin(y)-y;
    uprime(:,2,2) = -cos(x).*cos(y)-x;
    end   
   % Exact soltion for pressure p
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  pp(p)
    x = p(:,1); y = p(:,2);
    rhs =x.*x-y.*y;
    end
%
pde = struct('f',@f,'Df',@Df,'g_D',@g_D,'viscosity',@Diff1, ...
    'exactu',@exactu,'Du',@Du,'pp',@pp, 'KnownSol', KnownSol);
end