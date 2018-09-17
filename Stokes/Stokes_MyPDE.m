function pde = Stokes_MyPDE
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Stokes PDE information needed:
%
% 'KnownSol'  KnownSol: 
%             = 1 --> Know exact solution and would like errors be
%             computed. 
%             = 0 --> Do not know exact solution and no error computation
%
% 'f',        @f : load function data (right hand side function)
% 'Df',       @Df: 1st order partial derivative of the right hand side 
%                  function
%
%'g_D',      @g_D:  Dirichlet boundary data
%
%'viscosity',@Diff1: fluid viscosity entered as a function or scalar
%
%'exactu',   @exactu : exact solution of u (needed only if you want errors
%                      be computed.
%'Du',       @Du : 1st order partial derivatives of the exact solution u.
%                  Used for computing the gradient error.
%'pp',       @pp:  exact soltion of the pressure p (needed only if
%                  KnownSol=1.
%
%  Copyright (C)  Yujie LIU. Junping WANG. See COPYRIGHT.txt for details.
%  2018
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Set the solution indicator
%
   KnownSol=0;    % No exact solution is necessary.
%  KnownSol = 1;  % Yes, I know the exact solution and would like
%                   the error be computed.
% 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Enter viscosity 
%
    function viscosity =  Diff1(p)
    x = p(:,1); y = p(:,2);
    viscosity = 0*x+1;
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% load data (right hand side function)
%
    function rhs =  f(p)
% Get the point (x,y)
    x = p(:,1); y = p(:,2);
% Enter the first component of f=(f1,f2,f3) in the momentum equation
    rhs(:,1) =2*sin(x).*cos(y)+2*x;
% Enter the second component of f=(f1,f2,f3) in the momentum equation
    rhs(:,2) =-2*cos(x).*sin(y)-2*y;
% Enter the RHS function in the mass conservation equation. Should be
% zero for incompressible fluid. But the code takes any function. Make
% sure that Dirichlet boundary data is compatible with this function
    rhs(:,3) = y-x;
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Enter the 1st order partial derivatives of the right hand side function. 
% This is used to compute the load vector
%
   function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1,1) =  2*cos(x).*cos(y)+2;
    rhs(:,1,2) = -2*sin(x).*sin(y);
    rhs(:,2,1) =  2*sin(x).*sin(y);
    rhs(:,2,2) = -2*cos(x).*cos(y)-2;
    rhs(:,3,1) = -1+0.*x;
    rhs(:,3,2) =  1+0.*x;
   end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Enter the Dirichlet boundary data
%
    function u =  g_D(p)
    u =  exactu(p);
    end
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Enter the exact solution: 
%    IMPORTANT: Leave this part alone if you do not have 
%               the exact solution. Do not remove the following functions
% 
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u(:,1) = sin(x).*cos(y)+x.*y;
    u(:,2) =-cos(x).*sin(y)-x.*y;
    end
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%1st order derivative of the exact solution. Needed for computing the 
%gradient
%
%    IMPORTANT: Leave this function alone if you do not have 
%               the exact solution. Do not delete this function from the
%               code
%
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1,1) =  cos(x).*cos(y)+y;
    uprime(:,1,2) = -sin(x).*sin(y)+x;
    uprime(:,2,1) =  sin(x).*sin(y)-y;
    uprime(:,2,2) = -cos(x).*cos(y)-x;
    end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Enter the exact soltion for pressure p
%    IMPORTANT: Leave this function alone if you do not have 
%               the exact solution. But do not delete this function.
% 
    function rhs =  pp(p)
    x = p(:,1); y = p(:,2);
    rhs =x.*x-y.*y;
    end
%
pde = struct('f',@f,'Df',@Df,'g_D',@g_D,'viscosity',@Diff1, ...
    'exactu',@exactu,'Du',@Du,'pp',@pp, 'KnownSol', KnownSol);
end