function pde =  Poisson_MyPDE
% Poisson PDE information needed:
%
% 'KnownSol'  KnownSol: 
%             = 1 --> Know exact solution and would like errors be computed. 
%             = 0 --> Do not know exact solution and no error computation
%'f',        @f : Load data (right hand side function)
%'Df',       @Df: 1st order partial derivative of the right hand side function
%'g_D',      @g_D :  Dirichlet boundary condition
%'exactu',   @exactu : Exact solution of u, exact solution of u (needed only
%                      if you want errors be computed.
%'Du',       @Du : 1st order partial derivatives of the exact solution u.
%                  Used for computing the gradient error if KnownSol=1.
%
%  Copyright (C) Yujie LIU. Junping WANG. See COPYRIGHT.txt for details.
%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Set the solution indicator
    %
       KnownSol = 0;    % No exact solution is necessary.
    %  KnownSol = 1;  % Yes, I know the exact solution and would like
    %                   the error be computed.
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Enter load data (right hand side function)
    function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    rhs =  2*pi^2*sin(pi*x).*sin(pi*y);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Enter the derivative of load vector
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) =  2*pi^3*cos(pi*x).*sin(pi*y);
    rhs(:,2) =  2*pi^3*sin(pi*x).*cos(pi*y);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Enter the Dirichlet boundary condition
    function u =  g_D(p)
    u =  exactu(p);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %  Enter the exact solution
    %    IMPORTANT: Leave this part alone if you do not have 
    %               the exact solution. Do not remove the following functions
    % 
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u =  sin(pi*x).*sin(pi*y);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of the exact solution
    %    IMPORTANT: Leave this part alone if you do not have 
    %               the exact solution. Do not remove the following functions
    % 
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1) = pi*cos(pi*x).*sin(pi*y);
    uprime(:,2) = pi*sin(pi*x).*cos(pi*y);
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pde = struct('f',@f, 'Df',@Df,'exactu',@exactu,'g_D',@g_D,'Du',@Du,...
                 'KnownSol', KnownSol);
end