function pde = sinsindata
%% This is a pde data for the Poisson equation
%
%     f = 2*pi^2*sin(pi*x)*sin(pi*y);
%     u = sin(pi*x)*sin(pi*y);
%     Du = (pi*cos(pi*x)*sin(pi*y), pi*sin(pi*x)*cos(pi*y));
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%'f',        @f : Load data (right hand side function)
%'exactu',   @exactu : Exact solution of u
%'g_D',      @g_D :  Dirichlet boundary condition
%'Du',       @Du : Derivative of the exact solution u
%'Df',       @Df: Derivative of the right hand side function
%
%  Copyright (C) Yujie LIU. Junping WANG. See COPYRIGHT.txt for details.
%
    % Exact solution of u
    KnownSol=1; % set value 0 if do not know the exact solution.
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % load data (right hand side function)
    function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    rhs =  2*pi^2*sin(pi*x).*sin(pi*y);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %  Derivative of load vector
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) =  2*pi^3*cos(pi*x).*sin(pi*y);
    rhs(:,2) =  2*pi^3*sin(pi*x).*cos(pi*y);
    end   
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % exact solution
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u =  sin(pi*x).*sin(pi*y);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Dirichlet boundary condition
    function u =  g_D(p)
    u =  exactu(p);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of the exact solution
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1) = pi*cos(pi*x).*sin(pi*y);
    uprime(:,2) = pi*sin(pi*x).*cos(pi*y);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    pde = struct('f',@f, 'Df',@Df,'exactu',@exactu,'g_D',@g_D,'Du',@Du ,...
                 'KnownSol', KnownSol);
end