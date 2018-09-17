function pde = Diffusion_convection_sinsindata_test4
%% PDE data for model problem convection diffusion equation
%
%       -div(\alpha \nabla u) + \beta \cdot \nabla u + c u = f, 
%        u = g on \partial \Omega
% with
%        u     = sin(\pi x) sin(\pi y), 
%
%        \alpha= xy+1 & 0   
%                0    & 3xy+1,
%
%        \beta = x^3y+xy+1
%                3x^2y+xy+2,
%
%             c= x^4y^2+xy+1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%'beta',     @beta: vector coefficient beta 
%'Dbeta',    @Dbeta: Derivative of beta
%'c',        @c: coefficient c
%'Dc',       @Dc: Derivative of c
%'f',        @f : load data (right hand side function)
%'Df',       @Df: Derivative of the right hand side function
%'exactu',   @exactu : exact solution of u
%'g_D',      @g_D :  Dirichlet boundary condition
%'Du',       @Du : Derivative of the exact solution u
%'d',        @Diff: Diffusion coefficient tensor: alpha
%
% Copyright (C)  Yujie LIU. See COPYRIGHT.txt for details.

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Exact solution of u
    KnownSol=1; % set value 0 if do not know the exact solution.
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %Vector coefficient beta 
    function rhs =  beta(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) =  1.*x.*x.*x.*y + 1.*y.*x +1.;
    rhs(:,2) =  3.*x.*x.*y + 1.*x.*y +2.;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of beta
    function rhs =  Dbeta(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1,1) =  3.*x.*x.*y + 1.*y + 0.;
    rhs(:,1,2) =  1.*x.*x.*x + 1.*x + 0.;
    rhs(:,2,1) =  6.*x.*y + y + 0.;
    rhs(:,2,2) =  3.*x.*x + x + 0.;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %Coefficient c
    function rhs =  c(p)
    x = p(:,1); y = p(:,2);
    rhs =  1.*x.*x.*x.*x.*y.*y + x.*y + 1.;
    end 
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of c
    function rhs =  Dc(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) =  4.*x.*x.*x.*y.*y + y + 0.;
    rhs(:,2) =  2.*x.*x.*x.*x.*y + 1.*x + 0.;
    end 
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %load data (right hand side function)
    function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    coef_beta = beta(p);
    coef_c=c(p);
    rhs =  pi^2*sin(pi*x).*sin(pi*y).*(4*x.*y +2.)...
         - y.*(pi*cos(pi*x).*sin(pi*y)) ...
         - 3*x.*(pi*sin(pi*x).*cos(pi*y)) ...
         + pi*coef_beta(:,1).*cos(pi*x).*sin(pi*y)...
         + pi*coef_beta(:,2).*sin(pi*x).*cos(pi*y)...
         + coef_c.*sin(pi*x).*sin(pi*y);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of load data:f
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    coef_beta  = beta(p);
    coef_Dbeta = Dbeta(p);
    coef_c     = c(p);
    coef_Dc    = Dc(p);
    rhs(:,1) =  pi^3*cos(pi*x).*sin(pi*y).*(4*x.*y+2.) + y.*4*pi^2*sin(pi*x).*sin(pi*y)... 
              + pi*y.*(pi*sin(pi*x).*sin(pi*y)) ...
              - 3*(pi*sin(pi*x).*cos(pi*y)) - 3*x.*(pi^2*cos(pi*x).*cos(pi*y)) ...
              - pi^2*coef_beta(:,1).*sin(pi*x).*cos(pi*y)  + pi^2*coef_beta(:,2).*cos(pi*x).*cos(pi*y) ...
              + pi*coef_Dbeta(:,1,1).*cos(pi*x).*sin(pi*y) + pi*coef_Dbeta(:,2,1).*sin(pi*x).*cos(pi*y) ...
              + pi*coef_c.*cos(pi*x).*sin(pi*y) + coef_Dc(:,1).*sin(pi*x).*sin(pi*y);
    rhs(:,2) =  pi^3*sin(pi*x).*cos(pi*y).*(4*x.*y+2.) + x.*4*pi^2*sin(pi*x).*sin(pi*y)...
              - (pi*cos(pi*x).*sin(pi*y)) - y.*(pi^2*cos(pi*x).*cos(pi*y))...
              + 3*x.*(pi^2*sin(pi*x).*sin(pi*y))...
              + pi^2*coef_beta(:,1).*cos(pi*x).*cos(pi*y) - pi^2*coef_beta(:,2).*sin(pi*x).*sin(pi*y) ...
              + pi*coef_Dbeta(:,1,2).*cos(pi*x).*sin(pi*y) + pi*coef_Dbeta(:,2,2).*sin(pi*x).*cos(pi*y)...
              + pi*coef_c.*sin(pi*x).*cos(pi*y) + coef_Dc(:,2).*sin(pi*x).*sin(pi*y);
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
    % Diffusion coefficient tensor: alpha
    function coef =  Diff(p)
    x = p(:,1); y = p(:,2);
    coef(:,1,1) = 1.*x.*y+1;
    coef(:,1,2) = 0*x;
    coef(:,2,1) = 0*x;
    coef(:,2,2) = 3.*x.*y+1;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pde = struct('beta',@beta,'Dbeta',@Dbeta,'c',@c,'Dc',@Dc,'f',@f, 'Df',@Df,...
             'exactu',@exactu,'g_D',@g_D,'Du',@Du,'d',@Diff,'KnownSol', KnownSol);
end