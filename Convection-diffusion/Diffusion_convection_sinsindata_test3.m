function pde = Diffusion_convection_sinsindata_test3
%% PDE data for model problem convection diffusion equation
%               u  = \sin(\pi x)\sin(\pi y) + x^2-y^2, 
%            \alpha= 1 0
%                    0 1
%             \beta= 1
%                    2
%                 c= 1
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
% Copyright (C)  Yujie LIU.  Junping WANG. See COPYRIGHT.txt for details.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   % Exact solution of u
    KnownSol=1; % set value 0 if do not know the exact solution.
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %Vector coefficient beta 
    function rhs =  beta(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) =  0.*x + 1.;
    rhs(:,2) =  0.*x + 2.;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of beta
    function rhs =  Dbeta(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1,1) =  0.*x + 0.;
    rhs(:,1,2) =  0.*x + 0.;
    rhs(:,2,1) =  0.*x + 0.;
    rhs(:,2,2) =  0.*x + 0.;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %Coefficient c
    function rhs =  c(p)
    x = p(:,1); y = p(:,2);
    rhs =  0.*x + 1.;
    end 
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of c
    function rhs =  Dc(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) =  0.*x + 0.;
    rhs(:,2) =  0.*x + 0.;
    end 
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %load data (right hand side function)
    function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    coef_beta = beta(p);
    coef_c=c(p);
    rhs =  2*pi^2*sin(pi*x).*sin(pi*y)+ pi*coef_beta(:,1).*(cos(pi*x).*sin(pi*y))+ 2*x ...
          +pi*coef_beta(:,2).*(sin(pi*x).*cos(pi*y))-4*y ...
          +coef_c.*(sin(pi*x).*sin(pi*y)+x.*x-y.*y);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of load data:f
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    coef_beta  = beta(p);
    coef_Dbeta = Dbeta(p);
    coef_c =c(p);
    coef_Dc=Dc(p);
    rhs(:,1) =  2*pi^3*cos(pi*x).*sin(pi*y) ...
    - pi^2*coef_beta(:,1).*sin(pi*x).*cos(pi*y) + pi^2*coef_beta(:,2).*cos(pi*x).*cos(pi*y)...
    + pi*coef_Dbeta(:,1,1).*cos(pi*x).*sin(pi*y)+ pi*coef_Dbeta(:,2,1).*sin(pi*x).*cos(pi*y)...
    + pi*coef_c.*cos(pi*x).*sin(pi*y) + coef_Dc(:,1).*sin(pi*x).*sin(pi*y)+2+2*x;
    rhs(:,2) =  2*pi^3*sin(pi*x).*cos(pi*y)...
    + pi^2*coef_beta(:,1).*cos(pi*x).*cos(pi*y) - pi^2*coef_beta(:,2).*sin(pi*x).*sin(pi*y)...
    + pi*coef_Dbeta(:,1,2).*cos(pi*x).*sin(pi*y)+ pi*coef_Dbeta(:,2,2).*sin(pi*x).*cos(pi*y)...
    + pi*coef_c.*sin(pi*x).*cos(pi*y) + coef_Dc(:,2).*sin(pi*x).*sin(pi*y)-4-2*y;
    end   
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % exact solution
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u =  sin(pi*x).*sin(pi*y) + x.*x - y.*y;
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
    uprime(:,1) = pi*cos(pi*x).*sin(pi*y) + 2*x;
    uprime(:,2) = pi*sin(pi*x).*cos(pi*y) - 2*y;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Diffusion coefficient tensor: alpha
    function coef =  Diff(p)
    x = p(:,1); y = p(:,2);
    coef(:,1,1) = 0*x+1;
    coef(:,1,2) = 0*x;
    coef(:,2,1) = 0*x;
    coef(:,2,2) = 0*x+1;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of alpha % x
    function coef =  Diffx(p)
    x = p(:,1); y = p(:,2);
    coef(:,1,1) = 0.*y+0.;
    coef(:,1,2) = 0*x;
    coef(:,2,1) = 0*x;
    coef(:,2,2) = 0.*y+0.;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of alpha % y
    function coef =  Diffy(p)
    x = p(:,1); y = p(:,2);
    coef(:,1,1) = 0.*x+0.;
    coef(:,1,2) = 0*x;
    coef(:,2,1) = 0*x;
    coef(:,2,2) = 0.*x+0.;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Diffusion coefficient tensor in case of scalar
    function coef =  Diff1(p)
    x = p(:,1); y = p(:,2);
    coef = 0*x+1;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pde = struct('beta',@beta,'Dbeta',@Dbeta,'c',@c,'Dc',@Dc,'f',@f, 'Df',@Df,...
             'exactu',@exactu,'g_D',@g_D,'Du',@Du,'d',@Diff,'KnownSol', KnownSol);
end