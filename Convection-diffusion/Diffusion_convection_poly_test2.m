function pde = Diffusion_convection_poly_test2
%% PDE data for model problem convection diffusion equation
%              u    = 3x^2+2xy
%             \alpha= 2 0
%                     0 1
%             \beta = 1
%                     1
%                 c = 1;
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
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   % Exact solution of u
    KnownSol=1; % set value 0 if do not know the exact solution.
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %Vector coefficient beta 
    function rhs =  beta(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) =  0.*x + 1.;
    rhs(:,2) =  0.*x + 1.;
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
    rhs =  -12 + 2*y + 8*x + 2*x.*y + 3*x.*x;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of load data:f
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) = 8. + 2*y + 6*x;
    rhs(:,2) = 2. + 2*x;
    end   
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % exact solution
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u =  3*x.*x + 2*x.*y;
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
    uprime(:,1) = 6*x + 2*y;
    uprime(:,2) = 2*x;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Diffusion coefficient tensor: alpha
    function coef =  Diff(p)
    x = p(:,1); y = p(:,2);
    coef(:,1,1) = 0*x+2.;
    coef(:,1,2) = 0*x;
    coef(:,2,1) = 0*x;
    coef(:,2,2) = 0*x+1;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pde = struct('beta',@beta,'Dbeta',@Dbeta,'c',@c,'Dc',@Dc,'f',@f, 'Df',@Df,...
             'exactu',@exactu,'g_D',@g_D,'Du',@Du,'d',@Diff,'KnownSol', KnownSol);
end