function pde = Diffusion_convection_MyPDE
% Convection Diffusion PDE information needed:
%
% 'KnownSol'  KnownSol: 
%             = 1 --> Know exact solution and would like errors be
%             computed. 
%             = 0 --> Do not know exact solution and no error computation
%
% 'd',        @Diff: Diffusion coeficient tensor '\alpha' entered as a function
%
%'beta',      @beta: Vector coefficient beta entered as a function
%
%'Dbeta',     @Dbeta: Derivative of beta entered as a function
%
%'c',         @c: Coefficient c entered as a function
%
%'Dc',        @Dc: Derivative of c entered as a function
%
% 'f',        @f : load function data (right hand side function)
%
% 'Df',       @Df: 1st order partial derivative of the right hand side 
%                  function
%
%'g_D',       @g_D:  Dirichlet boundary data
%
%'exactu',    @exactu : exact solution of u (needed only if you want errors
%                      be computed.
%'Du',        @Du : 1st order partial derivatives of the exact solution u.
%                  Used for computing the gradient error if KnownSol=1.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Copyright (C)  Yujie LIU. Junping WANG. See COPYRIGHT.txt for details.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Set the solution indicator
%
   KnownSol=0;    % No exact solution is necessary.
%  KnownSol = 1;  % Yes, I know the exact solution and would like
%                   the error be computed.
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Enter the Diffusion coefficient tensor: alpha
    function coef =  Diff(p)
    x = p(:,1); y = p(:,2);
    coef(:,1,1) = abs(sin(8*x).*cos(10*y))+1;
    coef(:,1,2) = 0*x-0.5;
    coef(:,2,1) = 0*x-0.5;
    coef(:,2,2) = abs(cos(8*x).*cos(10*y))+1;
    end
    %Enter the vector coefficient beta 
    function rhs =  beta(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) =  0.*x + 1.;
    rhs(:,2) =  0.*x + 1.;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Enter the derivative of beta
    function rhs =  Dbeta(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1,1) =  0.*x + 0.;
    rhs(:,1,2) =  0.*x + 0.;
    rhs(:,2,1) =  0.*x + 0.;
    rhs(:,2,2) =  0.*x + 0.;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Enter the coefficient c
    function rhs =  c(p)
    x = p(:,1); y = p(:,2);
    rhs =  0.*x + 1.;
    end 
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Enter the derivative of c
    function rhs =  Dc(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) =  0.*x + 0.;
    rhs(:,2) =  0.*x + 0.;
    end 
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Enter the load data (right hand side function)
    function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    rhs =  y + x +x.*y+5*sin(5*x).*cos(5*y);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Enter the derivative of load data:f
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) = 1. + y;
    rhs(:,2) = 1. + x;
    end  
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Enter the Dirichlet boundary condition
    function u =  g_D(p)
    u =  exactu(p);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % exact solution
    %    IMPORTANT: Leave this part alone if you do not have 
    %               the exact solution. Do not remove the following functions
    % 
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u =  x.*y+3.0*sin(3*x).*cos(4*y);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Derivative of the exact solution
    %    IMPORTANT: Leave this function alone if you do not have 
    %               the exact solution. Do not delete this function from the
    %               code
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1) = y;
    uprime(:,2) = x;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pde = struct('beta',@beta,'Dbeta',@Dbeta,'c',@c,'Dc',@Dc,'f',@f, 'Df',@Df,...
             'exactu',@exactu,'g_D',@g_D,'Du',@Du,'d',@Diff ,'KnownSol', KnownSol);
end