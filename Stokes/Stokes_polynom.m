function pde = Stokes_polynom
%% This is a model problem pde data for the Stokes equation
%
% u1  = @(x,y)-256*x.*x.*(x-1).*(x-1).*y.*(y-1).*(2*y-1);
% u2  = @(x,y) 256*y.*y.*(y-1).*(y-1).*x.*(x-1).*(2*x-1);
% gu1x= @(x,y)-256*(4.*x.*x.*x-6*x.*x+2*x).*(2*y.*y.*y-3*y.*y+y);
% gu1y= @(x,y)-256*(x.*x.*x.*x-2*x.*x.*x+x.*x).*(6*y.*y-6*y+1);
% gu2x= @(x,y) 256*(y.*y.*y.*y-2*y.*y.*y+y.*y).*(6*x.*x-6*x+1);
% gu2y= @(x,y) 256*(4.*y.*y.*y-6*y.*y+2*y).*(2*x.*x.*x-3*x.*x+x);
% p   = @(x,y) 150*(x-0.5).*(y-0.5);
% f1  = @(x,y) 256.*(12*x.*x-12*x+2).*y.*(y-1).*(2*y-1)
%             +256*x.*x.*(x-1).*(x-1).*(12*y-6)+150*(y-0.5) ;
% f2  = @(x,y)-256.*(12*y.*y-12*y+2).*x.*(x-1).*(2*x-1)
%             -256*y.*y.*(y-1).*(y-1).*(12*x-6)+150*(x-0.5) ;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%'f',        @f : Load data (right hand side function)
%'exactu',   @exactu : Exact solution of u
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
% Exact solution of u
    KnownSol=1; % set value 0 if do not know the exact solution.

    % Load data (right hand side function)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) =  256.*(12*x.*x-12*x+2).*(2.*y.*y.*y-3*y.*y+y)...
               +256.*(x.*x.*x.*x -2.*x.*x.*x+x.*x).*(12*y-6)+150*(y-0.5);
    rhs(:,2) = -256.*(12*y.*y-12*y+2).*(2*x.*x.*x-3.*x.*x+x)...
               -256.*(y.*y.*y.*y-2.*y.*y.*y+y.*y).*(12*x-6)+150*(x-0.5);
    rhs(:,3) = 0.0*x; 
    end
    % Derivative of the right hand side function
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1,1) = 256.*(24*x-12).*(2.*y.*y.*y-3*y.*y+y)...
                +256*(4*x.*x.*x-6.*x.*x+2*x).*(12*y-6);          
    rhs(:,1,2) = 256.*(12*x.*x-12*x+2).*(6*y.*y-6*y+1)...
                +256.*(x.*x.*x.*x -2.*x.*x.*x+x.*x)*12 + 150;
    rhs(:,2,1) =-256.*(12*y.*y-12*y+2).*(6*x.*x-6*x+1)...
                -256*(y.*y.*y.*y-2.*y.*y.*y+y.*y)*12 + 150;
    rhs(:,2,2) =-256.*(24.*y -12).*(2*x.*x.*x-3.*x.*x+x)...
                -256*(4*y.*y.*y-6*y.*y+2*y).*(12*x-6);
    rhs(:,3,1) =0.*x ;
    rhs(:,3,2) =0.*x;
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u(:,1) = -256*x.*x.*(x-1).*(x-1).*y.*(y-1).*(2*y-1);
    u(:,2) =  256*y.*y.*(y-1).*(y-1).*x.*(x-1).*(2*x-1);
    end
    % Dirichlet boundary condition
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function u =  g_D(p)
    u =  exactu(p);
    end
    % Derivative of the exact solution u: Used for computing the gradient
    % error. No need to do anything with this part if no exact solution is
    % known.
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1,1) = -256*(4.*x.*x.*x-6*x.*x+2*x).*(2*y.*y.*y-3*y.*y+y);
    uprime(:,1,2) = -256*(x.*x.*x.*x-2*x.*x.*x+x.*x).*(6*y.*y-6*y+1);
    uprime(:,2,1) =  256*(y.*y.*y.*y-2*y.*y.*y+y.*y).*(6*x.*x-6*x+1);
    uprime(:,2,2) =  256*(4.*y.*y.*y-6*y.*y+2*y).*(2*x.*x.*x-3*x.*x+x);
    end
    % Diffusion coefficient tensor in case of scalar: THIS IS THE VISCOSITY
    % COEFFICIENT FOR STOKES
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function coef =  Diff1(p)
    x = p(:,1); y = p(:,2);
    coef = 0*x+1;
    end

    % Exact soltion of p
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  pp(p)
    x = p(:,1); y = p(:,2);
    rhs =  150*(x-0.5).*(y-0.5);
    end
%
pde = struct('f',@f,'Df',@Df,'exactu',@exactu,'g_D',@g_D,'Du',@Du, ...
             'viscosity',@Diff1,'pp',@pp, 'KnownSol', KnownSol);
end