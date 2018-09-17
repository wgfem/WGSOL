function pde = Stokes_sincos
%% SINCOSDATA trigonometric  data for Poisson equation
%
% u1 = @(x,y)sin(x).*sin(x).*cos(y).*sin(y);
% u2 = @(x,y)-cos(x).*sin(x).*sin(y).*sin(y);
% gu1x= @(x,y)2*sin(x).*cos(x).*cos(y).*sin(y);
% gu1y= @(x,y)sin(x).*sin(x).*(cos(y).*cos(y)-sin(y).*sin(y));
% gu2x= @(x,y)(sin(x).*sin(x)-cos(x).*cos(x)).*sin(y).*sin(y);
% gu2y= @(x,y)-2*cos(x).*sin(x).*sin(y).*cos(y);
% p  = @(x,y)cos(x).*cos(y);
% f1 = @(x,y)-sin(x).*cos(y)-cos(x).*cos(x).*sin(2*y)+3*sin(x).*sin(x).*sin(2*y);
% f2 = @(x,y)-cos(x).*sin(y)+cos(y).*cos(y).*sin(2*x)-3*sin(y).*sin(y).*sin(2*x);
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%'f',        @f : Load data (right hand side function)
%'exactu',   @exactu : Exact solution of u
%'g_D',      @g_D :  Dirichlet boundary condition
%'Du',       @Du : Derivative of the exact solution u 
%'viscosity',@Diff1 : Diffusion coefficient tensor in case of scalar
%'Df',       @Df: Derivative of the right hand side function
%'pp',       @pp: Exact solution of p
%
%  Copyright (C)  Yujie LIU. Junping WANG. See COPYRIGHT.txt for details.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
    % load data (right hand side function)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1) =-sin(x).*cos(y)-cos(x).*cos(x).*sin(2*y)+3*sin(x).*sin(x).*sin(2*y);
    rhs(:,2) =-cos(x).*sin(y)+cos(y).*cos(y).*sin(2*x)-3*sin(y).*sin(y).*sin(2*x);
    rhs(:,3) = x*0.0;
    end
% Dirichlet boundary condition
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function u =  g_D(p)
    u =  exactu(p);
    end
    % exact solution of u
 KnownSol= 1; % set this =0 if no exact solution is known
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u(:,1) = sin(x).*sin(x).*cos(y).*sin(y);
    u(:,2) =-cos(x).*sin(x).*sin(y).*sin(y);
    end
    
    % Derivative of the exact solution u
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1,1) = 2*sin(x).*cos(x).*cos(y).*sin(y);
    uprime(:,1,2) = sin(x).*sin(x).*(cos(y).*cos(y)-sin(y).*sin(y));
    uprime(:,2,1) = (sin(x).*sin(x)-cos(x).*cos(x)).*sin(y).*sin(y);
    uprime(:,2,2) =-2*cos(x).*sin(x).*sin(y).*cos(y);
    end
    % Diffusion coefficient tensor in case of scalar: This is the viscosity
    % for STOKES test
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function coef =  Diff1(p)
    x = p(:,1); y = p(:,2);
    coef = 0*x+1;
    end
    % Derivative of the right hand side function
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    rhs(:,1,1) =-cos(x).*cos(y)+2*cos(x).*sin(x).*sin(2*y)+6*sin(x).*cos(x).*sin(2*y);
    rhs(:,1,2) = sin(x).*sin(y)-2*cos(x).*cos(x).*cos(2*y)+6*sin(x).*sin(x).*cos(2*y);
    rhs(:,2,1) = sin(x).*sin(y)+2*cos(y).*cos(y).*cos(2*x)-6*sin(y).*sin(y).*cos(2*x);
    rhs(:,2,2) =-cos(x).*cos(y)-2*cos(y).*sin(y).*sin(2*x)-6*sin(y).*cos(y).*sin(2*x);
    rhs(:,3,1) = 0.*x;
    rhs(:,3,2) = 0.*x;
    end
    % Exact solution of p
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  pp(p)
    x = p(:,1); y = p(:,2);
    rhs =cos(x).*cos(y);
    end
%
pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du, ...
             'viscosity',@Diff1,'Df',@Df,'pp',@pp, 'KnownSol', KnownSol);
end