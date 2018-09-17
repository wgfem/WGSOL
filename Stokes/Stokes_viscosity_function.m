function pde = Stokes_viscosity_function
%% This is a model problem pde data for the stokes equation:
%
%    u1        = @(x,y) sin(x).*cos(y);
%    u2        = @(x,y)-cos(x).*sin(y);
%    gu1x      = @(x,y) cos(x).*cos(y);
%    gu1y      = @(x,y)-sin(x).*sin(y);
%    gu2x      = @(x,y) sin(x).*sin(y);
%    gu2y      = @(x,y)-cos(x).*cos(y);;
%    p         = @(x,y) x.*x-y.*y;
%    viscosity = @(x,y) xy+1
%    f1        = @(x,y) 2(1+xy)*sin(x).*cos(y)+2*x +xsin(x)sin(y) - ycos(x)cos(y);
%    f2        = @(x,y)-2(1+xy)*cos(x).*sin(y)-2*y -y sin(x)sin(y) +xcos(x)cos(y);
%    f3        = @(x,y) y-x;
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
    % load data (right hand side function)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  f(p)
    % Get the point
    x = p(:,1); y = p(:,2);
    % define the first component of f in the momentum equation
    rhs(:,1) =-y.*cos(x).*cos(y) +2*(x.*y+1).*sin(x).*cos(y)+2*x ...
        + x.*sin(x).*sin(y);
    % define the second component of f in the momentum equation
    rhs(:,2) = x.*cos(x).*cos(y) -2*(x.*y+1).*cos(x).*sin(y)-2*y ...
        - y.*sin(x).*sin(y);
    % define the RHS function in the mass conservation equation. Should be
    % zero for incompressible fluid. But the code takes any function. Make
    % sure that Dirichlet boundary data is compatible with this function
    rhs(:,3) = 0.*x;
    end
    % Derivative of the right hand side function
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  Df(p)
    x = p(:,1); y = p(:,2);
    %rhs(:,1) =-y.*cos(x).*cos(y) +2*(x.*y+1).*sin(x).*cos(y)+2*x ...
    %    + x.*sin(x).*sin(y);
    rhs(:,1,1) =  y.*sin(x).*cos(y) + 2*y.*sin(x).*cos(y) + 2*(x.*y+1).*cos(x).*cos(y) + 2. ...
        + sin(x).*sin(y) + x.*cos(x).*sin(y);
    rhs(:,1,2) = -cos(x).*cos(y)+y.*cos(x).*sin(y) + 2*x.*sin(x).*cos(y) ...
        - 2*(x.*y+1).*sin(x).*sin(y);
    %rhs(:,2) = x.*cos(x).*cos(y) -2*(x.*y+1).*cos(x).*sin(y)-2*y ...
    %    - y.*sin(x).*sin(y);
    rhs(:,2,1) = cos(x).*cos(y) - x.*sin(x).*cos(y) - 2*y.*cos(x).*sin(y) ...
        + 2*(x.*y+1).*sin(x).*sin(y);
    rhs(:,2,2) =-x.*cos(x).*sin(y)- 2*x.*cos(x).*sin(y) - 2*(x.*y+1).*cos(x).*cos(y) ...
        -2;
    rhs(:,3,1) =  0.*x;
    rhs(:,3,2) =  0.*x;
    end
    % exact solution of u
    KnownSol=1; % set value =0 if do not know the exact solution
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u(:,1) = sin(x).*cos(y);
    u(:,2) =-cos(x).*sin(y);
    end
    % Dirichlet boundary condition
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function u =  g_D(p)
    u =  exactu(p);
    end
    % Derivative of the exact solution: Need this in Stokes for Computing
    % Gradient Error. If no exact solution is known, then this part could
    % be anything.
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1,1) =  cos(x).*cos(y);
    uprime(:,1,2) = -sin(x).*sin(y);
    uprime(:,2,1) =  sin(x).*sin(y);
    uprime(:,2,2) = -cos(x).*cos(y);
    end
       % Diffusion coefficient tensor in case of scalar: VISCOSITY in STOKES
   % Please set the viscosity below
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function viscosity =  Diff1(p)
    x = p(:,1); y = p(:,2);
    viscosity = x.*y+1;
    end
   % Exact soltion for pressure p
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function rhs =  pp(p)
    x = p(:,1); y = p(:,2);
    rhs =x.*x-y.*y;
    end
%
pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du, ...
             'viscosity',@Diff1,'Df',@Df,'pp',@pp, 'KnownSol', KnownSol);
end