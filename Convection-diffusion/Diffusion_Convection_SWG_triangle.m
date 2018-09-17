function [u,AD,b,ers,info] = ...
    Diffusion_Convection_SWG_triangle(node,elem,pde,bdFlag,option)
% This is to solve the convection diffusion equation on triangle partition.
%
%   u = Diffusion_Convection_SWG_triangle(node,elem,pde,bdFlag,option) 
%   produces the linear SWG approximation of the Convection-Diffusion equation
% 
%       -div(\alpha \nabla u) + \beta \cdot \nabla u + c u = f, 
%        u = g on \partial \Omega
% 
%   The mesh is given by node and elem and the boundary edge is given by
%   bdFlag. See meshdoc, bddoc for details. The data is given by the
%   structure pde which contains function handles f, g_D or d.
%   
%   The function Diffusion_Convection_SWG_triangle assembes the matrix 
%   equation AD*u = b and solves it by the direct solver (small size <= 2e3).
%   The Dirichlet boundary condition is built into the matrix AD and is build into b.
%
% Copyright (C)  Yujie LIU.  Junping WANG. See COPYRIGHT.txt for details.
%
NT = size(elem,1);
N  = size(node,1);
h  = sqrt(2)/sqrt(NT);
tic;

elem_area=zeros(NT,1);
elem_center=zeros(NT,2);
elem_edge_normal=zeros(NT,3,2);
elem_edge_length=zeros(NT,3,1);

[elem2edge,edge] = dofedge(elem);
NE= size(edge,1);

Ndof = NE ;
A    = sparse(Ndof,Ndof);
b    = zeros(Ndof,1);
elem2elem =(1:NT)';

elem2dofu   = elem2edge;

for i=1:NT
    NV=size(elem(i,:),2);
    nodal=[];
    for j=2:NV
        nodal = [nodal; node(elem(i,j),:)];
    end 
    nodal = [nodal; node(elem(i,1),:)];
    Ndof_local    = 2*NV+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Get the element information
%
% Element = struct ('Edge', edge_info, 'Vertex', nodal, 'NumberEdge', ...
%                   num_edge, 'NumberVertex', num_vertex);
     num_edge=NV;
     num_vertex=NV;
% 
% get edge information
%
    for i_edge=1:num_edge-1
         MyEdge_info(i_edge,1)=i_edge;
         MyEdge_info(i_edge,2)=i_edge+1;
    end
         MyEdge_info(num_edge,1)=num_vertex;
         MyEdge_info(num_edge,2)=1;
%
% Get normal direction for each edge --- point outward of the domain
% Note that this information is supposed to be an input to the function
%
   for i_edge=1:num_edge
% Get the index for left and righ end points
        left_pt=MyEdge_info(i_edge,1);
        right_pt=MyEdge_info(i_edge,2);
% Get the tangent vector from left to right with lengh |e_i|
        edge_tangent_vec(i_edge,1)=nodal(right_pt, 1)-nodal(left_pt,1);
        edge_tangent_vec(i_edge,2)=nodal(right_pt,2)-nodal(left_pt,2);  
 % Get edge length
    mag_edge(i_edge) = sqrt(edge_tangent_vec(i_edge,1)^2 + edge_tangent_vec(i_edge,2)^2);
% 
%  Get outward normal vector (should be taken an an input to the function.
%  Need to make sure the vector points away from the polygon.
%
% safety check to make sure the edge is non-degenerate
%
   if mag_edge(i_edge) < 1.0e-14
    fprintf('WARMING --- WARNING ... WARNING \n');
    fprintf('One of the edges of the polygon is degenerate. Please check your \n');
    fprintf('polygonal partition ... \n');
   end
%
% The following assumes that the nodal points are in counterclockwise direction
% For clockwise nodal points, please change the direction of the
% following vector
%
    MyUnitNormal(i_edge,1)= edge_tangent_vec(i_edge,2)/mag_edge(i_edge);
    MyUnitNormal(i_edge,2)=-edge_tangent_vec(i_edge,1)/mag_edge(i_edge);
   end
%   
%
MyElement = struct('EdgeInfo', MyEdge_info, 'EdgeOrientation', MyUnitNormal,...
                   'VertexCoordinates', nodal, 'NumberOfEdge', num_edge);
%
% END OF Element Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
     elem_edge_normal(i,1:3,:)=MyUnitNormal(1:3,:);
     elem_edge_length(i,1:3,1)  =mag_edge(1:3);

     Ele_Stiffness = zeros(Ndof_local,Ndof_local);
     Ele_Load      = zeros(Ndof_local,1);
     
     [Ele_Stiffness, Ele_Load, Polygon_Area, Polygon_Center]=...
     Element_Stiffness_Matrix_Diffusion_Convection(MyElement,pde,option);
    
     elem_center(i,1) = Polygon_Center(1);
     elem_center(i,2) = Polygon_Center(2);  
     elem_area(i) = Polygon_Area;

     ii =double( elem2dofu(i,:))';
     A_localOne   = sparse(Ndof,Ndof);
     for iii=1:NV
         for jjj=1:NV
         A_localOne = A_localOne + sparse(ii(iii),ii(jjj),Ele_Stiffness(iii,jjj),Ndof,Ndof);
         end
     end
     A = A + A_localOne;
     b(ii) = b(ii) + Ele_Load(1:NV);
end


[AD,b,u,freeDof,isPureNeumann] = getbdCCFVcoef_triangle(A,b);

%% Record assembling time
assembleTime = toc;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeDof), return; end
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if NE <= 1e3  % Direct solver for small size systems
        option.solver = 'direct';
    else         
    fprintf('WARMING --- WARNING : Solve the system of linear equations \n');
    fprintf('size of the system too large !\n');
    fprintf(' ... \n');  
    end
end

solver = option.solver;
% solve
switch solver
    case 'direct'
        tic;
        %fprintf('Solver is Ok!\n')
        u(freeDof) = AD(freeDof,freeDof)\b(freeDof);
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);            
end


%% Output information
eqn = struct('A',AD,'b',b,'edge',edge,'freeDof',freeDof);
info.assembleTime = assembleTime;

%% Compute DU
DU = zeros(NT,2);
for i=1:3
DU(:,1)= DU(:,1) + u(elem2edge(:,i)).*elem_edge_normal(:,i,1).*elem_edge_length(:,i);
DU(:,2)= DU(:,2) + u(elem2edge(:,i)).*elem_edge_normal(:,i,2).*elem_edge_length(:,i);    
end

DU(:,1)= DU(:,1)./elem_area;
DU(:,2)= DU(:,2)./elem_area;

%% Show solution
middle_points = 0.5*(node(edge(:,1),:) + node(edge(:,2),:));

close all;

uI =zeros(N,1);
uI = pde.exactu(node);
DUIc= pde.Du(elem_center);
uIc = pde.exactu(middle_points);

info.Mid_Edge = middle_points;
info.Elem_center= elem_center;
info.NE=NE;
info.error=u-uIc;

if (option.exact==1)
figure(114);
subplot(1,3,1)
scatter3(middle_points(:,1),middle_points(:,2),u,30,u,'filled'), hold on;
title('Numerical solution u_h');
subplot(1,3,2)
showsolution(node,elem,uI ,[-62,58]);
title('Exact solution u');
subplot(1,3,3)
scatter3(middle_points(:,1),middle_points(:,2),u-uIc,30,u-uIc,'filled'); 
title('Error u_h-u ');

u1I =griddata(middle_points(:,1),middle_points(:,2),u(1:NE),node(:,1),node(:,2),'v4');

figure(115);
subplot(1,2,1)
h3 = patch('Faces', elem, 'Vertices', node);
set(h3,'FaceVertexCData',u1I,'FaceColor','interp','edgecolor','k');
view(2); axis equal; axis tight; axis off;
colorbar
hold off;
title('Numerical solution u_h');
subplot(1,2,2)
h3 = patch('Faces', elem, 'Vertices', node);
set(h3,'FaceVertexCData',uI(:,1),'FaceColor','interp','edgecolor','k');
view(2); axis equal; axis tight; axis off;
colorbar
hold off;
title('Exact solution u');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif (option.exact==0)
figure(114);
scatter3(middle_points(:,1),middle_points(:,2),u,30,u,'filled'), hold on;
title('Numerical solution u_h');

u1I =griddata(middle_points(:,1),middle_points(:,2),u(1:NE),node(:,1),node(:,2),'v4');
figure(115);
h3 = patch('Faces', elem, 'Vertices', node);
set(h3,'FaceVertexCData',u1I,'FaceColor','interp','edgecolor','k');
view(2); axis equal; axis tight; axis off;
colorbar
hold off;
title('Numerical solution u_h');    
end
error=h*h*((u-uIc).*(u-uIc));
err=sqrt(sum(error));

errordu=h*h*((DU(:,1)-DUIc(:,1)).*(DU(:,1)-DUIc(:,1)) ...
        + (DU(:,2)-DUIc(:,2)).*(DU(:,2)-DUIc(:,2)));
errdu=sqrt(sum(errordu));

ers.L2=err;
ers.H1=errdu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdWG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AD,b,u,freeDof,isPureNeumann]= getbdCCFVcoef_triangle(A,b)
    %% GETBDCR Boundary conditions for Poisson equation: WG element.
    
    u =zeros(Ndof,1);
    %% Initial check
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end

    %% Part 1: Modify the matrix for Dirichlet condition
    % Find Dirichlet boundary nodes: fixedNode
        fixedNode = []; freeNode = [];
    if ~isempty(bdFlag) % find boundary edges and boundary nodes
        % fprintf('Boundary node is Ok!\n')
        [fixedNode,bdEdge,isBdNode] = findboundary(elem,bdFlag);
        freeNode = find(~isBdNode);
        freeEdge = [];fixeEdge =[];
        idxD = (bdFlag(:) == 1);
        isDirichlet = false(NE,1);
        isDirichlet(elem2edge(idxD)) = true;
        fixEdge = edge(isDirichlet,:);  % Robin edges  
        freeEdge =edge(~isDirichlet,:);
        fixeEdgeIndex =find(isDirichlet);
        freeEdgeIndex =find(~isDirichlet);
    end
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
        % no bdFlag, only pde.g_D is given
        [fixedNode,bdEdge,isBdNode] = findboundary(elem);
        freeNode = find(~isBdNode);
    end
    isPureNeumann = false;    

% Modify the matrix
        if ~isempty(fixedNode)
        bdidx = zeros(Ndof,1); 
        bdidx(fixeEdgeIndex) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
    else
        AD = A;
    end
    
    %% Part 2: Find boundary edges and modify the right hand side b
     % Dirichlet boundary condition
    if isnumeric(pde.g_D) && all(pde.g_D == 0)   % zero g_D
        pde.g_D = [];
    end
    if ~isPureNeumann && ~isempty(fixedNode) && ~isempty(pde.g_D)
        if isnumeric(pde.g_D)  % pde.g_D could be a numerical array  
            
            %fprintf('Boundary node fixe is Ok!\n')
        else % pde.g_D is a function handle
            u_fix = pde.g_D(0.5*(node(fixEdge(:,1),:) + node(fixEdge(:,2),:)));
            u(fixeEdgeIndex) = u_fix(:,1); 
            %fprintf('Boundary node fixe is Ok!\n')
        end
        b = b - A*u;
    end
    if ~isPureNeumann % non-empty Dirichlet boundary condition
        b(fixeEdgeIndex) = u(fixeEdgeIndex);
    end
    
    freeDof = [freeEdgeIndex];
end

end