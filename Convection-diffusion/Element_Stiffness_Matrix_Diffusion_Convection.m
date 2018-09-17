function [Ele_Stiffness, Ele_Load, Polygon_Area, Polygon_Center]=...
          Element_Stiffness_Matrix_Diffusion_Convection(MyElement,pde,option)
% 
% Purpose:
%   Compute the element stiffness matrix for SWG applied to the convection- 
%   diffusion equation 
%              -div(\alpha \nabla u) + \beta \cdot \nabla u + c u = f, 
%                u = g on \partial \Omega
%
% Input:
%   MyElement(): a structure consisting the element information
%
% 1. MyElement.EdgeInfo(:)  edge information with MyElement.Edge(i,1) the 
% local label of the left end, and MyElement.Edge(i,2) the right end point.
%
% 2. MyElement.VertexCoordinates: length of (NumberOfVertex times 2). It 
% contains the coordinates of the i_th vertex point for i=1, 2, ... 
% NumberOfVertices. Note that NumerOfVertices is not a variable in the input.
%
% 3. MyElement.NumerOfEdge: It is the number of edge for this element
%
% 4. MyElement.EdgeOrientation(:) unit outward normal vector to the edge i.
%
% Copyright (C)  Junping WANG. Yujie LIU.  See COPYRIGHT.txt for details.
% 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Data for the current polygonal element
kappa=option.kappa; % kappa=4.0, this is the stabilization parameter in SWG
alpha=option.alpha; % alpha=-1, kappa * h^{alpha} is the factor of the stabilizer
visc_f =pde.d;      % coefficient tensor
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Number of edges % should be an input variable to the function
%
num_edge = MyElement.NumberOfEdge;
% 
% edge information: edge_info(i,1), edge_info(i,2) gives the local
% numbering of the left end and the righ end.
%
edge_info=MyElement.EdgeInfo;
%
%coordinates of nodal points
nodal=MyElement.VertexCoordinates;
%
% Get normal direction for each edge --- point outward of the domain
unitnormal_vec=MyElement.EdgeOrientation;
for i_edge=1:num_edge
    art = sqrt(unitnormal_vec(i_edge,1)^2 + unitnormal_vec(i_edge,2)^2);
% safety check to make sure the edge is non-degenerate
%
if abs(art-1) > 1.0e-10
    fprintf('WARMING --- WARNING ... WARNING \n');
    fprintf('Your edge orientation/normal is not a unit vector. \n');
    fprintf('Applied normailization ... \n');
    unitnormal_vec(i_edge,:)=unitnormal_vec(i_edge,:)/art;
end
end
%
for i_edge=1:num_edge
% Get the index for left and righ end points
    left_pt=edge_info(i_edge,1);
    right_pt=edge_info(i_edge,2);
% Get the tangent vector from left to right with lengh |e_i|
    edge_tangent_vec(i_edge,1)=nodal(right_pt, 1)-nodal(left_pt,1);
    edge_tangent_vec(i_edge,2)=nodal(right_pt,2)-nodal(left_pt,2);
% Get edge length
mag_edge(i_edge) = sqrt(edge_tangent_vec(i_edge,1)^2 + edge_tangent_vec(i_edge,2)^2);
%
% safety check to make sure the edge is non-degenerate
%
if mag_edge(i_edge) < 1.0e-12
    fprintf('WARMING --- WARNING ... WARNING \n');
    fprintf('One of the edges of the polygon is degenerate. Please check your \n');
    fprintf('polygonal partition ... \n');
end
%
% Safety check to make sure the user input normal is indeed normal to the
% edge
safety_normal=edge_tangent_vec(i_edge,1)*unitnormal_vec(i_edge,1);
safety_normal=safety_normal+edge_tangent_vec(i_edge,2)*unitnormal_vec(i_edge,2);
if abs(safety_normal) > 1.0e-10
    fprintf('WARMING --- WARNING ... WARNING \n');
    fprintf('Data mis-match: your edge orientation (normal) is NOT NORMAL to the edge \n');
    fprintf('Please check your Element Data ... \n');
end
%%
% Get edge center coordinates
%
   edge_mid(i_edge,1) = (nodal(left_pt,1)+nodal(right_pt,1))/2;
   edge_mid(i_edge,2) = (nodal(left_pt,2)+nodal(right_pt,2))/2;
end
%
% Compute the center of gravity of the polygon
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
poly_xx_c=@(x,y) x^2;
poly_xy_c=@(x,y) x*y;

center_x=0.0;
for ii=1:num_edge
    ii_left=edge_info(ii,1);
    ii_right=edge_info(ii,2);
   fleft=poly_xx_c(nodal(ii_left,1), nodal(ii_left,2))/2;
   fright=poly_xx_c(nodal(ii_right,1), nodal(ii_right,2))/2;
   fmid=poly_xx_c(edge_mid(ii,1), edge_mid(ii,2))/2;
   center_x=center_x+mag_edge(ii)*(fleft+4.0*fmid+fright)*unitnormal_vec(ii,1)/6;
end
center_y=0.0;
for ii=1:num_edge
    ii_left=edge_info(ii,1);
    ii_right=edge_info(ii,2);
   fleft=poly_xy_c(nodal(ii_left,1), nodal(ii_left,2));
   fright=poly_xy_c(nodal(ii_right,1), nodal(ii_right,2));
   fmid=poly_xy_c(edge_mid(ii,1), edge_mid(ii,2));
   center_y=center_y+mag_edge(ii)*(fleft+4.0*fmid+fright)*unitnormal_vec(ii,1)/6;
end
%
%
% compute the area and the center of gravity of the polygon
%
area=0.0;
for ii=1:num_edge
    for jj=1:2
    area = area + mag_edge(ii)*edge_mid(ii,jj)*unitnormal_vec(ii,jj)/2;
    end
end
%
center(1)=center_x/area;
center(2)=center_y/area;
Polygon_Center=center;
%
% safety check to make sure the correctness of the polygon data
%
if area < 0.0
    fprintf('WARMING --- WARNING ... WARNING \n');
    fprintf('Area for the polygon is negative. Either the polygon is not in \n');
    fprintf('counterclockwise orientation or the polygon is invalid... \n');
end
Polygon_Area=area;
%
%fprintf(' Polygon center = %d .\n', center);
%for i=1:N
%fprintf(' Polygon edge centers = %d .\n', edge_mid(i,:));
%end
%
%fprintf(' Polygon area = %d .\n', area);
%for i=1:N
%fprintf(' edge normal vector = %d .\n', unitnormal_vec(i,:));
%end
%
% Get the diagonal matrix E from the edge length.
%
EdgeMatrix=diag(mag_edge);
%
% Get the M-Matrix, as in the paper
%
for ii=1:num_edge
    M_Matrix(ii,:)=[1, edge_mid(ii,1)-center(1), edge_mid(ii,2)-center(2)];
end
%
% Get inv(M^tEM)
%
MtEM_inv=inv(M_Matrix'*EdgeMatrix*M_Matrix);
%
% Get A-Matrix
%
A_Matrix = EdgeMatrix-EdgeMatrix*M_Matrix*MtEM_inv*M_Matrix'*EdgeMatrix;
%
% Get D-Matrix
%
D_Matrix = MtEM_inv*M_Matrix'*EdgeMatrix;
%
% Now we compute the polynom integrations involved. The integration formula is derived
% by using the Gauss' law. It is critical to compute some weights for each
% polynomial function
%
poly_xx=@(x,y) (x-center(1))^2;
poly_xy=@(x,y) (x-center(1))*(y-center(2));
poly_xxx=@(x,y) (x-center(1))^3;
poly_xyy=@(x,y) (x-center(1))*(y-center(2))^2;
poly_xxy=@(x,y) (x-center(1))^2*(y-center(2));
%

poly_xxxx=@(x,y)(x-center(1))^4;
poly_xxxy=@(x,y)(x-center(1))^3*(y-center(2));
poly_xxyy=@(x,y)(x-center(1))^2*(y-center(2))^2;
poly_xyyy=@(x,y)(x-center(1))*(y-center(2))^3;
% find the integration of each basis function on the polygon
%
% For \int_T (x-xc) dT
%
coef12=0.0;
%for ii=1:num_edge
%    ii_left=edge_info(ii,1);
%    ii_right=edge_info(ii,2);
%   fleft=poly_xx(nodal(ii_left,1), nodal(ii_left,2))/2;
%   fright=poly_xx(nodal(ii_right,1), nodal(ii_right,2))/2;
%   fmid=poly_xx(edge_mid(ii,1), edge_mid(ii,2))/2;
%   coef12=coef12+mag_edge(ii)*(fleft+4.0*fmid+fright)*unitnormal_vec(ii,1)/6;
%end
%
% For \int_T (y-yc) dT
%
coef13=0.0;
%for ii=1:num_edge
%    ii_left=edge_info(ii,1);
%    ii_right=edge_info(ii,2);
%   fleft=poly_xy(nodal(ii_left,1), nodal(ii_left,2));
%   fright=poly_xy(nodal(ii_right,1), nodal(ii_right,2));
%   fmid=poly_xy(edge_mid(ii,1), edge_mid(ii,2));
%   coef13=coef13+mag_edge(ii)*(fleft+4.0*fmid+fright)*unitnormal_vec(ii,1)/6;
%end
%
%
% For \int_T (x-xc)^2 dT
%
coef22=0.0;
for ii=1:num_edge
    ii_left=edge_info(ii,1);
    ii_right=edge_info(ii,2);
   fleft=poly_xxx(nodal(ii_left,1), nodal(ii_left,2))/3;
   fright=poly_xxx(nodal(ii_right,1), nodal(ii_right,2))/3;
   fmid=poly_xxx(edge_mid(ii,1), edge_mid(ii,2))/3;
   coef22=coef22+mag_edge(ii)*(fleft+4.0*fmid+fright)*unitnormal_vec(ii,1)/6;
end
%
% For \int_T (x-xc)(y-yc) dT
%
coef23=0.0;
for ii=1:num_edge
    ii_left=edge_info(ii,1);
    ii_right=edge_info(ii,2);
   fleft=poly_xxy(nodal(ii_left,1), nodal(ii_left,2))/2;
   fright=poly_xxy(nodal(ii_right,1), nodal(ii_right,2))/2;
   fmid=poly_xxy(edge_mid(ii,1), edge_mid(ii,2))/2;
   coef23=coef23+mag_edge(ii)*(fleft+4.0*fmid+fright)*unitnormal_vec(ii,1)/6;
end
%
% For \int_T (y-yc)^2 dT
%
coef33=0.0;
for ii=1:num_edge
    ii_left=edge_info(ii,1);
    ii_right=edge_info(ii,2);
   fleft=poly_xyy(nodal(ii_left,1), nodal(ii_left,2));
   fright=poly_xyy(nodal(ii_right,1), nodal(ii_right,2));
   fmid=poly_xyy(edge_mid(ii,1), edge_mid(ii,2));
   coef33=coef33+mag_edge(ii)*(fleft+4.0*fmid+fright)*unitnormal_vec(ii,1)/6;
end
%
% Polynom od degree higher than 3:
%
% For \int_T (x-xc)^4 dT
coefx22=0.0;
for ii=1:num_edge
    ii_left=edge_info(ii,1);
    ii_right=edge_info(ii,2);
    tier1=(nodal(ii_left,:)/3 + 2*nodal(ii_right,:)/3);
    tier2=(2*nodal(ii_left,:)/3 + nodal(ii_right,:)/3);
   fleft =poly_xxxx(nodal(ii_left,1), nodal(ii_left,2))/4;
   fright=poly_xxxx(nodal(ii_right,1), nodal(ii_right,2))/4;
   ftier1=poly_xxxx(tier1(1), tier1(2))/4;
   ftier2=poly_xxxx(tier2(1), tier2(2))/4;  
   coefx22=coefx22+mag_edge(ii)*(fleft+3.0*ftier1+3.0*ftier2+fright)*unitnormal_vec(ii,1)/8;
end
%
% For \int_T (x-xc)^3(y-yc) dT
coefx23=0.0;
for ii=1:num_edge
    ii_left=edge_info(ii,1);
    ii_right=edge_info(ii,2);
    tier1=(nodal(ii_left,:)/3 + 2*nodal(ii_right,:)/3);
    tier2=(2*nodal(ii_left,:)/3 + nodal(ii_right,:)/3);
   fleft =poly_xxxy(nodal(ii_left,1), nodal(ii_left,2))/3;
   fright=poly_xxxy(nodal(ii_right,1), nodal(ii_right,2))/3;
   ftier1=poly_xxxy(tier1(1), tier1(2))/3;
   ftier2=poly_xxxy(tier2(1), tier2(2))/3;  
   coefx23=coefx23+mag_edge(ii)*(fleft+3.0*ftier1+3.0*ftier2+fright)*unitnormal_vec(ii,1)/8;
end
%
% For \int_T (x-xc)^2(y-xc)^2 dT
coefx33=0.0;
for ii=1:num_edge
    ii_left=edge_info(ii,1);
    ii_right=edge_info(ii,2);
    tier1=(nodal(ii_left,:)/3 + 2*nodal(ii_right,:)/3);
    tier2=(2*nodal(ii_left,:)/3 + nodal(ii_right,:)/3);
   fleft =poly_xxyy(nodal(ii_left,1), nodal(ii_left,2))/2;
   fright=poly_xxyy(nodal(ii_right,1), nodal(ii_right,2))/2;
   ftier1=poly_xxyy(tier1(1), tier1(2))/2;
   ftier2=poly_xxyy(tier2(1), tier2(2))/2;  
   coefx33=coefx33+mag_edge(ii)*(fleft+3.0*ftier1+3.0*ftier2+fright)*unitnormal_vec(ii,1)/8;
end
%
% For \int_T (x-xc)(y-yc)^3 dT
coefy33=0.0;
for ii=1:num_edge
    ii_left=edge_info(ii,1);
    ii_right=edge_info(ii,2);
    tier1=(nodal(ii_left,:)/3 + 2*nodal(ii_right,:)/3);
    tier2=(2*nodal(ii_left,:)/3 + nodal(ii_right,:)/3);
   fleft =poly_xyyy(nodal(ii_left,1), nodal(ii_left,2));
   fright=poly_xyyy(nodal(ii_right,1), nodal(ii_right,2));
   ftier1=poly_xyyy(tier1(1), tier1(2));
   ftier2=poly_xyyy(tier2(1), tier2(2));  
   coefy33=coefy33+mag_edge(ii)*(fleft+3.0*ftier1+3.0*ftier2+fright)*unitnormal_vec(ii,1)/8;
end
% Set up the weight matrix
%
coef=[area, coef12, coef13; coef12, coef22, coef23; coef13, coef23, coef33];
% coef * (x-xT)
coefx=[coef12, coef22, coef23; coef22, coefx22, coefx23; coef23, coefx23, coefx33];
% coef * (y-yT)
coefy=[coef13, coef23, coef33; coef23, coefx23, coefx33; coef33, coefx33, coefy33];
%
% We are now ready to get the local load vector as written in the paper
%
% assume rhsf_one = rhsf_one (xc, yc) + rhsf_one_x (x-xc) + rhsf_one_y
% (y-y_c). Then we have
%

% Get B-Matrix
%
% visc be scalar
%B_Matrix = visc*EdgeMatrix*unitnormal_vec*unitnormal_vec'*EdgeMatrix/area;
% visc be function tensor
alpha1=visc_f(center);
alpha1M=zeros(2,2);
%
alpha1M(1,1)=alpha1(1,1,1);
alpha1M(1,2)=alpha1(1,1,2);
alpha1M(2,1)=alpha1(1,2,1);
alpha1M(2,2)=alpha1(1,2,2);
%
alpha1n= zeros(num_edge,2);
for ii=1:num_edge
alpha1n(ii,:) =  (alpha1M*unitnormal_vec(ii,:)')';
end
B_Matrix   = EdgeMatrix*alpha1n*unitnormal_vec'*EdgeMatrix/area;
%
% Get element stiffness matrix
%
h=sqrt(area);
T_M = kappa*h^alpha*A_Matrix+B_Matrix;
Ele_Stiffness = [T_M];
%
% Get the load function f (called ff, and its gradient DF, as defined in
% pdf
%
beta =pde.beta(center);
Dbeta=pde.Dbeta(center);

epa=mag_edge/area;

R_matrix=zeros(num_edge,num_edge);
for i=1:num_edge
nepa1=unitnormal_vec(i,1).*epa(i);
nepa2=unitnormal_vec(i,2).*epa(i);
beta_one_vec=[beta(:,1)'*nepa1, Dbeta(:,1,1)'*nepa1, Dbeta(:,1,2)'*nepa1];
beta_two_vec=[beta(:,2)'*nepa2, Dbeta(:,2,1) *nepa2, Dbeta(:,2,2)'*nepa2];
beta_vec    = beta_one_vec + beta_two_vec;
R_matrix(i,:)=beta_vec*coef*D_Matrix;
end

cc   =pde.c(center);
Dc   =pde.Dc(center);

C_matrix = zeros(num_edge,num_edge);
% case c=1
% C_matrix=D_Matrix'*coef*D_Matrix;
% case c function
C_matrix=cc(:,1)*D_Matrix'*coef*D_Matrix +  Dc(:,1)*D_Matrix'*coefx*D_Matrix + Dc(:,2)*D_Matrix'*coefy*D_Matrix;

Ele_Stiffness = Ele_Stiffness + R_matrix' + C_matrix;

ff=pde.f(center);
DF=pde.Df(center);
%
% Compute the load vector
%
f_one_vec=[ff(:,1)', DF(:,1)', DF(:,2)'];
%
% The following was for the desining test case
%f_one_vec=[rhsf_one(xc,yc), rhsf_one_x(xc, yc), rhsf_one_y(xc,yc)];
%
load_Fone=f_one_vec*coef*D_Matrix;
%
Ele_Load=[load_Fone'];
end
