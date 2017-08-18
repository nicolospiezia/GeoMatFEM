function p = PSolid(p,T,X,pres,nip,cos)
%--------------------------------------------------------
% PSolid:
%   Creates and assembles the force vector
%   of solid elastic hexahedral elements or quadrilateral
%   element in plane strain.
%
% Syntax:
%   p = PSolid(p,T,X,pres,nip)
%
% Input:
%   p    : Existing force vector.
%   T    : System topology array.
%   X    : System nodal coordinate array (3D).
%   
%  pvet  : pressure vector.
%  nip   : Number of Gauss points.
%
% Output:
%   p    : Force vector.
%
% Date: 31/10/2014
%   Version 3.0   
% Updated by: Nico De Marchi
%--------------------------------------------------------

% Number of nodes per element and dof's per node
n = size(pres,1);
ne = size(T,2)-1;
dof = size(X,2);

if nargin == 5
    cos=zeros(dof,dof);
end

% Initialize arrays
pe = zeros(ne*dof,n);

% Loop load elements
parfor j = 1:n
    % element & load face
    i = pres(j,1);
    face = pres(j,2);
    
    % Define element arrays
    Te= T(i,1:ne);
    Xe = X(Te,:);
    pvet = pres(j,3:end)'; 
    
    % Define 2D element or 1D element
    [X2,T2] = Order(Xe,face);
    
    % Evaluate element force vector
    p2 = PeSolid(X2,pvet,nip,cos);
    pe(:,j)=Assem(pe(:,j),p2,T2,dof);
    
end

     % Assemble into global arrays
for j = 1:n
    i= pres(j,1);
    p = Assem(p,pe(:,j),T(i,:),dof);
end
end