function V = VSolid(V,T,X,G,av,nip,flag)
%--------------------------------------------------------
% VSolid:
%   Creates and assembles the volumetric force vector
%   of solid elastic hexahedral elements or quadrilateral 
%   elements in plane strain.
%
% Syntax:
%   V = VSolid(V,T,X,G,vol,nip)
%
% Input:
%   V    : Existing internal force vector.
%   T    : System topology array.
%   X    : System nodal coordinate array.
%   G    : Element property array.
%   av   : Acceleration vector.
%  nip   : Number of Gauss points.
%
% Output:
%   V    : Volumetric force vector.
%
% Date: 31/10/2014
%   Version 2.0   
% Updated by: Nico De Marchi
%--------------------------------------------------------

% Number of nodes per element and dof's per node
ne  = size(T,2)-1;
dof = size(X,2);
nelem = size(T,1);

% Initialize arrays
Ve = zeros(ne*dof,nelem);

% Loop over elements
parfor j = 1:nelem
    
    % Define element arrays
    Te = T(j,:); 
    
    if flag == Te(end)
        Xe = X(Te(1:ne),:);
        Ge = G(Te(end),:);
        
        % Evaluate element volumetric force vector
        Ve(:,j) = VeSolid(Xe,Ge,av,nip);
    end
    
end
 % Assemble into global arrays
for j = 1:nelem
   Te = T(j,:);
   V = Assem(V,Ve(:,j),Te,dof);
end
end
