function Ve = VeSolid(Xe,Ge,av,nip)
%--------------------------------------------------------
% VeSolid:
%   Creates the volumetric force vector of an 8/12/16/20 nodes elastic 
%   hexahedral element or 4/8 nodes quadrilateral element in plane strain.
%
% Syntax:
%   Ve = VeSolid(Xe,EPe,nip)
%
% Input:
%   Xe   : Element nodal coordinate array.
%   Ge   : Element property vector.
%   av   : Acceleration vector.
%  nip   : Number of Gauss points.
%
% Output:
%   Ve   : Volumetric force vector.
%
% Date: 31/10/2014
%   Version 3.0   
% Updated by: Nico De Marchi
%--------------------------------------------------------

% Number of nodes per element and dof's per node
ne  = size(Xe,1);
dof = size(Xe,2);

% Gauss abscissae and weights.
[r,w] = Gauss(nip);

% Material parameters
rhos = Ge(2);
n = Ge(15);
rhow = Ge(16);

% Element Effective density
rho=(1-n)*rhos+(n-1)*rhow;

% Initialize volumetric force vector
Ve = zeros(ne*dof,1);


% loop reduction for 2D integration
if dof == 2
    nipk=1;
    wk=ones(1,nip);
else
    nipk=nip;
    wk=w;
end

% Gauss integration of internal force vector:
for i = 1:nip
    for j = 1:nip
        for k = 1:nipk

            % Get element shape function matrix
            N = NSolid(ne,dof,dof,r(i),r(j),r(k));
            
            % Get element Jacobian matrix
            Jt = BSolid(Xe,r(i),r(j),r(k));
            
            % Element volumetric force vector
            Ve = Ve + w(i)*w(j)*wk(k)*(N'*av)*rho*det(Jt);
            
        end
    end
end

end
