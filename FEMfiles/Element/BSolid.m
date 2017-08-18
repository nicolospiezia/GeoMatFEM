function [Jt,B,dN,FT,fduT,G] = BSolid(X0e,xi,eta,zeta,ue,due)
% -------------------------------------------------------------------------
% BSolid:
%   Computes element Jacobian matrix Jt, element strain 
%   interpolation matrix B, shape function derivatives 
%   and deformation gradient at Gauss points for an hyperelastoplastic 
%   8/12/16/20 nodes solid hexahedral element or 4/8 nodes quadrilateral
%   element in plane strain.
%
% Syntax:
%   [Jt,B,dN,FT,fuT,G] = BSolid(X0e,xi,eta,zeta,ue,due)
%
% Input:
%   Xe   : Element nodal coordinate array.
%   xi   : First normalized coordinate.
%   eta  : Second normalized coordinate.
%   zeta : Third normalized coordinate.
%   ue   : Element nodal displacements.
%   due  : Element nodal incremental displacements.
%
% Output:
%   Jt   : Element Jacobian matrix.
%   B    : Element strain interpolation matrix.
%   dN   : Derivative of shape functions.
%   FT   : Transpose of deformation gradient tensor. 
%   fduT : Transpose of Update deformation gradient tensor.
%    G   : Derivative of shape functions (for large strain computation).
%
% Date: 31/10/2014
%   Version 3.0   
% Created by: Nicolò Spiezia
% Updated by: Nico De Marchi
% -------------------------------------------------------------------------

% Number of nodes per element and dof's per node
ne  = size(X0e,1);
dof = size(X0e,2);
%
if nargin == 4
    ue=zeros(ne*dof,1);
    due=ue;
end
% reshape displacement vectors
Ue = reshape(ue,size(X0e,2),size(X0e,1))';
DUe= reshape(due,size(X0e,2),size(X0e,1))';
Xe = (X0e + Ue + DUe);

% Gradient matrix
if dof == 2
    dNz = dN2Solid(xi,eta,ne);
else
    dNz = dN3Solid(xi,eta,zeta,ne);
end

% Element Jacobian matrix
Jt0= dNz*X0e;
Jt = dNz*Xe;
Jt1= dNz*(X0e+Ue);

% Transform to physical coordinates
dN0 = Jt0\dNz;
dN  = Jt\dNz;
dN1 = Jt1\dNz;

% Deformation gradient 

FT  = dN0*Xe;
fduT= dN1*Xe;

% Define unit matrix.
I = eye(dof);

if dof == 2 % Plate element
    
    % Add FT(3,3)= fduT(3,3) = 1
    FT = [FT(1,1)   FT(1,2)   0
          FT(2,1)   FT(2,2)   0
          0         0         1];
      
    fduT = [fduT(1,1)   fduT(1,2)   0
            fduT(2,1)   fduT(2,2)   0
            0           0           1];
        
    % Assemble strain interpolation matrix B 
    B = zeros(3,ne*dof);
    for i = 1:ne                % node number
        j = 2*i-1;              % column number
        B(:,j:j+1) = [ dN(1,i)    0
                        0        dN(2,i)
                       dN(2,i)   dN(1,i) ];
    end
    
    % Assemble strain interpolation matrix G for large strain
    if nargout == 6
        
        G = [ dN(1,1)*I  dN(1,2)*I  dN(1,3)*I  dN(1,4)*I
              dN(2,1)*I  dN(2,2)*I  dN(2,3)*I  dN(2,4)*I ]; 

        if ne >= 8
            Gcen = [ dN(1,5)*I  dN(1,6)*I  dN(1,7)*I  dN(1,8)*I 
                     dN(2,5)*I  dN(2,6)*I  dN(2,7)*I  dN(2,8)*I ]; 

            G =[G Gcen];
        end
        
        if ne == 9
            G9 = [ dN(1,9)*I   
                   dN(2,9)*I ]; 

            G =[G G9];
        end
    end
else
% solid element    
% Assemble strain interpolation matrix B at Gauss point
    B = zeros(6,ne*dof);
        for i = 1:ne                    % node number
            j = 3*i-2;                  % column number
            B(:,j:j+2) = [ dN(1,i)   0        0
                           0         dN(2,i)  0
                           0         0        dN(3,i)
                           dN(2,i)   dN(1,i)  0
                           0         dN(3,i)  dN(2,i)
                           dN(3,i)   0        dN(1,i)];
        end
   
        % Assemble strain interpolation matrix G for large strain    
        if nargout == 6

G = [ dN(1,1)*I dN(1,2)*I dN(1,3)*I dN(1,4)*I dN(1,5)*I dN(1,6)*I dN(1,7)*I dN(1,8)*I
      dN(2,1)*I dN(2,2)*I dN(2,3)*I dN(2,4)*I dN(2,5)*I dN(2,6)*I dN(2,7)*I dN(2,8)*I  
      dN(3,1)*I dN(3,2)*I dN(3,3)*I dN(3,4)*I dN(3,5)*I dN(3,6)*I dN(3,7)*I dN(3,8)*I ];

            if ne >= 12
                G12 = [ dN(1,9)*I dN(1,10)*I dN(1,11)*I dN(1,12)*I 
                        dN(2,9)*I dN(2,10)*I dN(2,11)*I dN(2,12)*I   
                        dN(3,9)*I dN(3,10)*I dN(3,11)*I dN(3,12)*I ];

                G =[G G12];
            end

            if ne >= 16
                G16 = [ dN(1,13)*I dN(1,14)*I dN(1,15)*I dN(1,16)*I 
                        dN(2,13)*I dN(2,14)*I dN(2,15)*I dN(2,16)*I   
                        dN(3,13)*I dN(3,14)*I dN(3,15)*I dN(3,16)*I ];

                G =[G G16];
            end

            if ne >= 20
                G20 = [ dN(1,17)*I dN(1,18)*I dN(1,19)*I dN(1,20)*I 
                        dN(2,17)*I dN(2,18)*I dN(2,19)*I dN(2,20)*I   
                        dN(3,17)*I dN(3,18)*I dN(3,19)*I dN(3,20)*I ];

                G =[G G20];
            end
        end
end

end
