function pe = PeSolid(Xe,pvet,nip,cos)
%--------------------------------------------------------
% PeSolid:
%   Creates the force vector of an 4/6/8 nodes of quadrilateral element or
%   2/3 node linear element
%
% Syntax:
%   pe = PeSolid(Xe,pvet,nip)
%
% Input:
%   Xe   : Element nodal coordinate array.
% pvet   : Pressure vector.
%  nip   : Number of Gauss points.
%
% Output:
%   pe   : Internal force vector.
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

% Initialize internal force vector
pe = zeros(ne*dof,1);
Ne=zeros(1,ne);

% loop reduction for 1D integration
if dof == 2
    nipj=1;
else
    nipj=nip;
end 

% Gauss integration of internal force vector:
for i = 1:nip
    for j = 1:nipj
        
            % Get element Shape function matrix and Jacobian
            N = NSolid(ne,dof-1,dof,r(i),r(j));
            
            for l= 1:ne
                    Ne(1,l)=N(1,dof*(l-1)+1);
            end
                
            Xgp=(Ne*Xe);
            
            if dof == 2
                % Jacobian matrix
                Jt = ASolid(Xe,r(i));
                if pvet(2,1)==0
                    dl=[Jt(1,2) -Jt(1,1)]';
                else  
                    dl=Jt';
                end
                % Element internal force vector
                pe = pe + w(i)*(N'*(pvet(1,1)+(Xgp*cos(:,1+pvet(2,1)))'))*dl;
            else
                % Determinant of Jacobian matrix
                [~,da] = ASolid(Xe,r(i),r(j));
                
                % Element internal force vector
                pe = pe + w(i)*w(j)*( N'*(pvet+(Xgp*cos)') )*da;
            end
    end
end
end