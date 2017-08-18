function A = Assem(A,Ae,Te,dof)
%--------------------------------------------------------
% Assem: 
%   Assembles system matrix or vector by adding 
%   element contributions from elements to an 
%   existing global matrix. The system matrix may be 
%   square like the stiffness mass or conductivity 
%   matrix, or may be a one-dimensional vector like 
%   the load vector.
%
% Syntax:
%   A = Assmk(A,Ae,Te)
%   A = Assmk(A,Ae,Te,dof)
%
% Input:
%   A  :  Global matrix.
%   Ae :  Element matrix.
%   Te :  Element topology vector.
%   dof:  Degrees of freedom per node.
%
% Output:
%   A  :  Updated global matrix. 
%
% Date:
%   Version 1.0   13.09.12
%--------------------------------------------------------

% Number of element nodes
ne = size(Te,2) - 1; 

% Define global address vector for element dofs
ig = zeros(1,ne*dof);
for i = 1:ne
  ig([(i-1)*dof+1:i*dof]) = [(Te(i)-1)*dof+1:Te(i)*dof]; 
end

% Add element matrix/vector to global matrix/vector
if size(A,2) ~= size(A,1) 
  A(ig) = A(ig) + Ae;
else
  A(ig,ig) = A(ig,ig) + Ae;
end
     
