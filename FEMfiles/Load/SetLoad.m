function p = SetLoad(p,P)
%--------------------------------------------------------
% SetLoad: 
%   Sets the nodal loads in global vector form. 
%
% Syntax:
%   p = SetLoad(p,P)
%
% Input:
%   p    : Initial global load vector.
%   P    : Loads given as [NodeNo Load(1) .. Load(dof)].
%
% Output:
%   p    : Loads in global column vector. 
%
% Date:
%   Version 1.0   15.10.12
%--------------------------------------------------------

% Determine number of dof per node
dof = size(P,2) - 1;

% Put loads into global load vector.
for i = 1:size(P,1)
  j = (P(i,1)-1)*dof; 
  p(j+1:j+dof) = P(i,1+1:1+dof)';
end
