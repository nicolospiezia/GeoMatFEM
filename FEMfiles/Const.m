function [iu,ic,id,uimp] = Const(C,dof,ndof)
%--------------------------------------------------------------------------
% IdenBC:
%   Identify constrained and free dof's
%
% Syntax:
%  [iu,ic] = IdenBC(BC,dof,ndof)
%
% Input:
%   C    :  Constraint matrix.
%   dof  :  Number of nodal dofs.
%   ndof :  total number of dofs.  
%
% Output:
%   iu   : Unconstrained dofs.
%   ic   : Constrained dofs.
%
% Date: 31/10/2014
%   Version 2.0   
% Created by: Nicolò Spiezia
% Updated by: Nico De Marchi
%--------------------------------------------------------------------------
dim=size(C,1);
id=ones(1,1);
uimp=zeros(1,1);
c = 0;

% Constrained dofs
ic = (C(:,1)-1)*dof + C(:,2);

for i=1:dim
    if C(i,3) ~= 0
     c=c+1;
        id(c,1) = (C(i,1)-1)*dof + C(i,2);
        uimp(c,1) = C(i,3);
    end
end

% Unconstrained dofs
iu = setdiff(1:ndof,ic);