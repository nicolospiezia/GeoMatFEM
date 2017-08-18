function [Jt,da] = ASolid(X0e,xi,eta)
%--------------------------------------------------------
% File: ASolid.m
%   Computes the element Jacobian Jt at Gauss points for a 4/8 node 
%   for elastic quadrilateral element or 2/3 node linear element
%
%  Nodes numbering:
%
%   4---(6)---3 
%   |         | 
%  (8)       (7)
%   |         | 
%   1---(5)---2 
%
%
% Syntax:
%   [Jt,da] = ASolid(X0e,xi,eta)
%
% Input:
%   X0e  : Initial element nodal coordinate array.
%   xi   : First normalized coordinate.
%   eta  : Second normalized coordinate.
%
% Output:
%   Jt   : Jacobian matrix
%   da   : Determinant of Jacobian matrix
%
% Date: 31/10/2014
%   Version 3.0   
% Updated by: Nico De Marchi
%--------------------------------------------------------

% Number of nodes per element and dof's per node
ne  = size(X0e,1);
dof = size(X0e,2);

if dof == 2

    % Gradient matrix
    dNz = dN1Solid(xi,ne);
    
    % Element Jacobian matrix
    Jt = dNz*X0e;

    % % Determinant of Jacobian matrix
    da = sqrt(Jt(1,:)*Jt(1,:)');
    
else

    % Gradient matrix
    dNz = dN2Solid(xi,eta,ne);
    
    % Element Jacobian matrix
    Jt = dNz*X0e;

    % Determinant of Jacobian matrix for surface in 3D space
    da = sqrt((Jt(1,:)*Jt(1,:)')*(Jt(2,:)*Jt(2,:)')-(Jt(1,:)*Jt(2,:)')^2);
end

end