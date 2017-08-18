function [dNz] = dN1Solid(xi,ne)
% -------------------------------------------------------------------------
% File: dN1Solid.m
%   Computes the Gradient matrix for a 2/3 nodes of linear element.
%
%  Nodes numbering:
%   1---(3)---2 
%
% Syntax:
%   [dNz] = dN1Solid(xi,ne)
%
% Input:
%   xi   : First normalized coordinate.
%    ne  : Nodes number.
%
% Output:
%   dNz  : Element gradient matrix
%
% Date: 31/10/2014
%   Version 3.0   
% Updated by: Nico De Marchi
% -------------------------------------------------------------------------

switch ne 
    case 2
        % Define shape function derivatives
        dNz1 = -1/2; 
        dNz2 =  1/2;
        % Gradient matrix     
        dNz = [dNz1 dNz2];
    case 3    
        % Define shape function derivatives
        dNz1 = xi-1/2; 
        dNz2 = xi+1/2;
        dNz3 = -2*xi;
        % Gradient matrix   
        dNz = [dNz1 dNz2 dNz3 ];
   
end

end