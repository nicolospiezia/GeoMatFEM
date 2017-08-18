function [r,w,p] = Gauss(n,a,b)
%--------------------------------------------------------
% Gauss: 
%   Evaluates Gauss points and weights for numerical 
%   integration using Gauss quadrature. 
%
% Syntax:
%   [r w p] = Gauss(n)
%   [r w p] = Gauss(n,a,b)
%
% Input:
%   n    : Number of integration points.
%   a    : Element material data.
%   b    : Element displacement vector.
%
% Output: 
%   r    : Vectors of gauss points.
%   w    : Vector of weights.
%   p    : Order of polynomium to be accurate 
%
% Date:
%   Version 2.0   31.10.14
%--------------------------------------------------------
% Order of polynomium which is integrated exact
p = 2*n-1;
% Gauss points and weights
if p <= 1 
    w = 2;
    r = 0; 
elseif p <= 3 
    w = [1 1];
    r = [-1 1]/sqrt(3);
elseif p <= 5 
    w = [5/9 8/9 5/9];
    r = [-sqrt(3/5) 0 sqrt(3/5)];
elseif p <= 7 
    w = [18+sqrt(30) 18+sqrt(30) 18-sqrt(30) 18-sqrt(30)]/36;
    r = [sqrt((3-2*sqrt(6/5))/7)  -sqrt((3-2*sqrt(6/5))/7) ...
         sqrt((3+2*sqrt(6/5))/7) -sqrt((3+2*sqrt(6/5))/7)      ];
else 
    fprintf('\n Gauss points and weights for polynomials of order 8\n or higher are not implemented yet \n')
end
% Change of interval
if nargin ==3
w = (b-a)/2*w;
r = (b-a)/2*r + (b+a)/2;
end