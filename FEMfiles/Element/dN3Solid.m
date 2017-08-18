function [dNz] = dN3Solid(xi,eta,zeta,ne)
% -------------------------------------------------------------------------
% File: dN3Solid.m
%   Computes the Gradient matrix for a 8/12/16/20 nodes of solid 
%   hexahedral element
%
% Syntax:
%   [dNz] = dN3Solid(xi,eta,zeta,ne)
%
% Input:
%   xi   : First normalized coordinate.
%   eta  : Second normalized coordinate.
%  zeta  : Third normalized coordinate.
%    ne  : Nodes number.
%
% Output:
%   dNz  : Element gradient matrix
%
% Date: 31/10/2014
%   Version 3.0   
% Updated by: Nico De Marchi
%--------------------------------------------------------

% Define shape function derivatives
dN1xi = -(eta-1)*(zeta-1);    dN2xi =  (eta-1)*(zeta-1);
dN3xi = -(eta+1)*(zeta-1);    dN4xi =  (eta+1)*(zeta-1);
dN5xi =  (eta-1)*(zeta+1);    dN6xi = -(eta-1)*(zeta+1);
dN7xi =  (eta+1)*(zeta+1);    dN8xi = -(eta+1)*(zeta+1);

dN1eta = -(xi-1)*(zeta-1);    dN2eta =  (xi+1)*(zeta-1);
dN3eta = -(xi+1)*(zeta-1);    dN4eta =  (xi-1)*(zeta-1);
dN5eta =  (xi-1)*(zeta+1);    dN6eta = -(xi+1)*(zeta+1);
dN7eta =  (xi+1)*(zeta+1);    dN8eta = -(xi-1)*(zeta+1);

dN1zeta = -(xi-1)*(eta-1);    dN2zeta =  (xi+1)*(eta-1);
dN3zeta = -(xi+1)*(eta+1);    dN4zeta =  (xi-1)*(eta+1);
dN5zeta =  (xi-1)*(eta-1);    dN6zeta = -(xi+1)*(eta-1);
dN7zeta =  (xi+1)*(eta+1);    dN8zeta = -(xi-1)*(eta+1);

% Define derivative matrix at gauss points
dNz = [dN1xi   dN2xi   dN3xi   dN4xi   dN5xi   dN6xi   dN7xi   dN8xi
       dN1eta  dN2eta  dN3eta  dN4eta  dN5eta  dN6eta  dN7eta  dN8eta
       dN1zeta dN2zeta dN3zeta dN4zeta dN5zeta dN6zeta dN7zeta dN8zeta ]/8;

if ne >= 12
    
    % Quadratic terms for midside nodes
    dN9xi  = -1/2*xi*(-1+zeta)*(-1+eta);
    dN10xi =  1/2*xi*(-1+zeta)*(1+eta);
    dN11xi = -1/2*xi*(1+zeta)*(1+eta);
    dN12xi =  1/2*xi*(1+zeta)*(-1+eta);
           
    dN9eta   = -1/4*(-1+xi)*(1+xi)*(-1+zeta);
    dN10eta  =  1/4*(-1+xi)*(1+xi)*(-1+zeta);
    dN11eta  = -1/4*(-1+xi)*(1+xi)*(1+zeta);
    dN12eta  =  1/4*(-1+xi)*(1+xi)*(1+zeta);
        
    dN9zeta   = -1/4*(-1+xi)*(1+xi)*(-1+eta);
    dN10zeta  =  1/4*(-1+xi)*(1+xi)*(1+eta);
    dN11zeta  = -1/4*(-1+xi)*(1+xi)*(1+eta);
    dN12zeta  =  1/4*(-1+xi)*(1+xi)*(-1+eta);
    
    dN12 = [dN9xi     dN10xi     dN11xi     dN12xi
            dN9eta    dN10eta    dN11eta    dN12eta
            dN9zeta   dN10zeta   dN11zeta   dN12zeta];
    
     % Modify corner nodes
     dNz(:,1) = dNz(:,1) - 1/2*dN12(:,1);
     dNz(:,2) = dNz(:,2) - 1/2*dN12(:,1);
     dNz(:,3) = dNz(:,3) - 1/2*dN12(:,2);
     dNz(:,4) = dNz(:,4) - 1/2*dN12(:,2);
     dNz(:,5) = dNz(:,5) - 1/2*dN12(:,4);
     dNz(:,6) = dNz(:,6) - 1/2*dN12(:,4);
     dNz(:,7) = dNz(:,7) - 1/2*dN12(:,3);
     dNz(:,8) = dNz(:,8) - 1/2*dN12(:,3);
     
     % Expand gradient matrix
    dNz = [dNz dN12];
end

if ne >= 16
    
    % Quadratic terms for midside nodes
    dN13xi =  1/4*(-1+eta)*(1+eta)*(-1+zeta);
    dN14xi = -1/4*(-1+eta)*(1+eta)*(-1+zeta);
    dN15xi =  1/4*(-1+eta)*(1+eta)*(1+zeta);    
    dN16xi = -1/4*(-1+eta)*(1+eta)*(1+zeta);
    
    dN13eta  =  1/2*eta*(-1+zeta)*(1+xi);
    dN14eta  = -1/2*eta*(-1+zeta)*(-1+xi);
    dN15eta  =  1/2*eta*(1+zeta)*(-1+xi);
    dN16eta  = -1/2*eta*(1+zeta)*(1+xi);
    
    dN13zeta  =  1/4*(-1+eta)*(1+eta)*(1+xi);
    dN14zeta  = -1/4*(-1+eta)*(1+eta)*(-1+xi);
    dN15zeta  =  1/4*(-1+eta)*(1+eta)*(-1+xi);
    dN16zeta  = -1/4*(-1+eta)*(1+eta)*(1+xi);
    
    dN16 = [dN13xi     dN14xi    dN15xi     dN16xi
            dN13eta    dN14eta   dN15eta    dN16eta
            dN13zeta   dN14zeta  dN15zeta   dN16zeta];
        
     % Modify corner nodes
     dNz(:,1) = dNz(:,1) - 1/2*dN16(:,2);
     dNz(:,2) = dNz(:,2) - 1/2*dN16(:,1);
     dNz(:,3) = dNz(:,3) - 1/2*dN16(:,1);
     dNz(:,4) = dNz(:,4) - 1/2*dN16(:,2);
     dNz(:,5) = dNz(:,5) - 1/2*dN16(:,3);
     dNz(:,6) = dNz(:,6) - 1/2*dN16(:,4);
     dNz(:,7) = dNz(:,7) - 1/2*dN16(:,4);
     dNz(:,8) = dNz(:,8) - 1/2*dN16(:,3);
     
     % Expand gradient matrix
    dNz = [dNz dN16];
    
end

if ne >= 20
    
    % Quadratic terms for midside nodes
    dN17xi = -1/4*(-1+zeta)*(1+zeta)*(-1+eta);
    dN18xi =  1/4*(-1+zeta)*(1+zeta)*(-1+eta);
    dN19xi = -1/4*(-1+zeta)*(1+zeta)*(1+eta);
    dN20xi =  1/4*(-1+zeta)*(1+zeta)*(1+eta);
    
    dN17eta  = -1/4*(-1+zeta)*(1+zeta)*(-1+xi);
    dN18eta  =  1/4*(-1+zeta)*(1+zeta)*(1+xi);
    dN19eta  = -1/4*(-1+zeta)*(1+zeta)*(1+xi);
    dN20eta  =  1/4*(-1+zeta)*(1+zeta)*(-1+xi);
    
    dN17zeta  = -1/2*zeta*(-1+xi)*(-1+eta);
    dN18zeta  =  1/2*zeta*(1+xi)*(-1+eta);
    dN19zeta  = -1/2*zeta*(1+xi)*(1+eta);
    dN20zeta  =  1/2*zeta*(-1+xi)*(1+eta);
    
    dN20 = [dN17xi    dN18xi     dN19xi     dN20xi
            dN17eta   dN18eta    dN19eta    dN20eta
            dN17zeta  dN18zeta   dN19zeta   dN20zeta];
                
     % Modify corner nodes
     dNz(:,1) = dNz(:,1) - 1/2*dN20(:,1);
     dNz(:,2) = dNz(:,2) - 1/2*dN20(:,2);
     dNz(:,3) = dNz(:,3) - 1/2*dN20(:,3);
     dNz(:,4) = dNz(:,4) - 1/2*dN20(:,4);
     dNz(:,5) = dNz(:,5) - 1/2*dN20(:,1);
     dNz(:,6) = dNz(:,6) - 1/2*dN20(:,2);
     dNz(:,7) = dNz(:,7) - 1/2*dN20(:,3);
     dNz(:,8) = dNz(:,8) - 1/2*dN20(:,4);
     
     % Expand gradient matrix
    dNz = [dNz dN20];
end
end