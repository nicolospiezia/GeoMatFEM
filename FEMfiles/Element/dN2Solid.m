function [dNz] = dN2Solid(xi,eta,ne)
% -------------------------------------------------------------------------
% File: dN2Solid.m
%   Computes the Gradient matrix for a 4/6/8 nodes of elastic 
%   quadrilateral element.
%
%  Nodes numbering:
%
%   4---(6)---3 
%   |         | 
%  (8)  (9)  (7)
%   |         | 
%   1---(5)---2 
%
% Syntax:
%   [dNz] = dN2Solid(xi,eta,ne)
%
% Input:
%   xi   : First normalized coordinate.
%   eta  : Second normalized coordinate.
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
    case 4
        
        % Define shape function derivatives
        dNz1 = [ -(1-eta)
                -(1-xi) ]/4; 
        dNz2 = [ (1-eta)  
               -(1+xi) ]/4;
        dNz3 = [ (1+eta) 
                (1+xi) ]/4;
        dNz4 = [ -(1+eta)
                 (1-xi)]/4;
        % Gradient matrix     
        dNz = [dNz1 dNz2 dNz3 dNz4];
        
    case 6
        
        % Define shape function derivatives
        dNz1 = [ 0.25*(-1+2*xi)*(1-eta)
                -0.25*xi*(-1+xi) ];
        dNz2 = [ 0.25*(1+2*xi)*(1-eta)
                -0.25*xi*(1+xi) ];
        dNz3 = [ 0.25*(1+2*xi)*(1+eta)
                 0.25*xi*(1+xi) ];
        dNz4 = [ 0.25*(-1+2*xi)*(1+eta)
                 0.25*xi*(-1+xi) ];
        dNz5 = [ -xi*(1-eta)   
                -0.5*(1-xi^2) ];
        dNz6 = [ -xi*(1+eta)   
                 0.5*(1-xi^2) ]; 
        % Gradient matrix   
        dNz = [dNz1 dNz2 dNz3 dNz4 dNz5 dNz6];

    case 8
        
        % Define shape function derivatives
        dNz4 = [ -0.25*(1+eta)*(eta-2*xi)
                 -0.25*(-1+xi)*(2*eta-xi) ]; 

        dNz3 = [  0.25*(1+eta)*(2*xi+eta)
                  0.25*(1+xi)*(xi+2*eta)  ];
  
        dNz2 = [  0.25*(-1+eta)*(eta-2*xi)  
                  0.25*(1+xi)*(2*eta-xi)  ];
         
        dNz1 = [ -0.25*(-1+eta)*(2*xi+eta) 
                 -0.25*(-1+xi)*(xi+2*eta) ]; 
        
        dNz6 = [ -xi*(1+eta)   
                  0.5*(1-xi^2) ]; 
         
        dNz7 = [  0.5*(1-eta^2) 
                 -eta*(1+xi)   ];     
 
        dNz5 = [ -xi*(1-eta)   
                 -0.5*(1-xi^2) ];
         
        dNz8 = [ -0.5*(1-eta^2)
                 -eta*(1-xi)   ];
        % Gradient matrix                 
        dNz = [dNz1 dNz2 dNz3 dNz4 dNz5 dNz6 dNz7 dNz8];
        
        case 9
        
        % Define shape function derivatives
        dNz4 = [ 0.25*(2*xi-1)*(eta^2+eta)
                 0.25*(xi^2-xi)*(2*eta+1) ]; 

        dNz3 = [ 0.25*(2*xi+1)*(eta^2+eta)
                 0.25*(xi^2+xi)*(2*eta+1) ];
  
        dNz2 = [ 0.25*(2*xi+1)*(eta^2-eta)
                 0.25*(xi^2+xi)*(2*eta-1) ];
         
        dNz1 = [ 0.25*(2*xi-1)*(eta^2-eta)
                 0.25*(xi^2-xi)*(2*eta-1) ];
        
        dNz6 = [ -xi*(eta^2+eta)   
                 0.5*(1-xi^2)*(2*eta+1) ]; 
         
        dNz7 = [ 0.5*(2*xi+1)*(1-eta^2)
                 -(xi^2+xi)*eta   ];     
 
        dNz5 = [ -xi*(eta^2-eta)   
                 0.5*(1-xi^2)*(2*eta-1) ];
         
        dNz8 = [ 0.5*(2*xi-1)*(1-eta^2)
                 -(xi^2-xi)*eta   ];
             
        dNz9 = [ -2*xi*(1-eta^2)
                 -2*eta*(1-xi^2)];
             
        % Gradient matrix                 
        dNz = [dNz1 dNz2 dNz3 dNz4 dNz5 dNz6 dNz7 dNz8 dNz9];

end

end