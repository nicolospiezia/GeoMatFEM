function N = NSolid(ne,dof,dmx,xi,eta,zeta)
%--------------------------------------------------------
% NSolid:
%   Computes the element shape function matrix at
%   coordinates (xi,eta,zeta) for an elastic 
%   8/12/16/20 nodes solid hexahedral element or 4/8 nodes quadrilateral
%   element in plane strain or 2/3 nodes of linear element.
%
% Syntax:
%   N = NSolid(Xe,xi,eta,zeta)
%
% Input:
%   Xe   : Element nodal coordinate array.
%   xi   : First normalized coordinate.
%   eta  : Second normalized coordinate.
%   zeta : Third normalized coordinate.
%
% Output:
%   N    : Element shape function matrix.
%
% Date: 31/10/2014
%   Version 3.0   
% Updated by: Nico De Marchi
%--------------------------------------------------------

% Number of nodes per element and dof's per node
%ne = size(Xe,1);
%dof = size(Xe,2);

% Unit matrix
I = eye(dmx);
if dof == 1

%  Nodes numbering:
%   1---(3)---2 
    switch ne 
        case 2
            % Define shape function derivatives
            N1 = (1-xi)/2; 
            N2 = (1+xi)/2;
            % Gradient matrix     
            N = [N1*I N2*I];
        case 3    
            % Define shape function derivatives
            N1 = xi*(xi-1)/2; 
            N2 = xi*(xi+1)/2;
            N3 = 1-xi^2;
            % Gradient matrix   
            N = [N1*I N2*I N3*I];
    end
    
elseif dof == 2  
    
%  Nodes numbering:
%
%   4---(6)---3 
%   |         | 
%  (8)  (9)  (7)
%   |         | 
%   1---(5)---2 
%    
    switch ne
        case 4
            % Vertex nodes
            N4 = 0.25*(1-xi+eta-xi*eta);
            N3 = 0.25*(1+xi+eta+xi*eta);
            N2 = 0.25*(1+xi-eta-xi*eta);
            N1 = 0.25*(1-xi-eta+xi*eta);
    
            % Define shape function matrix at Gauss points
            N = [ N1*I  N2*I  N3*I  N4*I];
            
        case 6
            % Vertex nodes
            N1 = 0.25*xi*(-1+xi)*(1-eta);
            N2 = 0.25*xi*(1+xi)*(1-eta);
            N3 = 0.25*xi*(1+xi)*(1+eta);
            N4 = 0.25*xi*(-1+xi)*(1+eta); 
            N5 = 0.5*(1-xi^2)*(1-eta);
            N6 = 0.5*(1-xi^2)*(1+eta);
    
            % Define shape function matrix at Gauss points 
            N = [ N1*I N2*I N3*I N4*I N5*I N6*I ];
        
        case 8
            % Vertex nodes
            N4 = 0.25*(1-xi)*(1+eta)*(-xi+eta-1);
            N3 = 0.25*(1+xi)*(1+eta)*(+xi+eta-1); 
            N2 = 0.25*(1+xi)*(1-eta)*(+xi-eta-1);        
            N1 = 0.25*(1-xi)*(1-eta)*(-xi-eta-1);       

            % Central nodes
            N6 = 0.5*(1-xi^2)*(1+eta);
            N7 = 0.5*(1+xi)*(1-eta^2);
            N5 = 0.5*(1-xi^2)*(1-eta);
            N8 = 0.5*(1-xi)*(1-eta^2);
    
            % Define shape function matrix at Gauss points
            N = [N1*I N2*I N3*I N4*I N5*I N6*I N7*I N8*I];
            
        case 9
            % Vertex nodes
            N4 = 0.25*(xi^2-xi)*(eta^2+eta);
            N3 = 0.25*(xi^2+xi)*(eta^2+eta); 
            N2 = 0.25*(xi^2+xi)*(eta^2-eta);        
            N1 = 0.25*(xi^2-xi)*(eta^2-eta);       

            % Central nodes
            N6 = 0.5*(1-xi^2)*(eta^2+eta);
            N7 = 0.5*(xi^2+xi)*(1-eta^2);
            N5 = 0.5*(1-xi^2)*(eta^2-eta);
            N8 = 0.5*(xi^2-xi)*(1-eta^2);
            
            N9 = (1-xi^2)*(1-eta^2);
    
            % Define shape function matrix at Gauss points
            N = [N1*I N2*I N3*I N4*I N5*I N6*I N7*I N8*I N9*I];
    end
else
    
    % Define shape functions for corner nodes
    N1 = -(xi-1)*(eta-1)*(zeta-1);  N2 =  (xi+1)*(eta-1)*(zeta-1);
    N3 = -(xi+1)*(eta+1)*(zeta-1);  N4 =  (xi-1)*(eta+1)*(zeta-1);
    N5 =  (xi-1)*(eta-1)*(zeta+1);  N6 = -(xi+1)*(eta-1)*(zeta+1);
    N7 =  (xi+1)*(eta+1)*(zeta+1);  N8 = -(xi-1)*(eta+1)*(zeta+1);

    % Define shape function matrix at Gauss points
    N = [N1*I N2*I N3*I N4*I N5*I N6*I N7*I N8*I ]/8;

    if ne >= 12
    
        % Quadratic terms for midside nodes (xi-dir)
        N9  = -1/4*(-1+xi)*(1+xi)*(-1+eta)*(-1+zeta);
        N10 =  1/4*(-1+xi)*(1+xi)*(1+eta)*(-1+zeta);
        N11 = -1/4*(-1+xi)*(1+xi)*(1+eta)*(1+zeta);
        N12 =  1/4*(-1+xi)*(1+xi)*(-1+eta)*(1+zeta);
    
        % Modify corner nodes
        N(:,1:3) = N(:,1:3) - 1/2*N9*I;             % N1
        N(:,4:6) = N(:,4:6) - 1/2*N9*I;             % N2
        N(:,7:9) = N(:,7:9) - 1/2*N10*I;            % N3
        N(:,10:12) = N(:,10:12) - 1/2*N10*I;        % N4
        N(:,13:15) = N(:,13:15) - 1/2*N12*I;        % N5
        N(:,16:18) = N(:,16:18) - 1/2*N12*I;        % N6
        N(:,19:21) = N(:,19:21) - 1/2*N11*I;        % N7
        N(:,22:24) = N(:,22:24) - 1/2*N11*I;        % N8
    
        N12 = [N9*I N10*I N11*I N12*I];
    
        % Expand interpolation matrix
        N = [N N12];
    end
    if ne >= 16
    
        % Quadratic terms for midside nodes (eta-dir)
        N13 =  1/4*(-1+eta)*(1+eta)*(1+xi)*(-1+zeta);
        N14 = -1/4*(-1+eta)*(1+eta)*(-1+xi)*(-1+zeta);
        N15 =  1/4*(-1+eta)*(1+eta)*(-1+xi)*(1+zeta);
        N16 = -1/4*(-1+eta)*(1+eta)*(1+xi)*(1+zeta);
    
        % Modify corner nodes
        N(:,1:3) = N(:,1:3) - 1/2*N14*I;            % N1
        N(:,4:6) = N(:,4:6) - 1/2*N13*I;            % N2
        N(:,7:9) = N(:,7:9) - 1/2*N13*I;            % N3
        N(:,10:12) = N(:,10:12) - 1/2*N14*I;        % N4
        N(:,13:15) = N(:,13:15) - 1/2*N15*I;        % N5
        N(:,16:18) = N(:,16:18) - 1/2*N16*I;        % N6
        N(:,19:21) = N(:,19:21) - 1/2*N16*I;        % N7
        N(:,22:24) = N(:,22:24) - 1/2*N15*I;        % N8
    
        N16 = [N13*I N14*I N15*I N16*I];
    
        % Expand interpolation matrix
        N = [N N16];
    end

    if ne >= 20
        % Quadratic terms for midside nodes (zeta-dir)
        N17 = -1/4*(-1+zeta)*(1+zeta)*(-1+xi)*(-1+eta);
        N18 =  1/4*(-1+zeta)*(1+zeta)*(1+xi)*(-1+eta);
        N19 = -1/4*(-1+zeta)*(1+zeta)*(1+xi)*(1+eta);
        N20 =  1/4*(-1+zeta)*(1+zeta)*(-1+xi)*(1+eta);
    
        % Modify corner nodes
        N(:,1:3) = N(:,1:3) - 1/2*N17*I;            % N1
        N(:,4:6) = N(:,4:6) - 1/2*N18*I;            % N2
        N(:,7:9) = N(:,7:9) - 1/2*N19*I;            % N3
        N(:,10:12) = N(:,10:12) - 1/2*N20*I;        % N4
        N(:,13:15) = N(:,13:15) - 1/2*N17*I;        % N5
        N(:,16:18) = N(:,16:18) - 1/2*N18*I;        % N6
        N(:,19:21) = N(:,19:21) - 1/2*N19*I;        % N7
        N(:,22:24) = N(:,22:24) - 1/2*N20*I;        % N8
        
        N20 = [N17*I N18*I N19*I N20*I];
    
        % Expand interpolation matrix
        N = [N N20];
    end
end
end

