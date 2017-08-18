function PlotStress2D(T,X,S,nip,comp)
%--------------------------------------------------------
% PlotStress3D.m:
%   Plots a 3D color contour displaying the stresses 
%   or strains of the hexahedral elements, held in the 
%   topology matrix T, using the coordinate matrix X.
%
% Syntax:
%   PlotStress3D(T,X,S,comp)
%   PlotStress3D(T,X,S)
%
% Input:
%   T    : System topology array.
%   X    : System nodal coordinate array (deformed).
%   S    : System stress/strain array.
%  comp  : Stress/strain component to plot 
%        (1 = x, 2 = y, 3 = z, 4 = yz, 5 = zx, 6 = xy)
%        (other values = Von Mises stress)
%
% Date:
%   Version 3.0   31.10.14
%--------------------------------------------------------
ne = size(T,2)-1;
nel = size(T,1);
ncomp = size(S,1);

% Initial graphics commands
colormap(jet);
hold on;

% Loop over elements to plot stresses/strains as contour 
for e = 1:nel
    
    % Extract stress component 'comp', if given 
    % Von Mises stress is default
    if nargin == 5 && comp <= ncomp
        Sc = S(comp,:,e);
    elseif nargin < 4 || comp > ncomp
        Sc = VMSolid(S(:,:,e));
    end
    
    % Extrapolate to nodes
    Snodes = ExtraPol(Sc,ne,nip);  

    % Order array
    if ne == 4
        order{1} = [1 2 3 4 1];
    elseif ne == 8
        order{1} = [1 5 2 7 3 6 4 8 1];
    elseif ne == 9
        order{1} = [1 5 9 8 1];
        order{2} = [5 2 7 9 5];
        order{3} = [9 7 3 6 9];
        order{4} = [8 9 6 4 8];
        order{5} = [1 5 2 7 3 6 4 8 1];
    end
    
    % face
    x1(1:length(order{1}),:) = X(T(e,order{1}),:);
    x1(length(order{1})+1,:) = X(T(e,1),:);
    col1(1:length(order{1})) = Snodes(order{1});
    col1(length(order{1})+1) = Snodes(1);
    
    % Plot the six element stress/strain face contours
    P1=patch(x1(:,1),x1(:,2),col1,'EdgeColor','none');
    
    if ne>=9
        %  face
    x2(1:length(order{2}),:) = X(T(e,order{2}),:);
    x2(length(order{2})+1,:) = X(T(e,5),:);
    col2(1:length(order{2})) = Snodes(order{2});
    col2(length(order{2})+1) = Snodes(5);       
    
    % Plot the six element stress/strain face contours
    P2=patch(x2(:,1),x2(:,2),col2,'EdgeColor','none');
        %  face
    x3(1:length(order{3}),:) = X(T(e,order{3}),:);
    x3(length(order{3})+1,:) = X(T(e,9),:);
    col3(1:length(order{3})) = Snodes(order{3});
    col3(length(order{3})+1) = Snodes(9);       
    
    % Plot the six element stress/strain face contours
    P3=patch(x3(:,1),x3(:,2),col3,'EdgeColor','none');
        %  face
    x4(1:length(order{4}),:) = X(T(e,order{4}),:);
    x4(length(order{4})+1,:) = X(T(e,8),:);
    col4(1:length(order{4})) = Snodes(order{4});
    col4(length(order{4})+1) = Snodes(8);       
    
    % Plot the six element stress/strain face contours
    P4=patch(x4(:,1),x4(:,2),col4,'EdgeColor','none');
    end
    if ne<=8
    plot(X(T(e,order{1}),1),X(T(e,order{1}),2),'k-');
    else
    plot(X(T(e,order{5}),1),X(T(e,order{5}),2),'k-');
    end
end

axis('equal');    % Set axis dimensions
axis('off');      % Turn off axes

axis(1.1*[min(X(:,1)),max(X(:,1)),min(X(:,2)),max(X(:,2))])



 
function Snodes = ExtraPol(Sc,ne,nip)
%% ------------------------------------------------------
% File: ExtraPol.m
%   Extrapolates stress/strain values at Gauss
%   point to nodes.
%
% Syntax:
%   Snodes = ExtraPol(Sc,ne)
%
%  Input:
%    Sc : Element stress/strain components.
%    ne : Number of element nodes.
%
% Output:
% Snodes : Element stresses/strains at nodes.
%
% Date:
%   Version 1.0    01.09.10
%--------------------------------------------------------
% Get Gauss points
[r w] = Gauss(nip);

% Get normalized node coordinates
Xp = [ -1   1   1  -1  
       -1  -1   1   1 ];

 % For 8 node elements
if ne == 8
    
 Xp = [ -1   1   1  -1   0   0   1  -1
        -1  -1   1   1  -1   1   0   0];
    
elseif ne >= 9
    
 Xp = [ -1   1   1  -1   0   0   1  -1  0
        -1  -1   1   1  -1   1   0   0  0];         
         
end

% Calculate nodal values
Snodes = zeros(ne,1);
for n = 1:ne
    L(:,1) = Lagrange(Xp(1,n),r);
    L(:,2) = Lagrange(Xp(2,n),r);
    
    for i= 1:nip
         for j = 1:nip

                Snodes(n,1) = Snodes(n,1) + L(i,1)*L(j,2)*Sc(1,(i-1)*nip+j);

        end
    end
end



function L = Lagrange(zeta,r)
%% ------------------------------------------------------
% File: Lagrange.m
%   Determines extrapolation coefficients. 
%
% Syntax:
%   L = Lagrange(zeta,r)
%
% Input:
% zeta   : Point with known function value.
%   r    : Number of element nodes.
% Output:
%   L    : Extrapolation coefficients.
%
% Date:
%   Version 1.0    01.09.10
% --------------------------------------------------------
nip = length(r);
L = ones(nip,1);
for j = 1:nip
    for i = 1:nip
        if i ~= j
            L(j) = L(j)*(zeta - r(i))/(r(j) - r(i));
        end
    end
end

