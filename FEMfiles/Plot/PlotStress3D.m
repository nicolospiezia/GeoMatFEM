function PlotStress3D(T,X,S,nip,comp)
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
    %if (T(e,end)==1 || T(e,end)==2)
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
        order{1} = [1 2 3 4];
        order{2} = [1 2 6 5];
        order{3} = [4 3 7 8];
        order{4} = [5 6 7 8];
        order{5} = [1 4 8 5];
        order{6} = [2 3 7 6];
        
    if ne >= 12
        order{1} = [1 9 2 3 10 4];
        order{2} = [1 9 2 6 12 5];
        order{3} = [4 10 3 7 11 8];
        order{4} = [5 12 6 7 11 8];
    elseif ne >= 16
        order{1} = [1 9 2 13 3 10 4 14];
        order{2} = [1 9 2 6 12 5];
        order{3} = [4 10 3 7 11 8];
        order{4} = [5 12 6 16 7 11 8 15];
        order{5} = [1 14 4 8 15 5];
        order{6} = [2 13 3 7 16 6];
    elseif ne >= 20
        order{1} = [1 9 2 13 3 10 4 14];
        order{2} = [1 9 2 18 6 12 5 17];
        order{3} = [4 10 3 19 7 11 8 20];
        order{4} = [5 12 6 16 7 11 8 15];
        order{5} = [1 14 4 20 8 15 5 17];
        order{6} = [2 13 3 19 7 16 6 18];
    end
        
    
    % Bottom face
    x1(1:length(order{1}),:) = X(T(e,order{1}),:);
    x1(length(order{1})+1,:) = X(T(e,1),:);
    col1(1:length(order{1})) = Snodes(order{1});
    col1(length(order{1})+1) = Snodes(1);
    
   % Side faces
    x2(1:length(order{2}),:) = X(T(e,order{2}),:);
    x2(length(order{2})+1,:) = X(T(e,1),:);
    col2(1:length(order{2})) = Snodes(order{2});
    col2(length(order{2})+1) = Snodes(1);
    x3(1:length(order{3}),:) = X(T(e,order{3}),:);
    x3(length(order{3})+1,:) = X(T(e,4),:);
    col3(1:length(order{3})) = Snodes(order{3});
    col3(length(order{3})+1) = Snodes(4);
    
    % Top face
    x4(1:length(order{4}),:) = X(T(e,order{4}),:);
    x4(length(order{4})+1,:) = X(T(e,5),:);
    col4(1:length(order{4})) = Snodes(order{4});
    col4(length(order{4})+1) = Snodes(5);
    
    % Back face
    x5(1:length(order{5}),:) = X(T(e,order{5}),:);
    x5(length(order{5})+1,:) = X(T(e,1),:);
    col5(1:length(order{5})) = Snodes(order{5});
    col5(length(order{5})+1) = Snodes(1);

    % Front face
    x6(1:length(order{6}),:) = X(T(e,order{6}),:);
    x6(length(order{6})+1,:) = X(T(e,2),:);
    col6(1:length(order{6})) = Snodes(order{6});
    col6(length(order{6})+1) = Snodes(2);
        
    % Plot the six element stress/strain face contours
    P1=patch(x1(:,1),x1(:,2),x1(:,3),col1);
    P2=patch(x2(:,1),x2(:,2),x2(:,3),col2);
    P3=patch(x3(:,1),x3(:,2),x3(:,3),col3);
    P4=patch(x4(:,1),x4(:,2),x4(:,3),col4);
    P5=patch(x5(:,1),x5(:,2),x5(:,3),col5);
    P6=patch(x6(:,1),x6(:,2),x6(:,3),col6);
    %end
end

axis('equal');    % Set axis dimensions
axis('off');      % Turn off axes

axis(1.1*[min(X(:,1)),max(X(:,1)),min(X(:,2)),max(X(:,2)),min(X(:,3)),max(X(:,3))])



 
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
Xp = [ -1   1   1  -1  -1   1   1  -1
       -1  -1   1   1  -1  -1   1   1
       -1  -1  -1  -1   1   1   1   1 ];

% For 12 node elements
if ne >= 12
    
    Xp12 = [  0   0   0   0
             -1   1   1  -1
             -1  -1   1   1];
    
    Xp = [Xp Xp12];
end

% For 16 node elements
if ne >= 16
    
    Xp16 = [  1  -1  -1   1
              0   0   0   0
             -1  -1   1   1];
    
    Xp = [Xp Xp16];
end

% For 20 node elements
if ne >= 20
    
    Xp20 = [ -1   1   1  -1
             -1  -1   1   1
              0   0   0   0];
    
    Xp = [Xp Xp20];
end

% Calculate nodal values
Snodes = zeros(ne,1);
for n = 1:ne
    L(:,1) = Lagrange(Xp(1,n),r);
    L(:,2) = Lagrange(Xp(2,n),r);
    L(:,3) = Lagrange(Xp(3,n),r);
    
    for i= 1:nip
         for j = 1:nip
            for k = 1:nip
                Snodes(n,1) = Snodes(n,1) + L(i,1)*L(j,2)*L(k,3)*Sc(1,(i-1)*nip^2+(j-1)*nip+k);
            end
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

