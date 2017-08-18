function PlotDisp3D(T,X,U,comp)
%--------------------------------------------------------
% PlotDisp3D.m:
%   Plots a 3D color contour displaying the stresses 
%   or strains of the hexahedral elements, held in the 
%   topology matrix T, using the coordinate matrix X.
%
% Syntax:
%   PlotDisp3D(T,X,U,comp)
%   PlotDisp3D(T,X,U)
%
% Input:
%   T    : System topology array.
%   X    : System nodal coordinate array (deformed).
%   U    : Nodal displacement vector [u1 u2 u3]
%  comp  : Displacement component
%
% Date:
%   Version 3.0   31.10.14
%--------------------------------------------------------

% Determine number of nodes per element
ne = size(T,2)-1;
nno = size(X,1);

% % Initial graphics commands
colormap(jet); 

hold on;

% Loop over elements to plot stresses/strains as contour 
for i = 1:size(T,1)    
    if (T(i,end)==1 || T(i,end)==2)  
    for j = 1:ne
        if nargin == 3
            val(i,j) = 0;
        else
            
            val(i,j) = U(T(i,j),comp);
        end
    end
    
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
    x1(1:length(order{1}),:) = X(T(i,order{1}),:);
    x1(length(order{1})+1,:) = X(T(i,1),:);
    col1(1:length(order{1})) = val(i,order{1});
    col1(length(order{1})+1) = val(i,1);
    
   % Side faces
    x2(1:length(order{2}),:) = X(T(i,order{2}),:);
    x2(length(order{2})+1,:) = X(T(i,1),:);
    col2(1:length(order{2})) = val(i,order{2});
    col2(length(order{2})+1) = val(i,1);
    x3(1:length(order{3}),:) = X(T(i,order{3}),:);
    x3(length(order{3})+1,:) = X(T(i,4),:);
    col3(1:length(order{3})) = val(i,order{3});
    col3(length(order{3})+1) = val(i,4);
    
    % Top face
    x4(1:length(order{4}),:) = X(T(i,order{4}),:);
    x4(length(order{4})+1,:) = X(T(i,5),:);
    col4(1:length(order{4})) = val(i,order{4});
    col4(length(order{4})+1) = val(i,5);
    
    % Back face
    x5(1:length(order{5}),:) = X(T(i,order{5}),:);
    x5(length(order{5})+1,:) = X(T(i,1),:);
    col5(1:length(order{5})) = val(i,order{5});
    col5(length(order{5})+1) = val(i,1);

    % Front face
    x6(1:length(order{6}),:) = X(T(i,order{6}),:);
    x6(length(order{6})+1,:) = X(T(i,2),:);
    col6(1:length(order{6})) = val(i,order{6});
    col6(length(order{6})+1) = val(i,2);
        
    % Plot the six element stress/strain face contours
    P1=patch(x1(:,1),x1(:,2),x1(:,3),col1);
    P2=patch(x2(:,1),x2(:,2),x2(:,3),col2);
    P3=patch(x3(:,1),x3(:,2),x3(:,3),col3);
    P4=patch(x4(:,1),x4(:,2),x4(:,3),col4);
    P5=patch(x5(:,1),x5(:,2),x5(:,3),col5);
    P6=patch(x6(:,1),x6(:,2),x6(:,3),col6);
    end
end

axis('equal');    % Set axis dimensions
axis('off');      % Turn off axes
axis(1.1*[min(X(:,1)),max(X(:,1)),min(X(:,2)),max(X(:,2)),min(X(:,3)),max(X(:,3))]);%max(X(:,3) max(X(1:7806,3)

