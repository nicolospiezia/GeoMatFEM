function PlotDisp2D(T,X,U,comp)
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
%nno = size(X,1);

% % Initial graphics commands
colormap(jet); 

hold on;

% Loop over elements to plot stresses/strains as contour 
for i = 1:size(T,1)    
        
    for j = 1:ne
        if nargin == 3
            val(i,j) = 0;
        else
            
            val(i,j) = U(T(i,j),comp);
        end
    end

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
    
    
    %  face
    x1(1:length(order{1}),:) = X(T(i,order{1}),:);
    x1(length(order{1})+1,:) = X(T(i,1),:);
    col1(1:length(order{1})) = val(i,order{1});
    col1(length(order{1})+1) = val(i,1);       
    
    % Plot the six element stress/strain face contours
    P1=patch(x1(:,1),x1(:,2),col1,'EdgeColor','none');
    
    if ne>=9
        %  face
    x2(1:length(order{2}),:) = X(T(i,order{2}),:);
    x2(length(order{2})+1,:) = X(T(i,5),:);
    col2(1:length(order{2})) = val(i,order{2});
    col2(length(order{2})+1) = val(i,5);       
    
    % Plot the six element stress/strain face contours
    P2=patch(x2(:,1),x2(:,2),col2,'EdgeColor','none');
        %  face
    x3(1:length(order{3}),:) = X(T(i,order{3}),:);
    x3(length(order{3})+1,:) = X(T(i,9),:);
    col3(1:length(order{3})) = val(i,order{3});
    col3(length(order{3})+1) = val(i,9);       
    
    % Plot the six element stress/strain face contours
    P3=patch(x3(:,1),x3(:,2),col3,'EdgeColor','none');
        %  face
    x4(1:length(order{4}),:) = X(T(i,order{4}),:);
    x4(length(order{4})+1,:) = X(T(i,8),:);
    col4(1:length(order{4})) = val(i,order{4});
    col4(length(order{4})+1) = val(i,8);       
    
    % Plot the six element stress/strain face contours
    P4=patch(x4(:,1),x4(:,2),col4,'EdgeColor','none');
    end
    if ne<=8
    plot(X(T(i,order{1}),1),X(T(i,order{1}),2),'k-');
    else
    plot(X(T(i,order{5}),1),X(T(i,order{5}),2),'k-');
    end
end

axis('equal');    % Set axis dimensions
axis('off');      % Turn off axes
axis(1.1*[min(X(:,1)),max(X(:,1)),min(X(:,2)),max(X(:,2))])

