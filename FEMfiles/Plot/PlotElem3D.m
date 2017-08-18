function PlotElem3D(T,X,strc,strl,nodenum)
%--------------------------------------------------------
% PlotElem3D.m:
%   Plots solid elements in topology matrix T with
%   coordinate matrix X. Uses linear line segment
%   between all nodes.
%
% Syntax:
%   PlotElem3D(T,X,strc,strl,nodenum)
%
% Input:
%   T       :  Element topology matrix.
%   X       :  Node coordinate matrix.
%   strc    :  String- or self defining color ('y',[0.5 0.6 0.1],...).
%   strl    :  String defining linetype ('-',':',...).
%   nodenum :  If nonum > 1 -> node numbers
%              will be plotted
% Date:
%   Version 3.0   31.10.14
%--------------------------------------------------------

if nodenum >= 1
    % Plot nodes
    plot3(X(:,1),X(:,2),X(:,3),[strc 'o'],'MarkerFaceColor',strc,'MarkerSize',4)%1-1.5-2-3-4
    hold on
    
    % determine offset for node numbers
    offset = min([max(abs(X(:,1))),max(abs(X(:,2))),max(abs(X(:,3)))])*1e-1;
        
    % Plot node numbers
    if nodenum > 1
        for i=1:size(X,1)
            text(X(i,1)+offset,X(i,2)+offset,X(i,3)+offset,int2str(i),'FontSize',14)
        end
    end
end

% Loop over 3D elements to plot connecting lines
for j = 1:size(T,1)
    if (T(j,end)==1 || T(j,end)==2)
    
    % Define number of nodes pr element
    nnodes = size(T,2)-1;

    % % Define order for H8 elements
    if nnodes == 8
        order = [1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4 ];
    elseif nnodes == 12
        order = [1 9 2 3 10 4 1 5 12 6 7 11 8 5 8 4 10 3 7 6 2];
    elseif nnodes == 16
        order = [1 9 2 13 3 10 4 14 1 5 12 6 16 7 11 8 15 5 15 8 4 10 3 7 16 6 2];
    elseif nnodes == 20
        order = [1 9 2 13 3 10 4 14 1 17 5 12 6 16 7 11 8 15 5 15 8 20 4 10 3 19 7 16 6 18 2];
    end

    % Plot element by calling function 'plot'
    plot3(X(T(j,order),1),X(T(j,order),2),X(T(j,order),3),strl,'color',strc),hold on
    end
end

axis('equal')    % Set axis dimensions
axis('off')      % Turn off axes
axis(1.05*[min(X(:,1)),max(X(:,1)),min(X(:,2)),max(X(:,2)),min(X(:,3)),max(X(:,3))]);%max(X(:,3))max(X(1:7806,3))])