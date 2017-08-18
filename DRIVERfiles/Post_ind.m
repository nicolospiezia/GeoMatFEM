%--------------------------------------------------------
% File: Postprocessore_node
%
% Driver for postprocessor of node values for coupled problems.
%
% Date: 31/10/2014
%   Version 3.0   
% Created by: Nicolò Spiezia
% Updated by: Nico De Marchi
%--------------------------------------------------------
clear F

 scl = 1;
UU = reshape(U(:,end),size(X,2),size(X,1))'*scl;
 % Add folder to Matlab search path
addpath(genpath('..\FEMfiles'));%

fig=figure;

    if dof==2
   
        PlotElem2D(T,X,'k','-',1);
        hold on
        %PlotNode2D(T,X,DOF,solution,1,1)
    else
        
        % Plot displaced configuration X+u
    
        PlotElem3D(T,X,'k',':',0)%,drawnow
        hold on
        PlotElem3D(T,X+UU,'k','-',0);
        hold on
        %PlotNode3D(T,X,DOF,solution,1,1)
        view(180,0)%(-25,25)
    end
    %title('Elements and nodes mesh')
    % Set min and max contour
    
    axis('equal');    % Set axis dimensions
    axis('off');      % Turn off axes
    axis([min(X(:,1)),max(X(:,1)),min(X(:,2)),max(X(:,2))])

    % Print output
    
    filename = ['mesh', '_', 'fsCC0'];
    
    set(gca, 'FontSize', 12); % Font size

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 15 10]);  
    set(gcf, 'PaperSize', [15 10]);
    
    print('-f1', '-r600', '-dpdf', '-painters', filename);%-dpdf-dpng
    
    close all
    clear title


