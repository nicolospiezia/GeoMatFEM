%--------------------------------------------------------
% File: Postprocessore
%
% Driver for 3D postprocessor of node values for coupled problems.
%
% Date: 31/10/2014
%   Version 3.0   
% Created by: Nicolò Spiezia
% Updated by: Nico De Marchi
%-------------------------------------------------------
% Time/Displacement curve under silos
%[~,plotdof] = max(abs(ft));
plotdof=782*3;%1329*3;%174*3;%
%-------------------------------------------------------
% Define load vector
loadV=[ 0 incrf];


fig=figure;clf
set(fig,'WindowStyle','docked')
set(gca, 'FontSize', 12); % Font size
set(gca, 'Box', 'on');
hold on
plot(-Us(plotdof,:),loadV,'k-s');%,'LineWidth',1,'MarkerSize',5) 
plot(-U(plotdof,:),loadV,'k-o');%,'LineWidth',1,'MarkerSize',5)
ylabel ('Load -q [MPa]');
xlabel ('Settlement -w [mm]');
%title(['Stress \sigma_{' num2str(1) '} / Height curve']);
legend('Small Strain','Finite Strain','Location','SouthEast')
%xlim([0,0.4])
%--------------------------------------------------------------------------
filename = ['image', '_', 'Pw_CC'];
    
    set(gca, 'FontSize', 12); % Font size

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0.2 0.2 15 10]);  
    set(gcf, 'PaperSize', [15 10]);
    
    print('-f1', '-r600', '-dpdf', '-painters', filename);

    
    close all
    clear title
%--------------------------------------------------------------------------

disp(' DONE!');