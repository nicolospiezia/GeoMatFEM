%--------------------------------------------------------------------------
% File: Postproc.m
%
% Driver for graphical rappresentation of the results 
%
% Date: 31/10/2014
%   Version 3.0   
% Created by: Nicolò Spiezia
% Updated by: Nico De Marchi
%--------------------------------------------------------------------------
% Plot parameters
plotdir=3;
plotcomp=3;
plotdof=2346;%522;%3987;
%[~,plotdof] = max(abs(f));
%plotdof=46773;%plotdof=47928;%destro 30048
plotb=3; % Plot brec variable
plotPQ=1; % Plot PLrec variable
nstep=1; % Time step convergence:
%--------------------------------------------------------------------------
% Add folder to Matlab search path
addpath(genpath('..\FEMfiles'));%path(path,'..\FEMfiles');%path(path,'FEMfiles');
%
% Postprocessing:
%
% Displacement scale factor
scl = 1;
UU = reshape(U(:,end),size(X,2),size(X,1))';
%
if dof==2
    % Plot elements and nodes of topology T
    fig1=figure(1);clf
    set(fig1,'WindowStyle','docked')
    PlotElem2D(T,X,'k','-',1),drawnow
    title('Plot elements and nodes of topology T');
    %
    % Plot displaced configuration X+u
    fig2=figure(2);clf
    set(fig2,'WindowStyle','docked')
    PlotElem2D(T,X,'k',':',0),drawnow
    hold on
    PlotElem2D(T,X+scl*UU,'b','-',1);
    title('Plot displaced configuration X+u');
    %
    % Plot nodal displacements
    fig3=figure(3);clf
    set(fig3,'WindowStyle','docked')
    PlotDisp2D(T,X+scl*UU,UU,plotdir)
    title(['Plot nodal displacements dir=U_{' num2str(plotdir) '}']);
    colorbar
    %
    % Plot stress
    fig4=figure(4);clf
    set(fig4,'WindowStyle','docked')
    PlotStress2D(T,X+scl*UU,S,nip,plotcomp)
    title(['Plot stress \sigma_{' num2str(plotcomp) '}']);
    colorbar
    %
    % Plot Left Cauchy-Green or Elastic strain componets
    fig5=figure(5);clf
    set(fig5,'WindowStyle','docked')
    PlotStress2D(T,X+scl*UU,b,nip,plotb)
    if nonlin==0
        title(['Plot Elastic strain: \epsilon_{' num2str(plotb) '}']);
    else
        title(['Plot Plastic variable: b_{' num2str(plotb) '}']);
    end
    colorbar
    %
    % Plot Plastic var
    fig6=figure(6);clf
    set(fig6,'WindowStyle','docked')
    PlotStress2D(T,X+scl*UU,PL,nip,plotPL)
    title(['Plot Plastic variable: PL_{' num2str(plotPL) '}']);
    colorbar
    %
    % Plot Tension invariant
    fig7=figure(7);clf
    set(fig7,'WindowStyle','docked')
    PlotStress2D(T,X+scl*UU,PQ,nip,plotPQ)
    if plotPQ ==1
        title('Plot tension invariant: P');
    else
        title('Plot tension invariant: Q');
    end
    colorbar
    %
else
    % Plot elements and nodes of topology T
    fig1=figure(1);clf
    set(fig1,'WindowStyle','docked')
    PlotElem3D(T,X,'k','-',0),drawnow
    title('Plot elements and nodes of topology T');
    view(-45,20)
    %
    % Plot displaced configuration X+u
    fig2=figure(2);clf
    set(fig2,'WindowStyle','docked')
    PlotElem3D(T,X,'k',':',0),drawnow
    hold on
    PlotElem3D(T,X+scl*UU,'b','-',0);
    title('Plot displaced configuration X+u');
    view(-45,20)
    %
    % Plot nodal displacements
    fig3=figure(3);clf
    set(fig3,'WindowStyle','docked')
    PlotDisp3D(T,X+scl*UU,UU,plotdir)
    title(['Plot nodal displacements dir=U_{' num2str(plotdir) '}']);
    colorbar
    view(-45,20)
    %
    % Plot stress
    fig4=figure(4);clf
    set(fig4,'WindowStyle','docked')
    PlotStress3D(T,X+scl*UU,S,nip,plotcomp)
    title(['Plot stress \sigma_{' num2str(plotcomp) '}']);
    colorbar
    view(-45,20)
    %
    % Plot Left Cauchy-Green or Elastic strain componets
    fig5=figure(5);clf
    set(fig5,'WindowStyle','docked')
    PlotStress3D(T,X+scl*UU,b,nip,plotb)
    if nonlin==0
        title(['Plot Elastic strain: \epsilon_{' num2str(plotb) '}']);
    else
        title(['Plot Plastic variable: b_{' num2str(plotb) '}']);
    end
    colorbar
    view(-45,20)
    %
    % Plot Plastic var
    fig6=figure(6);clf
    set(fig6,'WindowStyle','docked')
    PlotStress3D(T,X+scl*UU,PL,nip,3)
    title(['Plot Plastic variable: PL_{' num2str(3) '}']);
    colorbar
    view(-45,20)%view(0,0)%
    %
    % Plot Plastic var
    fig7=figure(7);clf
    set(fig7,'WindowStyle','docked')
    PlotStress3D(T,X+scl*UU,PL,nip,1)
    title(['Plot Plastic variable: PL_{' num2str(1) '}']);
    colorbar
    view(-45,20)
    %
    %
    % Plot Plastic var
    fig8=figure(8);clf
    set(fig8,'WindowStyle','docked')
    PlotStress3D(T,X+scl*UU,PL,nip,2)
    title(['Plot Plastic variable: PL_{' num2str(2) '}']);
    colorbar
    view(-45,20)
    %
    % Plot Tension invariant
    fig9=figure(9);clf
    set(fig9,'WindowStyle','docked')
    PlotStress3D(T,X+scl*UU,PQ,nip,plotPQ)
    if plotPQ ==1
        title('Plot tension invariant: P');
    else
        title('Plot tension invariant: Q');
    end
    colorbar
    view(-45,20)
    %
end
%
% Load/displacement curve
fig10=figure(10);clf
set(fig10,'WindowStyle','docked')
axis ([0 abs(U(plotdof,end)) 0 abs(F(plotdof,end))]);
%axis ([0 max(abs(U(plotdof,:))) 0 max(abs(F(plotdof,:)))]);
plot(abs(U(plotdof,1:n+1)),abs(F(plotdof,1:n+1)),'r-*',  ...
     abs(U(plotdof,n+1)),abs((F(plotdof,n+1))),'ko','MarkerSize',5);
hold on
 xlabel ('u [mm]');
ylabel ('F [N]');
title (sprintf('Load/Displacement curve of dof=%d', plotdof));
%
ord=[];
asci=[];
count=0;
for i=1:imax
    if RESrec(i,nstep)~= 0;
        count=count+1;
        ord(count)=RESrec(i,nstep);
    else
    end
end
asci=(1:1:count);
% Residual/iter
fig11=figure(11);clf
set(fig11,'WindowStyle','docked')
semilogy(asci,ord,'b-*')
%plot(asci,ord,'b-*')
hold on
xlabel ('Iter [n]');
ylabel ('Log-res');
title (sprintf('Iteraction/Reridual curve of load step=%d', nstep));
%
% Visualize sparsity pattern
fig12=figure(12);clf
set(fig12,'WindowStyle','docked')
spy(K)
ncz=nnz(K)*100/numel(K);
xlabel (['Non-zero element {' num2str(ncz) '}%']);
hold on
title('Plot Sparsity pattern');
%
%
elelem2=1;
nstep2=1;
iter2=1;
ord=[];
asci=[];
count=0;
for i=1:imax
    if RMAPrec(i,elelem2,iter2,nstep2)~= 0;
        count=count+1;
        ord(count)=RMAPrec(i,elelem2,iter2,nstep2);
    else
    end
end
asci=(1:1:count);
% Residual/iter
fig13=figure(13);clf
set(fig13,'WindowStyle','docked')
semilogy(asci,ord,'b-*')
%plot(asci,ord,'b-*')
hold on
xlabel ('Iter [n]');
ylabel ('Log-res');
title (sprintf('Iteraction/Reridual curve of load step=%d', nstep2));
%
disp('Plot DONE! :-)');