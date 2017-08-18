%--------------------------------------------------------------------------
nstep=[1 8 15 25];
asci=1:12;
ord1=[];
asci1=[];
count1=0;
for i=1:imax
    if RESrec(i,nstep(1))~= 0;
        count1=count1+1;
        ord1(count1)=RESrec(i,nstep(1));
    else
    end
end
asci1=(1:1:count1);

ord2=[];
asci2=[];
count2=0;
for i=1:imax
    if RESrec(i,nstep(2))~= 0;
        count2=count2+1;
        ord2(count2)=RESrec(i,nstep(2));
    else
    end
end
asci2=(1:1:count2);

ord3=[];
asci3=[];
count3=0;
for i=1:imax
    if RESrec(i,nstep(3))~= 0;
        count3=count3+1;
        ord3(count3)=RESrec(i,nstep(3));
    else
    end
end
asci3=(1:1:count3);

ord4=[];
asci4=[];
count4=0;
for i=1:imax
    if RESrec(i,nstep(4))~= 0;
        count4=count4+1;
        ord4(count4)=RESrec(i,nstep(4));
    else
    end
end
asci4=(1:1:count4);
%plot()

% Residual/iter
ytol=ones(12,1)*epsr;
fig=figure;
%set(fig,'WindowStyle','docked')

%axis('manual');
semilogy(asci1,ord1,'k-s',asci2,ord2,'k-o',asci3,ord3,'k-v',asci4,ord4,'k-*',asci,ytol,'k:')
%semilogy(asci1,ord1,'k-s')
%semilogy(asci2,ord2,'k-o')
%semilogy(asci3,ord3,'k-v')
%semilogy(asci4,ord4,'k-*')
%semilogy(asci4,ytol,'k:')
%

axis([1,10,1e-20,1e5])
xlabel ('Iteration','FontSize', 15);
ylabel ('Normalized residual','FontSize', 15);
%title (sprintf('Iteraction/Reridual curve of time step=%d', nstep));

legend(['n_{' num2str(nstep(1),'%d') '}'],['n_{' num2str(nstep(2),'%d') '}'],['n_{' num2str(nstep(3),'%d') '}'],['n_{' num2str(nstep(4),'%d') '}'],['Toll'],'Location','NorthEast')
filename = ['image', '_', 'residuo_fsCC0'];
    
    set(gca, 'FontSize', 12); % Font size

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0.2 0.2 15 10]);  
    set(gcf, 'PaperSize', [15 10]);
    
    print('-f1', '-r600', '-dpdf', '-painters', filename);

    
    close all
    clear title