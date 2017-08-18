%--------------------------------------------------------------------------
elem=73;
nt=15;
%--------------------------------------------------------------------------
nstep=[2 3 4 5 6 7];
asci=1:15;
ord1=[];
asci1=[];
count1=0;
for i=1:imax
    if RMAPrec(i,elem,nstep(1),nt)~= 0;
        count1=count1+1;
        ord1(count1)=RMAPrec(i,elem,nstep(1),nt);
    else
    end
end
if count1==0
    ord1=[1];asci1=[1];
else
asci1=(1:1:count1);
end

ord2=[];
asci2=[];
count2=0;
for i=1:imax
    if RMAPrec(i,elem,nstep(2),nt)~= 0;
        count2=count2+1;
        ord2(count2)=RMAPrec(i,elem,nstep(2),nt);
    else
    end
end
if count2==0
    ord2=[1];asci2=[1];
else
asci2=(1:1:count2);
end

ord3=[];
asci3=[];
count3=0;
for i=1:imax
    if RMAPrec(i,elem,nstep(3),nt)~= 0;
        count3=count3+1;
        ord3(count3)=RMAPrec(i,elem,nstep(3),nt);
    else
    end
end
if count3==0
    ord3=[1];asci3=[1];
else
asci3=(1:1:count3);
end

ord4=[];
asci4=[];
count4=0;
for i=1:imax
    if RMAPrec(i,elem,nstep(4),nt)~= 0;
        count4=count4+1;
        ord4(count4)=RMAPrec(i,elem,nstep(4),nt);
    else
    end
end
if count4==0
    ord4=[1];asci4=[1];
else
asci4=(1:1:count4);
end

ord5=[];
asci5=[];
count5=0;
for i=1:imax
    if RMAPrec(i,elem,nstep(5),nt)~= 0;
        count5=count5+1;
        ord5(count5)=RMAPrec(i,elem,nstep(5),nt);
    else
    end
end
if count5==0
    ord5=[1];asci5=[1];
else
asci5=(1:1:count5);
end

ord6=[];
asci6=[];
count6=0;
for i=1:imax
    if RMAPrec(i,elem,nstep(6),nt)~= 0;
        count6=count6+1;
        ord6(count6)=RMAPrec(i,elem,nstep(6),nt);
    else
    end
end
if count6==0
    ord6=[1];asci6=[1];
else
asci6=(1:1:count6);
end

% Residual/iter
ytol=ones(size(asci,2),1)*1e-6;%<--------- tolleranza locale rm
fig=figure;
%set(fig,'WindowStyle','docked')

%axis('manual');
semilogy(asci1,ord1,'k-s',asci2,ord2,'k-o',asci3,ord3,'k-v',asci4,ord4,'k-*',asci5,ord5,'k-d',asci5,ord5,'k-^',asci,ytol,'k:')%
%semilogy(asci1,ord1,'k-s')
%semilogy(asci2,ord2,'k-o')
%semilogy(asci3,ord3,'k-v')
%semilogy(asci4,ord4,'k-*')
%semilogy(asci4,ytol,'k:')
%

axis([1,10,1e-20,1e5])
xlabel ('Local iteration','FontSize', 15);
ylabel ('Normalized residual','FontSize', 15);
%title (sprintf('Iteraction/Reridual curve of time step=%d', nstep));

legend(['i_{' num2str(nstep(1),'%d') '}'],['i_{' num2str(nstep(2),'%d') '}'],['i_{' num2str(nstep(3),'%d') '}'],['i_{' num2str(nstep(4),'%d') '}'],['i_{' num2str(nstep(5),'%d') '}'],['i_{' num2str(nstep(6),'%d') '}'],['Toll'],'Location','NorthEast')%['i_{' num2str(nstep(4),'%d') '}'],
filename = ['image', '_', 'rloc_el73_n15_fsCC1'];
    
    set(gca, 'FontSize', 12); % Font size

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0.2 0.2 15 10]);  
    set(gcf, 'PaperSize', [15 10]);
    
    print('-f1', '-r600', '-dpdf', '-painters', filename);

    
    close all
    clear title