% ------------------------------------------------------------------------
% File: es41.m
%
% Driver for import input from Straus7 Release 2.3.3 text file
%
% Warning: this script dont read pressure dof "Rz"
%
% Date: 31/10/2014
%   Version 1.0   
% Created by: Nico De Marchi
% -------------------------------------------------------------------------
clc,clear all,close all

% Open text file
dati = fopen ('es413D.txt','r');

% dim -> 2D Plane Strain or 3D problem
dim=3;

% nonlin= 0 -> Small strain; =  1 -> Large strain
nonlin = 1;

% Save command: 
filename = ['es42', '_', num2str(dim),'D', '_', num2str(nonlin)];

% Material properties
% Material properties
P0=-0.01;                        pi=-0.01; %/exp(-1);
OCR=1.0;                         Pc0=OCR*pi;
nu=0.261;                        n=0.7024;%e=2.36;
k=0.05;                          ks=k/(1+k);
lamda=0.2;                       lamdas=lamda/(1+lamda);         
epsev0=k*log(pi/P0);             %K=-pi*(1+e)/ks; 
epsev0s=ks*log(pi/P0);           %mu0=K*(3/2)*(1-2*nu)/(2*(1+nu));
%P0=-0.01/exp(epsev0s/ks);
   % Modified Cam-clay                                    % Soil & Fluid properties
if nonlin == 0 
   %G= [CC rhos   mu0 alfa k   lamda   M   Po Pc0 epsev0  - - - D n rhow   muw kx/rg ky/rg kz/rg]
    G= [4  2.7e-9 0.2 0.0  ks  lamdas  1.0 P0 Pc0 epsev0s 0 0 0 0 n 1.0e-9 0.0 1.0   1.0   1.0];  
else
   %G= [CC rhos   mu0 alfa k  lamda  M   Po Pc0 epsev0 - - - D n rhow   muw kx/rg ky/rg kz/rg]
    G= [4  2.7e-9 0.2 0.0  k  lamda  1.0 P0 Pc0 epsev0 0 0 0 0 n 1.0e-9 0.0  1.0   1.0   1.0]; 
end
clear pi P0 nu e k ks K mu0 lamdas lamda epsev0 epsev0s OCR Pc0 n
% Volumetric load av=[ax ay az]' 
% WARNING:write gravity in (ax) or (ay) in 2D Plane Strain!
av=[0 0 -1.0e4]';

% Number of integration points
nip = 2;

% Iteration parameters
imax = 15;                 % Max number of iterations
epsr = 1.0e-9;             % Error tolerence on residual 
nmax = 20;                 % Max number of load steps

% dim -> 2D Plane Strain or 3D problem
dim=3;

% Load & Gravity increments      
incrf = 0;%1.0/nmax:1/nmax:1.0;  
incrv = 1;%ones(1,nmax);

% -------------------------------------------------------------------------
 
w= fgetl(dati);
[i,j,k,l,m,n,o,p]=deal(0);
P = [ 1  0  0  0];
pres=[1  1  0  0  0];
prese=[1  1  0  ];
shear=[1  1  0  ];

while ischar(w)
    w= fscanf(dati,'%s',1);
    
    node= strfind(w,'Node');
    flagn = length(node);
    
    quad = strfind(w,'Quad');
    flagq = length(quad);
    
    hexa = strfind(w,'Hexa');
    flagh = length(hexa);
    
    cons = strfind(w,'NdFreedom');
    flagc = length(cons);
    %---------------------------------------
    x= strfind(w,'DX'); 
    y= strfind(w,'DY');
    z= strfind(w,'DZ');
    nx = length(x); 
    ny = length(y);
    nz = length(z);
    %---------------------------------------
    force = strfind(w,'NdForce');
    flagf = length(force);
    
    sp = strfind(w,'BkGlobalLoad');
    flags = length(sp);
    
    se = strfind(w,'PlEdgePressure');
    flagse = length(se);
    
    ss = strfind(w,'PlEdgeShear');
    flagss = length(ss);
    
    if flagn~=0
        % Nodal coordinates
        i=i+1;
        X(i,:) = fscanf(dati,'%*d %f %f %f',3);
        
    elseif flagq~=0
        %Topology matrix:
        j=j+1;
        nq = sscanf(w,'Quad%d',1);
        if nq == 4
            let=('%*d %*d %d %f %f %f %f');
        elseif nq==8
            let=('%*d %*d %d %f %f %f %f %f %f %f %f');
        end
        T1(j,1:nq+1) = fscanf (dati,let,nq+1);
   
    elseif flagh~=0
        %Topology matrix:
        k=k+1;
        ne = sscanf(w,'Hexa%d',1);
        switch ne
            case 8
            let=('%*d %*d %d %f %f %f %f %f %f %f %f');
            case 16
            let=('%*d %*d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
            case 20
            let=('%*d %*d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');    
        end
        T0(k,1:ne+1) = fscanf (dati,let,ne+1);
    
    elseif flagc~=0
        % Constraints, C = [ node dof u ]
        c = fscanf (dati,'%*d %d');
        w = fgetl(dati);
        x= strfind(w,'DX'); 
        y= strfind(w,'DY');
        z= strfind(w,'DZ');
        nx = length(x); 
        ny = length(y);
        nz = length(z);
        if nx~=0
            l=l+1;
            num = sscanf(w,'DX %f',1);
            nn = length(num);
            if nn~=0
            C(l,1:3)=[c 1 num];
            else
            C(l,1:3)=[c 1 0]; 
            end
        end
        if ny~=0
            l=l+1;
            num = sscanf(w,'DY %f',1);
            nn = length(num);
            if nn~=0
            C(l,1:3)=[c 2 num];
            else
            C(l,1:3)=[c 2 0]; 
            end
        end
        if nz~=0 && dim==3
            l=l+1;
            num = sscanf(w,'DZ %f',1);
            nn = length(num);
            if nn~=0
            C(l,1:3)=[c 3 num];
            else
            C(l,1:3)=[c 3 0]; 
            end
        end
    %----------------------------------
    elseif nx~=0
            l=l+1;
            num = fscanf(dati,'%f',1);
            C(l,1:3)=[c 1 num];
    elseif ny~=0
            l=l+1;
            num = fscanf(dati,'%f',1);
            C(l,1:3)=[c 2 num];
    elseif nz~=0 && dim==3
            l=l+1;
            num = fscanf(dati,'%f',1);
            C(l,1:3)=[c 3 num];
    %----------------------------------
    elseif flagf~=0
        m=m+1;
        % Nodal loads, P = [ node Fx Fy Fz ]
        P(m,1:4) = fscanf (dati,'%*s %d %f %f %f',4);           
     
    elseif flags~=0
        n=n+1;
        % Pressure load, pres = [ brick face px py pz ]
        pres(n,1:5) = fscanf (dati,'%*s %d %d %f %f %f',5);     
    
    elseif flagse~=0
        o=o+1;
        % plate edge pressure, prese = [ brick face pe ]
        prese(o,1:3) = fscanf (dati,'%*s %d %d %f',3);
    
    elseif flagss~=0
        p=p+1;
        % plate edge shear, press = [ brick face pe ]
        shear(p,1:3) = fscanf (dati,'%*s %d %d %f',3);
        
    else
        w= fgetl(dati);
    end
    
end
%
w=fclose(dati);
% -------------------------------------------------------------------------
if dim ==2
    X(:,3)=[];
    av(3,:)=[];
    P(:,4)=[];
    pres=[prese zeros((size(prese,1)),1)
          shear  ones((size(shear,1)),1)];
    
%   Re-organize the element topology matrix for 2d element
%   Nodes numbering:
%   4---(6)---3 
%   |         | 
%  (8)       (7)
%   |         | 
%   1---(5)---2 

    if nq == 4
        T=[T1(:,2:5) T1(:,1)];
    elseif nq == 8
        T=[T1(:,2:6) T1(:,8) T1(:,7) T1(:,9) T1(:,1)];
    end
else
    
%   Re-organize the element topology matrix for 3D solid element
    switch ne
        case 8
    T=[T0(:,2:9) T0(:,1)];
        case 16
    %T= [T0(:,1:9) T0(:,11) T0(:,15) T0(:,13) T0(:,10) T0(:,12) T0(:,16) T0(:,15) T0(:,17)];
    T= [T0(:,2:10) T0(:,12) T0(:,16) T0(:,14) T0(:,11) T0(:,13) T0(:,17) T0(:,16) T0(:,1)];
        case 20
    %T= [T0(:,1:9) T0(:,11) T0(:,19) T0(:,17) T0(:,10) T0(:,12) T0(:,20) T0(:,18) T0(:,13:16) T0(:,21)];
    T= [T0(:,2:10) T0(:,12) T0(:,20) T0(:,18) T0(:,11) T0(:,13) T0(:,21) T0(:,19) T0(:,14:17) T0(:,1)];
    end
end
% -------------------------------------------------------------------------

clear w i j k l m n o p c x y z nn nx ny nz nq sp se ss sv ne T0 T1 let num  ...
      node quad flagn flags hexa cons flagh flagq flagc force flagse flagss  ...
      flagf dati dim prese shear  XU ...

disp(' DONE! :-)');