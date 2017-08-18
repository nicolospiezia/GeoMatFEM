%--------------------------------------------------------
% File: Driver.m
%
% Driver for nonlinear analysis of Hyper-elasto-plastic
% problem with 4/8 nodes quadrilateral elements or 8/12/16/20 nodes 
% solid hexahedral element.
%
% Input: 
%
%   X    : Node coordinate array.
%   T    : Topology array.
%   G    : Material property array.
%   C    : Constraints.
%   P    : Prescribed nodal loads.
%  nip   : Number of Gauss points.
%  imax  : Maximum number of iterations.
%  epsr  : Error tolerence on residual.
%  incr  : Coefficients for load increment.
% nonlin : Flag for non-linear geometry
%
% Output:
%
%   U    : Nodal displacement array.
%   F    : Nodal force array.
%   S    : System stress array.
%
%  -----------------------------------------------------------------------  *
%   COPYRIGHT STATEMENT 													*
%																			*	
%   Copyright (C) 2014  Nicolo Spiezia                                      *
%																			*
%   This program is free software: you can redistribute it and/or modify    *
%   it under the terms of the GNU General Public License as published by    *
%   the Free Software Foundation, either version 3 of the License, or       *
%   (at your option) any later version.                                     *
%																			*
%   This program is distributed in the hope that it will be useful,         *
%   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
%   GNU General Public License for more details.                            *
%																			*
%   You should have received a copy of the GNU General Public License       *
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
%																			*
%   E-mail: nicolospiezia@gmail.com											*
% ------------------------------------------------------------------------  *
% Date: 31/10/2014
%   Version 3.0   
% Updated by: Nico De Marchi 
% E-mail: nicodm9@hotmail.it
% ------------------------------------------------------------------------

disp('Load initial variables from elatic analysis: no=0; yes=1');
flagin=input('flagin=' );

tic;
parpool ('local',8);
%matlabpool ('open',2);

% Add folder and subfolders to Matlab search path
addpath(genpath('..\FEMfiles'));%path(path,'..\FEMfiles'); %path(path,'FEMfiles');

% System parameters
nmax= size(incrv,2);
dof = size(X,2);
nno = size(X,1);
ndof = nno*dof;
ne = size(T,2)-1;
nelem = size(T,1);

% Initialize load and displacement vectors
[u,Du,du,f,p,sup,sup1,sup2,V] = deal(zeros(ndof,1));
[U,F] = deal(zeros(ndof,nmax+1));
[brec,Srec] = deal(zeros(2*dof,nip^dof,nelem,nmax+1));
[PLrec] = zeros(5,nip^dof,nelem,nmax+1);
[PQrec] = zeros(3,nip^dof,nelem,nmax+1);
[FLAGrec] = false(nip^dof+1,nelem,nmax+1);
[RESrec] = zeros(imax,nmax+1);

% Line search variables
f1=1.0;
f2=0.0;
alfa=1.0;

%--------------------------------------------------------------------------

% Inizialize elastic and plastic variable record
if flagin==0
    [brec(:,:,:,1),PLrec(3,:,:,1)] = InizbPL(nip,dof,ne,nelem,nonlin,T,G);
else
    load var0_iniz_es44.mat;
    PLrec(3,:,:,1)=PL(3,:,:);
    PQrec(:,:,:,1)=PQ(:,:,:);
    Srec(:,:,:,1)=S(:,:,:);
    brec(:,:,:,1)=leftc(b,S,nip,dof,nelem,nonlin,T,G);
end
%--------------------------------------------------------------------------

% Identify unconstrained/constrained dofs
[iu,ic,id,uimp] = Const(C,dof,ndof);
%u(id,1) = uimp; U(:,1)=u(:,1);

% Evaluate volumetric load --->  WARNING <--->  WARNING <--->  WARNING <---

V = VSolid(V,T,X,G,av,nip,1);
V = VSolid(V,T,X,G,av,nip,2);
V = VSolid(V,T,X,G,av,nip,3);
V = VSolid(V,T,X,G,av,nip,4);

% Set Face load
sup1 = PSolid(sup,T,X,pres1,nip);
sup2 = PSolid(sup,T,X,pres2,nip);

% Set Nodal load
%p = SetLoad(p,P);

for n = 1:nmax
    
    fprintf('\n Load step %d \n',n);
    
    % initialize displacement increment
    Du = zeros(ndof,1);
    
    % update load vector
    %f = (sup+p)*incrf(n)+V*incrv(n);
    f = sup1*incrf1(n)+sup2*incrf2(n)+V*incrv(n); 
    % Recall plastic variable
    PLn = PLrec(:,:,:,n);
    bn  = brec (:,:,:,n);
       
    % equilibrium iterations
    for iter = 1:imax
               
        % evaluate internal forces
        g = zeros(ndof,1);
        [g,K,S,b,PL,PQ,FLAG] = gKSolid(g,T,X,G,u,Du,nip,bn,PLn,nonlin); 
        % residual
        r = f - g;
        
                % check for convergence
                % Total force norm
                nrmf=norm(f(iu,1));
                if nrmf==0
                    nrmf=1;
                end
                f2 = norm(r(iu,1));
                RESrec(iter,n) = f2/nrmf;
                fprintf('\n Residual norm %d \n',RESrec(iter,n));
                if RESrec(iter,n) < epsr || RESrec(iter,n) > 1.0e15
                    break
                else
                % ------- The Line Search ---------
                if  f2<f1 || iter==1
                alfa=1.0;  
                else
                alfa=alfa^2*f1/(2*(f2+f1*alfa-f1));
                end
                f1=f2;
                % ---------------------------------
                  
                % solve for displacement increment
                du(iu) = K(iu,iu)\r(iu);
                
                Du = Du + alfa*du;
                
                end
    end
    %
    if iter == imax
        fprintf('\n No convergence in load step %d \n',n);
        % Stop program!!!
        %break
    end
     
    % update displacements
    u = u + Du;
    F(:,n+1) = f;
    U(:,n+1) = u;
    
    % record plastic variables
    PLrec(:,:,:,n+1) = PL;
    PQrec(:,:,:,n+1) = PQ;
    FLAGrec(:,:,n+1) = FLAG;
    Srec(:,:,:,n+1) = S;
    brec(:,:,:,n+1) = b;
end

delete(gcp)
%matlabpool ('close');
tempo=toc;
save(filename);
disp('Save elastoplastic variables?? no=0; yes=1');
flagout=input('flagout=' );
if flagout~=0 
save('var0_iniz_es44.mat','PL','PQ','b','S');
end
disp(' DONE! :-)');
load handel.mat;
nBits = 16;
sound(y,Fs,nBits);
