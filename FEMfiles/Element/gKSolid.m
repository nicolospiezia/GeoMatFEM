function [g,K,S,b,PL,PQ,FLAG,RMAP] = gKSolid(g,T,X,G,u,Du,nip,bn,PLn,nonlin)
%--------------------------------------------------------
%   gKSolid: 
%   Creates and assembles internal force vector, stress/strain matrix
%   and tangent stiffness matrix for a group of
%   8/12/16/20 nodes solid hexahedral elements or 4/8 nodes quadrilateral 
%   elements in plane strain.
%
% Syntax:
%   [g,K,S,b,PL,FLAG] = gKSolid(g,T,X,G,u,Du,nip,bn,PLn,nonlin)
%
% Input:
%   g    : Existing global internal force vector.
%   X    : System nodal coordinate array.
%   T    : System topology array.
%   G    : Element property array.
%   u    : System nodal displacement vector.
%   Du   : System nodal increment displacement vector.
%  nip   : Number of Gauss points.
%   bn   : Left Cauchy-Green tensor at step n 
%  PLn   : System variable for plastic algorithm at step n.
% nonlin : Flag for non-linear geometry
%
% Output:
%   g    : New global internal force vector.
%   S    : Global stress matrix.
%   K    : Tangent stiffness compact matrix
%   b    : Left Cauchy-Green tensor
%   PL   : System variable for plastic algorithm at time n+1.
%  FLAG  : System flag for plasticity.
%
% Date: 31/10/2014
%   Version 3.0   
% Created by: Nicolò Spiezia
% Updated by: Nico De Marchi
%--------------------------------------------------------
%load gong.mat;
% Determine number of nodes per element
ne = size(T,2)-1;
dof = size(X,2);
nelem = size(T,1);
nc=(ne*dof)^2;

% reshape global displacement vector
U = reshape(u,size(X,2),size(X,1))';
DU = reshape(Du,size(X,2),size(X,1))';

% Initialize arrays
[S,b] = deal(zeros(2*dof,nip^dof,nelem)); 
[PL] = zeros(5,nip^dof,nelem);
[PQ] = zeros(3,nip^dof,nelem);
[FLAG] = false(nip^dof+1,nelem);
[RMAP] = zeros(15,nelem);
ge = zeros(ne*dof,nelem);
[r,c,q]=deal(zeros(nc,nelem));

% Loop over elements
parfor j = 1:nelem
     
% Define element arrays 
  Te = T(j,:);
  Xe =  X(Te(1:ne),:);              
  Ue =  U(Te(1:ne),:);
  ue = reshape(Ue',ne*dof,1);
  
  DUe =  DU(Te(1:ne),:);
  due = reshape(DUe',ne*dof,1);
  
  Ge = G(Te(ne+1),:);
  bne = bn(:,:,j);
  PLne = PLn(:,:,j);
  
  % Evaluate internal forces for element
  [ge(:,j),Ke,Se,be,PLe,PQe,FLAGe,RMAPe] = gKeSolid(Xe,Ge,ue,due,nip,bne,PLne,nonlin);
  [r(:,j),c(:,j),q(:,j)] = rcq2(Ke,Te,dof);
  
  if FLAGe(nip^dof+1,1) == 1
  fprintf('\n Warning: Det(F)<0 - Element %d \n',j)
  %sound(y);
  end
  
  % Assemble into global arrays 
  S(:,:,j) = Se;
  b(:,:,j) = be;
  PL(:,:,j) = PLe;
  PQ(:,:,j) = PQe;
  FLAG(:,j) = FLAGe;
  RMAP(:,j) = RMAPe;
end

 % Assemble into global arrays
for j = 1:nelem
   Te = T(j,:);
   g = Assem(g,ge(:,j),Te,dof);
end

% Assemble global stiffness matrices
r=reshape(r,nc*nelem,1);
c=reshape(c,nc*nelem,1);
q=reshape(q,nc*nelem,1);
%xx=r==0;
%r(xx==1)=[];
%c(xx==1)=[];
%q(xx==1)=[];
K=sparse(r,c,q);

end
           