function [P,Q,De,flag,fac] = tens_EL(Ge,epsev,epses)
%--------------------------------------------------------
% tensVM: 
%   Compute the Kirchhoff principal tension according to 
%   associative or non associative Drucker-Prager yield criterion.
%
% Syntax:
%   [tenspr,P,Q,epse,epspv,aep,flag] = tensDP(Ge,defepr,epspvn)
%
% Input:
% Ge      : Material property.
% defepr  : Elastic principal deformation.
% epspvn  : Plastic volumetric deformation at step n.
%
% Output: 
% tenspr : Vector of Kirchhoff principal tension. [beta(1) beta(2) beta (3)]'
% P      : Mean normal Kirchhoff stress
% Q      : Deviatoric Kirchhoff stress      
% epse   : Vector of elastic principal strain. [epse(1) epse(2) epse (3)]'
% epspv  : Plastic volumetric deformation 
% aep    : Algorithmic stress-strain matrix in principal direction (a_AB)
% flag   : Flag for plasticity
%
% Date:
%   Version 1.0   25.05.15
%--------------------------------------------------------

% Set isotropic elastic parameter

% 3  4	  
% E  ni   

E = Ge(3);  
ni = Ge(4);

G = E/(2*(1+ni));
K = E/(3*(1-2*ni));

% Compute stress    
flag = 0;

[P,Q] = PQ_EL(epsev,epses,K,G);

De = DE_EL(K,G);

fac = G;
  
end

