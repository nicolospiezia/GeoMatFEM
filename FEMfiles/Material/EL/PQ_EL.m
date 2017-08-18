function [P,Q] = PQ_EL(epsev,epses,K,G)
%--------------------------------------------------------
% PQ Compute the mean and deviatoric invariants of Kirchhoff stress
%
% Date: 25/05/2015
%   Version 1.0   
% 
% Created by: Nicolò Spiezia
%--------------------------------------------------------

P = K*epsev;
Q = 3*G*epses;
end

