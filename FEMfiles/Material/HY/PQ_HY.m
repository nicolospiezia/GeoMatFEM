function [P,Q] = PQ_HY(epsev, epses, P0, alfa, kappa, epsev0, mu0)
%--------------------------------------------------------
% PQ Compute the mean and deviatoric invariants of Kirchhoff stress
%
% Reference:
% Borja R.I., Tamagnini C., Cam-Clay plasticity, part III: Estension of the
% infinitesimal model to include finite strain, CMAME 155, (1998), 73-95.
%
% Date: 25/10/2013
%   Version 1.0   
% 
% Created by: Nicolò Spiezia
%--------------------------------------------------------

OMEGA = -(epsev-epsev0)/kappa;

P = P0*exp(OMEGA)*(1+(3*alfa)*(epses^2)/(2*kappa));
Q = 3*(mu0-alfa*P0*exp(OMEGA))*epses;

end

