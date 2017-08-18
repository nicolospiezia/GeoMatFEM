function [P,Q,epsev,epses,epspv,epsps,Pi,Dep] = rm_DP(epsevTR,epsesTR,Ptr,Qtr,...
                                                       epspvn,epspsn,kappa,mu0,alfa,P0,epsev0,m,c0,mbar,Pin)

                                                   
% Parameter
imax = 20;
toll = 10e-6;
NORMErec = zeros(imax);

% Inizialize variables
x = zeros(3,1);
x(1,1) = epsevTR;   % epsev
x(2,1) = epsesTR;   % epses
x(3,1) = 0;         % dgamma


P = Ptr;
Q = Qtr;

for iter = 1:imax

    % evaluate residual                

    F = Q-m*P-c0;

    dPG = -mbar;
    dQG = 1;

    r = [x(1) - epsevTR+x(3)*dPG
         x(2) - epsesTR+x(3)*dQG
               F                ];

       if iter == 1
       r0 = norm(r);
       end

    NORMErec(iter) = norm(r)/norm(r0);

    % check for convergence
    if norm(r) < toll*r0

        break
    else

        % evaluate tangent matrix
    Atang = Atang_CCDP1( x(1,1),x(2,1),P,kappa,mu0,alfa,P0,epsev0,m,mbar);

        % solve for displacement increment
        dx = - (Atang\r);

        x = x + dx;                 

        % Update 

        [P,Q]=PQ_HY(x(1,1),x(2,1),P0,alfa,kappa,epsev0,mu0);

        epspv = epsevTR-x(1)+epspvn;
        epsps = epsesTR-x(2)+epspsn;


    end
end

    if iter == imax
    fprintf('No convergence of RM')
    end


%Update variable
epsev = x(1,1);
epses = x(2,1);
dgamma = x(3,1);

Pi = Pin;

% Compute matrix Dep for DP non-associative

Dep = DEP_CCDP1( epsev,epses,P,kappa,mu0,alfa,P0,epsev0,m,mbar);
end

