function [P,Q,epsev,epses,epspv,epsps,Pi,Dep,flag] = rm_CC(epsevTR,epsesTR,Ptr,Qtr,...
                                                       epspvn,epspsn,kappa,mu0,alfa,P0,epsev0,Pi0,A,pp,epstar,m,c0,mbar,Pin)

% Parameter
imax = 20;
toll = 10e-6;
NORMErec = zeros(imax);

% --- CAP CC Return Mapping
RM = 1;  
                                               
                                                 
% Inizialize variables
x = zeros(3,1);
x(1,1) = epsevTR;   % epsev
x(2,1) = epsesTR;   % epses
x(3,1) = 0;         % dgamma 

Pi = Pin;
P = Ptr;
Q = Qtr;

for iter = 1:imax

    % evaluate residual                
    B2 = m^2*Pi^2-m^2*A^2+2*m*c0*Pi+c0^2;

    F = B2*(P-Pi)^2+(A*Q)^2-B2*A^2;

    dPF = 2*B2*(P-Pi);
    dQF = 2*(A^2)*Q;

    r = [x(1) - epsevTR+x(3)*dPF
         x(2) - epsesTR+x(3)*dQF
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
    Atang = Atang_CCDP2( x(1,1),x(2,1),x(3,1),P,Q,Pi,kappa,mu0,...
                    alfa,P0,epsev0,pp,Pi0,epstar,A,B2,m,c0,epspvn,epsevTR);

        % solve for displacement increment
        dx = - (Atang\r);

        x = x + dx;                 
    end
    
        [P,Q]=PQ_HY(x(1,1),x(2,1),P0,alfa,kappa,epsev0,mu0);

        epspv = epsevTR-x(1)+epspvn;
        epsps = epsesTR-x(2)+epspsn;


        Pi = Pi0*(epstar/(epstar-epspv))^pp;  
    

end

if iter == imax
    fprintf('No convergence of CC-RM')

end

% Compute ptilde
b = sqrt(m^2*Pi^2-m^2*A^2+2*m*c0*Pi+c0^2); 
ptilde = Pi-(mbar*A*b)/(sqrt((b^4/A^2)+mbar^2*b^2));     


% --- CAP CCDP Return Mapping
if P > ptilde  

RM = 2;

% Inizialize variables
x = zeros(3,1);

x(1,1) = epsevTR;   % epsev
x(2,1) = epsesTR;   % epses
x(3,1) = 0;         % dgamma 

Pi = Pin;
P = Ptr;
Q = Qtr;   
    
for iter = 1:imax 
        
    B2 = m^2*Pi^2-m^2*A^2+2*m*c0*Pi+c0^2;
    F = B2*(P-Pi)^2+(A*Q)^2-B2*A^2;

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
        Atang = Atang_CCDP21( x(1,1),x(2,1),P,Q,Pi,kappa,mu0,...
                        alfa,P0,epsev0,pp,Pi0,epstar,A,B2,m,mbar,c0,epspvn,epsevTR);

            % solve for displacement increment
            dx = - (Atang\r);

            x = x + dx;  

        end

        [P,Q]=PQ_HY(x(1,1),x(2,1),P0,alfa,kappa,epsev0,mu0);

        epspv = epsevTR-x(1)+epspvn;
        epsps = epsesTR-x(2)+epspsn;


        Pi = Pi0*(epstar/(epstar-epspv))^pp;     
        
        
   end

end 

if iter == imax
    fprintf('No convergence of RM')
end


%Update variable
epsev = x(1,1);
epses = x(2,1);
dgamma = x(3,1);

% ------------------------------------------------------------------

% Compute matrix Dep  
if RM == 1;
Dep = DEP_CCDP2( epsev,epses,dgamma,P,Q,Pi,kappa,mu0,alfa,...
             P0,epsev0,pp,Pi0,epstar,A,B2,m,c0,epspvn,epsevTR);
flag = 22;
else
Dep = DEP_CCDP21( epsev,epses,P,Q,Pi,kappa,mu0,alfa,...
             P0,epsev0,pp,Pi0,epstar,A,B2,m,mbar,c0,epspvn,epsevTR);
flag = 21;
end


end

