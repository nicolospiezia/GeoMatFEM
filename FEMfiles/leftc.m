function [ b ] = leftc( e,s,nip,dof,nelem,nonlin,T,GG )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

b=zeros(2*dof,nip^dof,nelem);
imax=15;
toll=1.0e-18;
n0=[1.0 1.0 1.0]';

if nonlin==0
    b(:,:,:)=e(:,:,:);
else
    for i=1:nip^dof
        for j=1:nelem
            
            if dof == 2
                 S = [s(1,i,j)  s(3,i,j)  0.0
                      s(3,i,j)  s(2,i,j)  0.0
                      0.0          0.0      s(4,i,j)];
            else
                
                 S = [s(1,i,j)  s(4,i,j)  s(6,i,j)
                      s(4,i,j)  s(2,i,j)  s(5,i,j)
                      s(6,i,j)  s(5,i,j)  s(3,i,j)];
            end
            % Principal stresses
            [eigve,eigva] = eig(S);
            tensp=[eigva(1,1);eigva(2,2);eigva(3,3)];
            
            % P,Q stress invariants
            P=sum(tensp)/3;
            xsi=tensp-P*n0;
            normxsi=norm(xsi);
            Q=normxsi*sqrt(3/2);
            %
            prop=T(j,end);
            Ge=GG(prop,:);
            
            %  Volumetric and deviatoric strain invariants
            if Ge(1)==1 ||  Ge(1)==2
                E = Ge(3);  
                ni = Ge(4);
                G = E/(2*(1+ni));
                K = E/(3*(1-2*ni));
                epsev=P/K;
                epses=Q/(3*G);
            else
                mu0  = Ge(3);  
                alfa = Ge(4);
                kappa = Ge(5);
                P0 = Ge(8);
                epsev0 = Ge(10);
                epsev1=epsev0;
                % 
                for iter=1:imax   
                OMEGA = -(epsev1-epsev0)/kappa;
                epses=Q /(3*(mu0-alfa*P0*exp(OMEGA)));
                epsev=epsev0-kappa*log(P/(P0*(1+(3*alfa)*(epses^2)/(2*kappa))));
                sol=norm(epsev-epsev1);
                if sol<=toll*norm(epsev1);
                    break
                else
                    epsev1=epsev;
                end
                if iter == imax
                        fprintf('\n No convergence  \n');
                        break 
                end
                end
            end
            
            %
            if normxsi == 0
            n(:,1) = xsi; 
            else
            n(:,1) = xsi/ normxsi;  
            end
            
            % Compute principal elastic strain
            epse(:,1) = (1/3)*epsev*n0+sqrt(3/2)*epses*n(:,1);
            Ee =[(exp(epse(1,1)))^2 (exp(epse(2,1)))^2 (exp(epse(3,1)))^2]';
              
            Be = zeros(3);
            
            for k=1:3

                Be(1:3,1:3) = Be(1:3,1:3)+Ee(k)*eigve(:,k)*eigve(:,k)';    

            end
            
            if dof == 2
                b(:,i,j) = [Be(1,1); Be(2,2); Be(1,2); Be(3,3)];
            else
                b(:,i,j) = [Be(1,1); Be(2,2); Be(3,3); Be(1,2); Be(2,3); Be(1,3)];
            end
            
        end
    end
end
end