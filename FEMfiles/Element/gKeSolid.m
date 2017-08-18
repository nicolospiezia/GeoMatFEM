function [ge,Ke,Se,be,PLe,PQe,FLAGe,RMAPe] = gKeSolid(Xe,Ge,ue,due,nip,bne,PLne,nonlin)
% -------------------------------------------------------------------------
% gKeSolid: 
%   Evaluates internal force vector, stresses, strains
%   and tangent stiffness for an hyperelastoplastic 8/12/16/20 nodes 
%   solid hexahedral element or 4/8 nodes quadrilateral
%   element in plane strain.
%
% Syntax: 
%   [ge,Ke,Se,be,PLe,FLAGe] = gKeSolid(Xe,Ge,ue,due,nip,bne,PLne,nonlin)
%
% Input:
%   Xe   : Element nodal coordinate array.
%   Ge   : Element property vector.
%   ue   : Element nodal displacements.
%  due   : Element nodal incremental displacements.
%  nip   : Number of Gauss points.
%  bne   : Left Cauchy-Green tensor at step n.
%  PLne  : Element variable for plastic algorithm at step n.
% nonlin : Flag for non-linear geometry
%
% Output: 
%   ge    : Internal force vector.
%   Ke    : Tangent stiffness matrix
%   Se    : Element stress (Cauchy) vector.
%   be    : Left Cauchy-Green tensor at stp n+1.
%   PLe   : Element variable for plastic algorithm at time n+1.
%  FLAGe  : Element flag for plasticity at GP.
%
% Date: 31/10/2014
%   Version 3.0   
% Created by: Nicolò Spiezia
% Updated by: Nico De Marchi
% -------------------------------------------------------------------------

% Number of nodes per element and dof's per node
ne  = size(Xe,1);
dof = size(Xe,2);

% Gauss point coordinates
[r,w] = Gauss(nip);

% Initialize vectors and matrix
[Se,be] = deal(zeros(2*dof,nip^dof));
ge = zeros(ne*dof,1);
PLe   = zeros(5,nip^dof);
PQe   = zeros(3,nip^dof);
FLAGe = false(nip^dof+1,1);
RMAPe   = zeros(15,1);
Ke = zeros(ne*dof);
P0=0.0;
%
%theta2 = 0:1:180;
%determ = zeros(length(theta2),1);

% loop reduction for 2D integration
if dof == 2
    nipk=1;
    wk=ones(1,nip);
else
    nipk=nip;
    wk=w;
end 

% Gauss integration of internal force vector.
ip = 0;
for i = 1:nip
    for j = 1:nip
        for k = 1:nipk
            
            % Gauss point number
            ip = ip + 1;
            
            % Plastic variables
            epspvn = PLne(1,ip);
            epspsn = PLne(2,ip);
            csin   = PLne(3,ip);
            
            if nonlin == 0
            %--------------------------------------------------------------               
            % Element Jacobian matrix,strain interpolation matrix,
            % total deformation gradient and Incremental deformation gradient
            [Jt,B,~,~,~,G] = BSolid(Xe,r(i),r(j),r(k));
            detJt=det(Jt);
            Jnew=1;
            
            % Evaluate strain
            deps=B*due;
            
            if dof == 2
                epsn = [bne(1,ip)     bne(3,ip)   0.0 
                        bne(3,ip)     bne(2,ip)   0.0 
                        0.0         0.0     bne(4,ip)];
                  
                Deps = [deps(1)     deps(3)/2   0.0 
                        deps(3)/2   deps(2)     0.0 
                        0.0         0.0         0.0  ];
                  
            else
                
                epsn = [bne(1,ip)     bne(4,ip)   bne(6,ip)
                        bne(4,ip)     bne(2,ip)   bne(5,ip) 
                        bne(6,ip)     bne(5,ip)   bne(3,ip)];
                    
                Deps = [deps(1)     deps(4)/2   deps(6)/2 
                        deps(4)/2   deps(2)     deps(5)/2
                        deps(6)/2   deps(5)/2   deps(3)  ];
            end
            
            % Infinitesimal trial strain tensor
            EeTr= epsn+Deps;
            
            [dirpr,eigva] = eig(EeTr);
            epseTr=[eigva(1,1) eigva(2,2) eigva(3,3)]';
            
            % Compute principal Kirchhoff tension and plastic variables
            [tenspr,P,Q,epse,epspv,epsps,csi,aep,FLAGe(ip,1),fac,rmap] = tens(Ge,epseTr,epspvn,epspsn,csin);
            
            % Compute Elastic strain tensor at n+1
            Ee = zeros(3);
            
            for ii=1:3

                Ee(1:3,1:3) = Ee(1:3,1:3)+(epse(ii))*dirpr(:,ii)*dirpr(:,ii)';    

            end
            
            Be=Ee;
            
            else
            %--------------------------------------------------------------    
            % Element Jacobian matrix,strain interpolation matrix,
            % total deformation gradient and Incremental deformation gradient
            [Jt,B,~,FT,fduT,G] = BSolid(Xe,r(i),r(j),r(k),ue,due);
            detJt=det(Jt);
            Jnew=det(FT);
            
            % Control Jacobian 
            if Jnew<0    
            FLAGe(nip^dof+1,1) = 1; % Warning for negative determinant of F
            end
            
            % Recall variables at time step n
            if dof == 2
                Ben = [bne(1,ip)  bne(3,ip)  0.0
                       bne(3,ip)  bne(2,ip)  0.0
                       0.0          0.0      bne(4,ip)];
            else
                
                Ben = [bne(1,ip)  bne(4,ip)  bne(6,ip)
                       bne(4,ip)  bne(2,ip)  bne(5,ip)
                       bne(6,ip)  bne(5,ip)  bne(3,ip)];
                % ---------------------------------------------------------   
                % clean incremental f - ATTENZIONE 
                for ii=1:3
                    for jj=1:3
                        r0=abs(fduT(ii,jj));
                        r1=abs(1-r0);
                        if r0<1e-10
                            fduT(ii,jj)=0.0;
                        end
                        if r1<1e-10 && fduT(ii,jj)>0
                            fduT(ii,jj)=1.0;
                        elseif r1<1e-10 && fduT(ii,jj)<0
                            fduT(ii,jj)=-1.0;
                        end
                    end
                end
                % ---------------------------------------------------------
            end
            
            % Compute Trial left cauchy-Green
            BeTr = fduT'*Ben*fduT;
            
            % Compute principal deformation and direction
            [dirpr,eigva] = eig(BeTr);
            
            lambar = [eigva(1,1) eigva(2,2) eigva(3,3)]';
            epseTr = 0.5*log(lambar); 
            
            % Compute principal Kirchhoff tension and plastic variables
            [tenspr,P,Q,epse,epspv,epsps,csi,aep,FLAGe(ip,1),fac,rmap] = tens(Ge,epseTr,epspvn,epspsn,csin);
            
            % Compute left Cauchy Green tensor at n+1
            stretpr(:,1) = exp(epse(:,1));
            
            Be = zeros(3);
            
            for ii=1:3
                
                Be(1:3,1:3) = Be(1:3,1:3)+(stretpr(ii))^2*dirpr(:,ii)*dirpr(:,ii)';    

            end
            %--------------------------------------------------------------
            end
            
            % Compute Kirchhoff tension tensor = Cauchy tension tensor in
            % small strain
            TTe = zeros(3);
            
            for ii=1:3
                    
                TTe(1:3,1:3) = TTe(1:3,1:3)+tenspr(ii)*dirpr(:,ii)*dirpr(:,ii)';    

            end
            
            if dof == 2
                
                    Te = [TTe(1,1) TTe(2,2) TTe(1,2) TTe(3,3) ]';
                
                % Compute Cauchy tension
                Se(:,ip) = Te/Jnew;
                
                % Evaluate internal force
                ge = ge + w(i)*w(j)*wk(k)*B'*Se(1:3,ip)*detJt;
                
                % Allocate elastic left Cauchy-Green component
                % or elastic strain tensor in small stain
                be(:,ip) = [Be(1,1); Be(2,2); Be(1,2); Be(3,3)];
            
            else
                
                Te = [TTe(1,1) TTe(2,2) TTe(3,3) TTe(1,2) TTe(2,3) TTe(1,3)]';
           
                % Compute Cauchy tension
                Se(:,ip) = Te/Jnew;
           
                % Evaluate internal force
                ge = ge + w(i)*w(j)*wk(k)*B'*Se(:,ip)*detJt;
            
                % Allocate elastic left Cauchy-Green component
                % or elastic strain tensor in small stain
                be(:,ip) = [Be(1,1); Be(2,2); Be(3,3); Be(1,2); Be(2,3); Be(1,3)];
                
            end
            
            % Allocate internal state variables for plastic algorithm
            if Ge(1)==4 || Ge(1)==6
                P0=P+Q^2/(Ge(7)^2*P);
            end
            PLe(:,ip)= [epspv; epsps; csi; 0 ;Jnew];
            PQe(:,ip) = [P;Q;P0];
            if ip==1
                % save return mapping residual r for gauss point 1
                RMAPe=rmap;
            end
           
            %  Compute tangent operator
            
            if nonlin == 0
            %--------------------------------------------------------------    
                % Recall consistent tangent algorithm 
                [a,~] = alfa(aep,tenspr,epseTr,dirpr,dof,nonlin,fac);
                   %alfatens
            else
            %--------------------------------------------------------------    
                % Compute tangent operator alfa
                [ALFAe,~] = alfa(aep,tenspr,lambar,dirpr,dof,nonlin,fac);
                       %alfatens
                % Recall tensor sigma
                [sigmaI]= prodton(Se(:,ip),dof);
            
                % Recall consistent tangent algorithm
                a = ALFAe/Jnew-sigmaI;
            end            
            %--------------------------------------------------------------
            
            Ke = Ke + w(i)*w(j)*wk(k)*(G'*a*G)*detJt;
            
%--------------------------------------------------------------------------            
            % Find minimum determinant of Eulerian acustic tensor
            
            %for nn=1:length(theta2)
            %    nvet = [cosd(theta2(nn)); sind(theta2(nn)); 0 ];

            %    ae = AEEu(alfatens,TTe,nvet);

            %    determ(nn)= det(ae);
            %end

            %BIF   = min(determ);
            
            %PLe(4,ip)= BIF;
            
%--------------------------------------------------------------------------
        end
    end
end
end
