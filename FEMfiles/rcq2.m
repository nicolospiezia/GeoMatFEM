function [r,c,q] = rcq2(Ae,Te,dof)
%--------------------------------------------------------
%rcq:
%
% Syntax:
%  [r,c,q] = rcq(Ae,Te,dof)
%
% Input:
%   A  :  Global matrix.
%   Ae :  Element matrix.
%   Te :  Element topology vector.
%   dof:  Degrees of freedom per node.
%
% Output:
%   r  :  Row index vector.
%   c  :  Column index vector.
%   q  :  Non-zero element vector. 
%
% Date: 31/10/2014
%   Version 2.0   
% Created by: Nico De Marchi
% 
%--------------------------------------------------------

% Number of element nodes
ne = size(Te,2) - 1; %   
na = size(Ae,1);
[r,c,q]=deal(zeros(numel(Ae),1));

% Define global address vector for element dofs

ig = zeros(1,ne*dof);

for i = 1:ne
    
  ig((i-1)*dof+1:i*dof) = ((Te(i)-1)*dof+1:Te(i)*dof); 
  
end

 count=0;
 for i=1:na
     for j=1:na
         count=count+1;
         r(count,1)=j;
         c(count,1)=i;
         q(count,1)=Ae(j,i);
     end
 end
 
 %[r,c,q] = find(Ae);
 
 %if length(c)<numel(Ae)
     
 %    z=zeros(numel(Ae)-length(c),1);
 %    c=[c;z];
 %    r=[r;z];
 %    q=[q;z];
 %end
 for j=1:length(r)
     for i=1:length(ig)   
         if c(j)==i
             c(j)=ig(i);
             break
         end 
     end
     for i=1:length(ig)
         if r(j)==i
             r(j)=ig(i);
             break
         end 
     end 
 end
 % attenzione: non è possibile fare un unico ciclo "for i=1:length(ig)"
 % perchè c'è il break dentro i cicli if
 