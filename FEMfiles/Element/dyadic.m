function C = dyadic(A,B)
%DYADIC product of two second order tensor
%--------------------------------------------------------------------------
m = size(A,1);
n = size(A,2);
o = size(B,1);
p = size(B,2);
C=zeros(m,n,o,p);

for i=1:m
    for j =1:n
        for k=1:o
            for l=1:p
       
                C(i,j,k,l)= A(i,j)*B(k,l);
            
            end
        end
    end
end

