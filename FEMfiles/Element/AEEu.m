function ae = AEEu(alfa,tau,n)
%AEEU Summary of this function goes here
%   Detailed explanation goes here

I = eye(3);
ae=zeros(3,3);

for i=1:3
    for j =1:3
        
        sum = 0;
        
        for k = 1:3
        for l = 1:3
        
            sum = sum+n(k)*(alfa(i,k,j,l)-tau(i,l)*I(j,k))*n(l);
        
        end
        end
        
        ae(i,j) = sum;
    end
end
end



