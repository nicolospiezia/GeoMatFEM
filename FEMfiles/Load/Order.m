function [X2,T2] = Order(X,face)
%--------------------------------------------------------
% Order:
%   Creates local matrix X2 and local topology array
%
% Syntax:
%   [X2,T2] = Order(X,face)
%
% Input:
%   X    : System nodal coordinate array
%  face  : Hexaedral face load
%
% Output:
%  X2    : Local matrix 
%  T2    : Local topology array.
%
% Date: 31/10/2014
%   Version 3.0   
% Updated by: Nico De Marchi
%--------------------------------------------------------
%  Nodes numbering:
%
%   4---(6)---3 
%   |         | 
%  (8)       (7)
%   |         | 
%   1---(5)---2 
%--------------------------------------------------------

node= size(X,1);
dof = size(X,2);

if dof ==2
    if node >= 4
        switch face
            case 1
                T2 = [1 2];
            case 2
                T2 = [2 3];
            case 3
                T2 = [3 4];
            case 4
                T2 = [4 1];
        end
    end
    
    if node >= 6
        switch face
            case 1
                T2 = [T2 5];
            case 3
                T2 = [T2 6];
        end
    end
    
    if node >= 8
        switch face
            case 2
                T2 = [T2 7];
            case 4
                T2 = [T2 8];
        end
    end
    
else
    if node >= 8
        switch face
            case 1
                T2 = [1 2 3 4];
            case 2
                T2 = [2 3 7 6];
            case 3
                T2 = [5 6 7 8];
            case 4
                T2 = [1 4 8 5];
            case 5
                T2 = [4 3 7 8];
            case 6
                T2 = [1 2 6 5];
        end
    end

    if node >= 12
        switch face
            case 1 
                T2 = [T2 9 10];
            case 3
                T2 = [T2 12 11];
            case 5
                T2 = [T2 10 11]; 
            case 6
                T2 = [T2 9 12];
        end
    end
    
    if node >= 16
        switch face
            case 1 
                T2 = [T2 13 14];
            case 2
                T2 = [T2 13 16];
            case 3
                T2 = [T2 16 15];
            case 4
                T2 = [T2 14 15];
        end
    end
    
    if node >= 20
        switch face
            case 2 
                T2 = [T2 19 18];
            case 4
                T2 = [T2 20 17];
            case 5
                T2 = [T2 19 20];
            case 6
                T2 = [T2 18 17];
        end
    end
end

X2=X(T2,:); 
T2=[T2 0];

end