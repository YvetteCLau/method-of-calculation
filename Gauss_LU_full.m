%3.b)
function [x] = Gauss_LU_full(B, p, q, b)
    n = size(B, 1); 

    %根据行置换 p 调整右端项 b
    Pb = b(p);

    %前向替代求解 Ly = Pb
    L = eye(n);  
    for i = 2:n
        for j = 1:i-1
            L(i, j) = B(i, j); 
        end
    end
    y = LowerTri(L, Pb); 

    %后向替代求解 Uz = y
    U = triu(B);  
    z = UpperTri(U, y); 

    %根据列置换向量 q 调整解 z
    x = zeros(n, 1);
    x(q) = z; 
end