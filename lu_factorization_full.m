%3.a)
function [A, p, q] = lu_factorization_full(A)
    n = size(A, 1); 
    p = 1:n;         
    q = 1:n;         

    for k = 1:n-1
        %寻找全主元
        subMatrix = abs(A(k:n, k:n));
        [maxVal, maxIdx] = max(subMatrix(:));  %最大值及其线性索引
        [rowOffset, colOffset] = ind2sub(size(subMatrix), maxIdx); 

        row = rowOffset + k - 1; 
        col = colOffset + k - 1;  

        if maxVal == 0
            error('矩阵奇异，无法进行 LU 分解');
        end

        %行交换
        if row ~= k
            A([k, row], :) = A([row, k], :); 
            p([k, row]) = p([row, k]);        
        end

        %列交换
        if col ~= k
            A(:, [k, col]) = A(:, [col, k]);  
            q([k, col]) = q([col, k]);       
        end

        %计算消元因子并更新矩阵A
        for i = k+1:n
            A(i, k) = A(i, k) / A(k, k);  
            A(i, k+1:n) = A(i, k+1:n) - A(i, k) * A(k, k+1:n); 
        end
    end
    if A(n, n) == 0
        error('矩阵奇异，无法进行 LU 分解');
    end
end