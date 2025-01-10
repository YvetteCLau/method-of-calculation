%2
function [A, p] = lu_factorization(A)
    n = size(A, 1);  
    p = (1:n)';      
    for k = 1:n-1
        %寻找列主元并记录行交换信息
        [~, idx] = max(abs(A(k:n, k)));  
        idx = idx + k - 1;              

        if A(idx, k) == 0
            error('矩阵奇异，无法进行 LU 分解');
        end
        %若主元不在当前行，则进行行交换
        if idx ~= k
            A([k, idx], :) = A([idx, k], :);
            p([k, idx]) = p([idx, k]);
        end

        %计算消元因子，并存储于 A 的下三角部分
        for i = k+1:n
            A(i, k) = A(i, k) / A(k, k);  
            
            %消元：更新 A(i, j)
            A(i, k+1:n) = A(i, k+1:n) - A(i, k) * A(k, k+1:n);
        end
    end

    if A(n, n) == 0
        error('矩阵奇异，无法进行 LU 分解');
    end
end