%第三题
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


%%
%3.a)函数检验
clear
clc

A = [2, 3, 1;
     4, 1, 3;
     5, 2, 1];

[A_factored, p, q] = lu_factorization_full(A);

n = size(A, 1);
P = eye(n);
P = P(p, :);  %根据行置换向量生成 P

Q = eye(n);
Q = Q(:, q);  %根据列置换向量生成 Q

L = eye(n);
U = triu(A_factored);  %提取上三角部分
for i = 2:n
    for j = 1:i-1
        L(i, j) = A_factored(i, j);  %提取下三角部分
    end
end

disp('验证 PAQ = LU：');
disp(P * A * Q - L * U);  %结果应为零矩阵

%输出
disp('置换矩阵 P：');
disp(P);
disp('置换矩阵 Q：');
disp(Q);
disp('下三角矩阵 L：');
disp(L);
disp('上三角矩阵 U：');
disp(U);


%%
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


%%
%3.b)用该函数求解2.b)中的方程
clear
clc

A = [51, 49, 23, 4, 68, 71, 74;
     9, 58, 49, 88, 14, 90, 50;
     26, 24, 62, 91, 72, 89, 48;
     80, 46, 68, 79, 11, 34, 90;
     3, 96, 40, 10, 65, 70, 61;
     92, 55, 37, 26, 49, 20, 62;
     73, 52, 98, 34, 78, 4, 86];

b1 = [340; 358; 412; 408; 345; 341; 425];
b2 = [1518; 1584; 1854; 1581; 1527; 1216; 1623];

[B, p, q] = lu_factorization_full(A);

%求解 Ax = b1
x1 = Gauss_LU_full(B, p, q, b1);

%求解 Ax = b2
x2 = Gauss_LU_full(B, p, q, b2);

%输出
disp('方程组 Ax = b1 的解 x1：');
disp(x1);

disp('方程组 Ax = b2 的解 x2：');
disp(x2);

%验证解的正确性
disp('验证 A * x1 ≈ b1：');
disp(A * x1 - b1);  % 结果应接近零向量

disp('验证 A * x2 ≈ b2：');
disp(A * x2 - b2);  % 结果应接近零向量