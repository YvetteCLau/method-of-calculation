%3.a)
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