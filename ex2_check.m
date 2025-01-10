%2
clear
clc

A = [2, 2, 0;
     1, 1, 2;
     2, 1, 1];
b = [6; 9; 7];

[A_factored, p] = lu_factorization(A);


n = size(A, 1);
P = eye(n);
P = P(p, :);

L = eye(n);
U = triu(A_factored);  %提取上三角部分
for i = 2:n
    for j = 1:i-1
        L(i, j) = A_factored(i, j);  %提取下三角部分
    end
end

disp('置换矩阵 P：');
disp(P);
disp('下三角矩阵 L：');
disp(L);
disp('上三角矩阵 U：');
disp(U);

disp('验证 PA = LU：');
disp(P * A - L * U);  %结果为零矩阵

%解方程组 PAx = Pb
Pb = P * b;

%前向替代求解 Ly = Pb
y = zeros(n, 1);
for i = 1:n
    y(i) = Pb(i) - L(i, 1:i-1) * y(1:i-1);
end

%后向替代求解 Ux = y
x = zeros(n, 1);
for i = n:-1:1
    x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
end

disp('方程组的解 x：');
disp(x);