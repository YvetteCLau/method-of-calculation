%第二题
%2函数
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


%%
%2函数检验 p122例6
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


%%
%2.a)
clear
clc

A = [51, 49, 23, 4, 68, 71, 74;
     9, 58, 49, 88, 14, 90, 50;
     26, 24, 62, 91, 72, 89, 48;
     80, 46, 68, 79, 11, 34, 90;
     3, 96, 40, 10, 65, 70, 61;
     92, 55, 37, 26, 49, 20, 62;
     73, 52, 98, 34, 78, 4, 86];

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
%输出
disp('自定义函数分解的置换矩阵 P：');
disp(P);
disp('自定义函数分解的下三角矩阵 L：');
disp(L);
disp('自定义函数分解的上三角矩阵 U：');
disp(U);

disp('验证 PA = LU（自定义函数）：');
disp(P * A - L * U);  %结果应为零矩阵

%使用 MATLAB 自带的 lu 函数进行分解
[L_matlab, U_matlab, P_matlab] = lu(A);
%输出
disp('MATLAB 自带的置换矩阵 P：');
disp(P_matlab);
disp('MATLAB 自带的下三角矩阵 L：');
disp(L_matlab);
disp('MATLAB 自带的上三角矩阵 U：');
disp(U_matlab);

disp('验证 PA = LU（MATLAB 自带函数）：');
disp(P_matlab * A - L_matlab * U_matlab);  %结果应为零矩阵

%比较自定义函数和 MATLAB 的分解结果
disp('比较置换矩阵 P：');
disp(P - P_matlab);
disp('比较下三角矩阵 L：');
disp(L - L_matlab);
disp('比较上三角矩阵 U：');
disp(U - U_matlab);


%%
%2.b)
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

[L, U, P] = lu(A);

%求解 Ax = b1
Pb1 = P * b1;
y1 = LowerTri(L, Pb1);
x1 = UpperTri(U, y1);

%求解 Ax = b2
Pb2 = P * b2;
y2 = LowerTri(L, Pb2);
x2 = UpperTri(U, y2);

%输出
disp('方程组 Ax = b1 的解 x1：');
disp(x1);

disp('方程组 Ax = b2 的解 x2：');
disp(x2);