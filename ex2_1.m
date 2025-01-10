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