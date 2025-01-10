%1.b)
clear
clc
n = 5;

a = [1; 2; 3; 4];         
b = [4; 5; 6; 7; 8];      
c = [1; 2; 3; 4];         
d = [7; 8; 15; 22; 30];  

x = tridiagonal(a, b, c, d);

disp('解向量 x：');
disp(x);

%构造原矩阵 A
A = diag(b) + diag(a, -1) + diag(c, 1);
disp('验证 Ax = d:');
disp(A*x);  %结果接近 d