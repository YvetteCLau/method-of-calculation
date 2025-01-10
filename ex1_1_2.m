%1.a)
clear
clc
U = [2, 2, 3;
     0, 3, 1;
     0, 0, 6];
b = [3; -5; 6];

x = UpperTri(U, b);
disp('1a)例2解：');
disp(x);