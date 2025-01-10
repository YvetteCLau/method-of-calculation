%1.a)
clear
clc
L = [1, 0, 0;
     2, 1, 0;
    -1, 2, 1];
b = [3; 1; -7];

x = LowerTri(L, b);
disp(['1a)例1解：']);
disp(x);