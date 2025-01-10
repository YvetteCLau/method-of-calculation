clear;
clc;

f1 = @(x) 1./x;
f2 = @(x) 100*sin(10./x)./x.^2;
f3 = @(x) sin(x)./x;
f4 = @(x) exp(-x.^2);

tol = 1e-6;

[I1, n1] = simpson_adaptive(f1, 1, 100, tol);
fprintf('积分1 结果: %.10f, 积分点的个数: %d\n', I1, n1);
[Q1, fcnt1] = quad(f1, 1, 100, tol);
fprintf('quad积分1 结果: %.10f, 积分点的个数: %d\n\n', Q1, fcnt1);

[I2, n2] = simpson_adaptive(f2, 1, 3, tol);
fprintf('积分2 结果: %.10f, 积分点的个数: %d\n', I2, n2);
[Q2, fcnt2] = quad(f2, 1, 3, tol);
fprintf('quad积分2 结果: %.10f, 积分点的个数: %d\n\n', Q2, fcnt2);

% 用eps避免0点
[I3, n3] = simpson_adaptive(f3, eps, pi/2, tol);  
fprintf('积分3 结果: %.10f, 积分点的个数: %d\n', I3, n3);
[Q3, fcnt3] = quad(f3, eps, pi/2, tol);
fprintf('quad积分3 结果: %.10f, 积分点的个数: %d\n\n', Q3, fcnt3);

[I4, n4] = simpson_adaptive(f4, 0, 100, tol);
fprintf('积分4 结果: %.10f, 积分点的个数 %d\n', I4, n4);
[Q4, fcnt4] = quad(f4, 0, 100, tol);
fprintf('quad积分4 结果: %.10f, 积分点的个数: %d\n', Q4, fcnt4);

% 验证积分4结果
theoretical4 = sqrt(pi)/2; %e^(-x^2)从0到无穷的值
fprintf('\n积分4的理论值: %.10f\n', theoretical4);
fprintf('相对误差: %.2e\n', abs(I4-theoretical4)/theoretical4);