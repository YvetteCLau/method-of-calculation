clear;
clc;

f1 = @(x) 1./x;
f2 = @(x) 100.*sin(10./x)./(x.^2);
f3 = @(x) sin(x)./x;
f4 = @(x) exp(-x.^2);
n = 1000;  
tol = 1e-6; 

%四个函数的图像
figure('Position', [100, 100, 1200, 800]);

subplot(2,2,1);
x1 = linspace(1, 100, 1000);
plot(x1, f1(x1), 'b-', 'LineWidth', 1.5);
title('f(x) = 1/x');
xlabel('x'); ylabel('f(x)');
grid on;

subplot(2,2,2);
x2 = linspace(1, 3, 1000);
plot(x2, f2(x2), 'r-', 'LineWidth', 1.5);
title('f(x) = 100sin(10/x)/x^2');
xlabel('x'); ylabel('f(x)');
grid on;

subplot(2,2,3);
x3 = linspace(0.01, pi/2, 1000);
plot(x3, f3(x3), 'g-', 'LineWidth', 1.5);
title('f(x) = sin(x)/x');
xlabel('x'); ylabel('f(x)');
grid on;

subplot(2,2,4);
x4 = linspace(0, 5, 1000); 
plot(x4, f4(x4), 'm-', 'LineWidth', 1.5);
title('f(x) = e^{-x^2}');
xlabel('x'); ylabel('f(x)');
grid on;

%分析图像
figure('Position', [100, 100, 1200, 400]);
n_values = [10 20 40 80 160 320 640];
errors = zeros(length(n_values), 4);
%计算各个函数在不同n值下的误差
for i = 1:length(n_values)
    I1 = Simpson_Integration1(f1,1,100,n_values(i));
    I2 = Simpson_Integration1(f2,1,3,n_values(i));
    I3 = Simpson_Integration1(f3,eps,pi/2,n_values(i));
    I4 = Simpson_Integration1(f4,0,100,n_values(i));    
    I1_exact = integral(f1, 1, 100);
    I2_exact = integral(f2, 1, 3);
    I3_exact = integral(f3, eps, pi/2);
    I4_exact = integral(f4, 0, 100); 
    errors(i,1) = abs(I1 - I1_exact);
    errors(i,2) = abs(I2 - I2_exact);
    errors(i,3) = abs(I3 - I3_exact);
    errors(i,4) = abs(I4 - I4_exact);
end
%收敛性分析
subplot(1,2,1);
loglog(n_values, errors, 'o-', 'LineWidth', 1.5);
legend('1/x', '100sin(10/x)/x^2', 'sin(x)/x', 'exp(-x^2)', 'Location', 'southwest');
title('收敛性分析');
xlabel('分段数 n'); ylabel('误差');
grid on;
%收敛阶分析
conv_order = log2(errors(1:end-1,:)./errors(2:end,:));
subplot(1,2,2);
plot(2:length(n_values), conv_order, 'o-', 'LineWidth', 1.5);
legend('1/x', '100sin(10/x)/x^2', 'sin(x)/x', 'exp(-x^2)', 'Location', 'best');
title('收敛阶分析');
xlabel('细化次数'); ylabel('收敛阶');
grid on;

%数值积分结果
fprintf('数值积分结果：\n');
fprintf('积分1（1/x）: %.10f\n', Simpson_Integration1(f1, 1, 100, n));
fprintf('积分2（100sin(10/x)/x^2）: %.10f\n', Simpson_Integration1(f2, 1, 3, n));
fprintf('积分3（sin(x)/x）: %.10f\n', Simpson_Integration1(f3, eps, pi/2, n));
fprintf('积分4（e^(-x^2)）: %.10f\n', Simpson_Integration1(f4, 0, 100, n));