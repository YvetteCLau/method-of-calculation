%1.a)
clear
clc

f1 = @(x) exp(x);
a = 0; b = 1;n = 100;
result1 = int_Trapezoidal(f1,a,b,n);
exact1 = exp(1) - 1;
fprintf('测试exp(x)在[0,1]上的积分:\n');
fprintf('复化梯形公式的结果: %.10f\n', result1);
fprintf('精确值: %.10f\n', exact1);
fprintf('误差: %.10e\n\n', abs(result1-exact1));

%可视化积分过程
figure
x_plot = linspace(a,b,1000);
y_plot = f1(x_plot);
x_trap = linspace(a,b,n+1);
y_trap = f1(x_trap);
plot(x_plot, y_plot, 'b-', 'LineWidth', 1.5);
hold on;
plot(x_trap, y_trap, 'ro-');
fill([x_trap fliplr(x_trap)], [y_trap zeros(size(y_trap))], ...
    'g', 'FaceAlpha', 0.3);
title('exp(x)的复化梯形积分示意图');
xlabel('x'); ylabel('y');
grid on;
