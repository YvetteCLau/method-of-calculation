%第一题
%复化梯形公式函数
function [Tn] = int_Trapezoidal(f, a, b, n)
    h = (b-a)/n;
    x = a: h: b;
    y = f(x);
    Tn = h * (sum(y) - (y(1) + y(end))/2);
end


%%
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


%%
%1.b)exp(x)
f = @(x) exp(x);    
a = 0;              
b = 1;              
m = 6;             
n = 7;              

err = zeros(1, n);  
h = zeros(1, n);    

for k = 1:n
    h(k) = (b-a) / (m * 2^k);  
    err(k) = abs(integral(f, a, b) - int_Trapezoidal(f, a, b, m * 2^k));  
end

figure;
loglog(h, err, 'r-*', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;

loglog(h, h.^2, 'b--', 'LineWidth', 1.5); 

p = polyfit(log(h), log(err), 1);  
loglog(h, exp(polyval(p, log(h))), 'g-.', 'LineWidth', 1.5); 

grid on;  
title('非周期函数收敛阶', 'FontSize', 14); 
xlabel('步长 h (log scale)', 'FontSize', 12); 
ylabel('误差 |E(h)| (log scale)', 'FontSize', 12); 
legend('误差', '理论收敛阶 (h^2)', sprintf('拟合误差 y = %.2fh^{%.2f}',...
    exp(p(2)), p(1)), 'Location', 'SouthEast', 'FontSize', 10); 
hold off;


%%
%1.b)sin(x).^2 .* x.^2
f = @(x) sin(x).^2 .* x.^2;  
a = -pi;                    
b = pi;                    
m = 6;                      
n = 7;                     

err = zeros(1, n);         
h = zeros(1, n);         

for k = 1:n
    h(k) = (b - a) / (m * 2^k);  
    err(k) = abs(integral(f, a, b) - int_Trapezoidal(f, a, b, m * 2^k)); 
end

figure;
loglog(h, err, 'r-*', 'LineWidth', 1.5, 'MarkerSize', 8); 
hold on;

loglog(h, h.^4, 'b--', 'LineWidth', 1.5);

p = polyfit(log(h), log(err), 1);  
loglog(h, exp(polyval(p, log(h))), 'g-.', 'LineWidth', 1.5); 

grid on; 
title('周期函数收敛阶', 'FontSize', 14);
xlabel('步长 h (log scale)', 'FontSize', 12); 
ylabel('误差 |E(h)| (log scale)', 'FontSize', 12); 
legend('误差', '理论收敛阶 (h^4)', ...
       sprintf('拟合误差 y = %.2fh^{%.0f}', exp(p(2)), p(1)), ...
       'Location', 'SouthEast', 'FontSize', 10); 
hold off;


%%
%1.c)
% f = @(x) sin(x); % 如果测试函数用sin(x)，误差太小，很快就达到机器精度，不能看出收敛速度
f = @(x) exp(cos(x)); 
a = 0; b = 2*pi;      
m = 1;              
n = 7;               

err = zeros(1, n);    
h = zeros(1, n);     

for k = 1:n
    h(k) = (b - a) / (m + 2*k); 
    err(k) = abs(integral(f, a, b) - int_Trapezoidal(f, a, b, m + 2*k));
end

figure;
loglog(h, err, 'r-*', 'LineWidth', 1.5, 'MarkerSize', 8); 
hold on;

loglog(h, h.^12, 'b--', 'LineWidth', 1.5); 

grid on; 
title('解析的周期函数收敛阶', 'FontSize', 14);
xlabel('步长 h (log scale)', 'FontSize', 12); 
ylabel('误差 |E(h)| (log scale)', 'FontSize', 12);
legend('误差', '理论收敛阶 (h^{12})', ...
       'Location', 'southeast', 'FontSize', 10);
hold off;

figure;
semilogy(m + 3*[1:n], err, '-b*', 'LineWidth', 1.5, 'MarkerSize', 8); 
grid on;
title('误差与插值节点数的关系', 'FontSize', 14); 
xlabel('插值节点数：n', 'FontSize', 12); 
ylabel('误差', 'FontSize', 12); 
legend('误差', 'FontSize', 10); 
hold off;