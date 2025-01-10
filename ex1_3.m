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