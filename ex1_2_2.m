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