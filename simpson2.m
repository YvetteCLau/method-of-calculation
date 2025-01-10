clear;
clc;

f1 = @(x) 1./x;
f2 = @(x) 100.*sin(10./x)./(x.^2);
f3 = @(x) sin(x)./x;
f4 = @(x) exp(-x.^2);
intervals = {[1,100], [1,3], [eps,pi/2], [0,100]};
functions = {f1, f2, f3, f4};
func_names = {'1/x', '100sin(10/x)/x^2', 'sin(x)/x', 'exp(-x^2)'};

tol_values = logspace(-12, -3, 10);
results = zeros(length(tol_values), length(functions));
points_used = zeros(length(tol_values), length(functions));
quad_results = zeros(length(tol_values), length(functions));
quad_points = zeros(length(tol_values), length(functions));

%精确值
exact_values = zeros(1, length(functions));
for i = 1:length(functions)
    exact_values(i) = integral(functions{i}, intervals{i}(1), intervals{i}(2), 'AbsTol', 1e-15);
end

for i = 1:length(tol_values)
    tol = tol_values(i);
    for j = 1:length(functions)
        %Simpson
        [I_simpson, n_simpson] = Simpson_Integration2(functions{j}, intervals{j}(1), intervals{j}(2), tol);
        points_used(i,j) = n_simpson;
        results(i,j) = abs(I_simpson - exact_values(j));
        
        %quad
        [Q,fcnt] = quad(functions{j}, intervals{j}(1), intervals{j}(2), tol);
        quad_results(i,j) = abs(Q - exact_values(j));
        quad_points(i,j) = fcnt;
    end
end

%误差比较
figure
for j = 1:length(functions)
    loglog(tol_values, results(:,j), 'o-', 'LineWidth', 1.5);
    hold on;
    loglog(tol_values, quad_results(:,j), 's--', 'LineWidth', 1.5);
end
loglog(tol_values, tol_values, 'k--', 'LineWidth', 1);
hold off;
grid on;
title('误差比较');
xlabel('误差容限');
ylabel('实际误差');
legend([strcat('Simpson-', func_names), strcat('quad-', func_names), '理想误差线'], ...
    'Location', 'northwest','FontSize',7.5);

%所需分点数
figure
subplot(1,2,1)
for j = 1:length(functions)
    semilogx(tol_values, points_used(:,j), 'o-', 'LineWidth', 1.5);
    hold on;
end
hold off;
grid on;
title('Simpson方法所需分点数');
xlabel('误差容限');
ylabel('分点数');
legend(func_names, 'Location', 'best');

subplot(1,2,2)
for j = 1:length(functions)
    semilogx(tol_values, quad_points(:,j), 'o-', 'LineWidth', 1.5);
    hold on;
end
hold off;
grid on;
title('quad方法所需分点数');
xlabel('误差容限');
ylabel('分点数');
legend(func_names, 'Location', 'best');

%收敛率分析
figure
subplot(1,2,1)
for j = 1:length(functions)
    loglog(points_used(:,j), results(:,j),'o-','LineWidth', 1.5);
    hold on;
end
hold off;
grid on;
title('Simpson方法收敛率分析');
xlabel('分点数');
ylabel('实际误差');
legend(func_names, 'Location','best','FontSize',6);

subplot(1,2,2)
for j = 1:length(functions)
    loglog(quad_points(:,j), quad_results(:,j),'o-','LineWidth', 1.5);
    hold on;
end
hold off;
grid on;
title('quad方法收敛率分析');
xlabel('分点数');
ylabel('实际误差');
legend(func_names, 'Location','best');

%效率分析 (误差×分点数)
figure
subplot(1,2,1)
efficiency_simpson=results .* points_used;
for j = 1:length(functions)
    loglog(tol_values, efficiency_simpson(:,j), 'o-','LineWidth', 1.5);
    hold on;
end
hold off;
grid on;
title('Simpson方法效率分析 (误差×分点数)');
xlabel('误差容限');
ylabel('效率指标');
legend(func_names, 'Location', 'best');

subplot(1,2,2)
efficiency_quad=quad_results .* quad_points;
for j = 1:length(functions)
    loglog(tol_values, efficiency_quad(:,j), 'o-','LineWidth', 1.5);
    hold on;
end
hold off;
grid on;
title('quad方法效率分析 (误差×分点数)');
xlabel('误差容限');
ylabel('效率指标');
legend(func_names, 'Location', 'best');

%数值比较结果
fprintf('数值比较结果：\n');
fprintf('容许误差\t函数\t\tSimpson误差\tquad误差\tSimpson分点数\tquad分点数\n');
for i = 1:length(tol_values)
    for j = 1:length(functions)
        fprintf('%.2e\t%s\t%.2e\t%.2e\t%d\t\t%d\n', ...
            tol_values(i), func_names{j}, results(i,j), ...
            quad_results(i,j), points_used(i,j), quad_points(i,j));
    end
end

%收敛率
fprintf('\nSimpson方法收敛率：\n');
for j = 1:length(functions)
    p = polyfit(log(points_used(end-4:end,j)), ...
        log(results(end-4:end,j)), 1);
    fprintf('%s: %.2f\n', func_names{j}, -p(1));
end

fprintf('\nquad方法收敛率：\n');
for j = 1:length(functions)
    p = polyfit(log(quad_points(end-4:end,j)), ...
        log(quad_results(end-4:end,j)), 1);
    fprintf('%s: %.2f\n', func_names{j}, -p(1));
end