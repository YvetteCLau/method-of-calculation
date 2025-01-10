%第二题
%前期准备1 函数Simpson_Integration1
function [I] = Simpson_Integration1(f,a,b,n)
    if mod(n,2) ~= 0
        n = n + 1;
    end
    h = (b-a)/n;
    x = a:h:b;
    y = f(x);
    I = h/3 * (y(1)+y(end)+4*sum(y(2:2:end-1))+2*sum(y(3:2:end-2)));
end


%%
%前期准备1 simpson1.m：通过函数Simpson_Integration1得出的结果
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


%%
%前期准备2 函数adaptive_step
function [I, n] = adaptive_step(f,a,b, fa,fc, fb,S, tol,depth, maxIter)
    h = (b-a)/2;
    c = a + h;
    d = a + h/2;
    e = c + h/2;
    fd = f(d);
    fe = f(e);
    S1 = h/6*(fa + 4*fd + fc);
    S2 = h/6*(fc + 4*fe + fb);
    if depth>=maxIter
        I = S1 + S2;
        n = 2^depth + 1;
        return;
    end
    if abs(S1 + S2 - S) < 15*tol
        I = S1 + S2;
        n = 2^depth + 1;
    else
        [I1, n1] = adaptive_step(f,a,c, fa, fd, fc, S1, tol/2, depth+1, maxIter);
        [I2, n2] = adaptive_step(f,c,b, fc, fe, fb, S2, tol/2, depth+1, maxIter);
        I = I1 + I2;
        n = n1+n2-1;
    end
end


%%
%前期准备2 函数simpson_adaptive
function [I, n] = simpson_adaptive(f,a,b,tol)
    maxIter = 1000; 
    I = 0;
    n = 1;
    h = (b-a)/2;
    fa = f(a);
    fc = f(a+h);
    fb = f(b);
    S = h/3*(fa+4*fc+fb);
    [I, n] = adaptive_step(f,a,b, fa, fc, fb, S, tol, 1, maxIter);
end


%%
%前期准备2 simpson.m：通过函数adaptive_step和函数simpson_adaptive得出的结果。
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


%%
%完整结果 函数Simpson_Integration2
function [I,n] = Simpson_Integration2(f,a,b,tol)
    n = 10;
    err = 1;
    I_old = 0;  
    while err > tol
        h = (b-a)/n;
        x = a:h:b;
        y = f(x);
        I = h/3 * (y(1) + y(end) + 4*sum(y(2:2:end-1)) ...
            + 2*sum(y(3:2:end-2)));
        if n > 10
            err = abs(I-I_old);
        end
        I_old = I;
        n = n*2;
        if n>1e6
            warning('达到最大迭代次数');
            break;
        end
    end
    n = n/2;
end


%%
%完整结果 simpson2.m：通过函数Simpson_Integration2得出的结果
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