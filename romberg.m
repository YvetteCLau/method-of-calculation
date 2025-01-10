clear;
clc;

f1 = @(x) 1./x;
f2 = @(x) 100.*sin(10./x)./(x.^2);
f3 = @(x) sin(x)./x;
intervals = {[1,100], [1,3], [eps,pi/2]};
functions = {f1, f2, f3};
func_names = {'1/x', '100sin(10/x)/x^2', 'sin(x)/x'};

tol_values = logspace(-8,-3,6);
rom_results = zeros(length(tol_values), length(functions));
rom_points = zeros(length(tol_values), length(functions));
simpson_results = zeros(length(tol_values), length(functions));
simpson_points = zeros(length(tol_values), length(functions));

%精确值
exact_values = zeros(1, length(functions));
for i = 1:length(functions)
    exact_values(i) = integral(functions{i},intervals{i}(1),...
        intervals{i}(2), 'AbsTol', 1e-15);
end

for i = 1:length(tol_values)
    tol = tol_values(i);
    for j = 1:length(functions)
        % Romberg
        [I_rom, n_rom]=Romberg_Integration(functions{j}, ...
            intervals{j}(1), intervals{j}(2), tol);
        rom_results(i,j) = abs(I_rom - exact_values(j));
        rom_points(i,j) = n_rom;
        
        % Simpson
        [I_simpson, n_simpson]=Simpson_Integration2(functions{j}, ...
            intervals{j}(1), intervals{j}(2), tol);
        simpson_results(i,j) = abs(I_simpson - exact_values(j));
        simpson_points(i,j) = n_simpson;
    end
end

% 误差比较
figure('Position', [100,100,600,500]);
for j = 1:length(functions)
    loglog(tol_values, rom_results(:,j), 'o-','LineWidth', 1.5);
    hold on;
    loglog(tol_values, simpson_results(:,j), 's--','LineWidth', 1.5);
end
loglog(tol_values, tol_values, 'k--', 'LineWidth', 1);
hold off;
grid on;
title('误差比较');
xlabel('误差容限');
ylabel('实际误差');
legend([strcat('Romberg-', func_names), ...
    strcat('Simpson-', func_names), '理想误差线'], 'Location', 'best');

% 所需分点数比较
figure('Position', [100, 100, 600, 500]);
for j = 1:length(functions)
    semilogx(tol_values,rom_points(:,j), 'o-', 'LineWidth', 1.5);
    hold on;
    semilogx(tol_values,simpson_points(:,j), 's--', 'LineWidth', 1.5);
end
hold off;
grid on;
title('所需分点数比较');
xlabel('误差容限');
ylabel('分点数');
legend([strcat('Romberg-', func_names),...
    strcat('Simpson-', func_names)],'Location', 'best');

% 收敛率分析
figure('Position', [100, 100, 600, 500]);
for j = 1:length(functions)
    loglog(rom_points(:,j), rom_results(:,j), 'o-', 'LineWidth', 1.5);
    hold on;
    loglog(simpson_points(:,j), ...
        simpson_results(:,j), 's--', 'LineWidth', 1.5);
end
hold off;
grid on;
title('收敛率分析');
xlabel('分点数');
ylabel('实际误差');
legend([strcat('Romberg-', func_names), ...
    strcat('Simpson-', func_names)],'Location', 'southwest','FontSize',4);

% 计算效率分析 (误差×分点数)
figure('Position', [100, 100, 600, 500]);
rom_efficiency = rom_results .* rom_points;
simpson_efficiency = simpson_results .* simpson_points;
for j = 1:length(functions)
    loglog(tol_values, rom_efficiency(:,j), 'o-', 'LineWidth', 1.5);
    hold on;
    loglog(tol_values, simpson_efficiency(:,j), 's--', 'LineWidth', 1.5);
end
hold off;
grid on;
title('计算效率分析 (误差×分点数)');
xlabel('误差容限');
ylabel('效率指标');
legend([strcat('Romberg-', func_names), ...
    strcat('Simpson-', func_names)], 'Location', 'best');

% 数值比较结果
fprintf('数值比较结果：\n');
fprintf('容许误差\t函数\t\tRomberg误差\tSimpson误差\tRomberg点数\tSimpson点数\n');
for i = 1:length(tol_values)
    for j = 1:length(functions)
        fprintf('%.2e\t%s\t%.2e\t%.2e\t%d\t\t%d\n', ...
            tol_values(i), func_names{j}, rom_results(i,j), ...
            simpson_results(i,j), rom_points(i,j), simpson_points(i,j));
    end
end

% 计算收敛率
fprintf('\n收敛率分析：\n');
for j = 1:length(functions)
    % Romberg
    x = log(rom_points(:,j));
    y = log(rom_results(:,j));
    x_mean = mean(x);
    x_std = std(x);
    y_mean = mean(y);
    y_std = std(y);
    x_norm = (x - x_mean) / x_std;
    y_norm = (y - y_mean) / y_std;
    p = polyfit(x_norm, y_norm, 1);
    rom_slope = p(1) * (y_std / x_std);
    
    % Simpson
    x = log(simpson_points(:,j));
    y = log(simpson_results(:,j));
    x_mean = mean(x);
    x_std = std(x);
    y_mean = mean(y);
    y_std = std(y);
    x_norm = (x - x_mean) / x_std;
    y_norm = (y - y_mean) / y_std;
    p = polyfit(x_norm, y_norm, 1);
    simpson_slope = p(1) * (y_std / x_std);
    
    fprintf('%s:\n',func_names{j});
    fprintf('Romberg收敛率: %.2f\n',-rom_slope);
    fprintf('Simpson收敛率: %.2f\n\n',-simpson_slope);
end

%Romberg过程展示
fprintf('\nRomberg积分过程展示:\n');
for func_idx = 1:length(functions)
    fprintf('\n函数 %s 在区间 [%.2f, %.2f] 上的Romberg积分过程:\n', ...
        func_names{func_idx}, intervals{func_idx}(1), ...
        intervals{func_idx}(2));
    
    tol = 1e-6;
    a = intervals{func_idx}(1);
    b = intervals{func_idx}(2);
    f = functions{func_idx};
    max_iter = 20;
    R = zeros(max_iter);
    h = b-a;
    R(1,1) = h/2*(f(a)+f(b));
    %Romberg迭代过程
    for i = 2:max_iter
        h = h/2;
        sum = 0;
        for k = 1:2^(i-2)
            sum = sum + f(a + (2*k-1)*h);
        end
        R(i,1) = 0.5*R(i-1,1) + h*sum;       
        %Richardson外推
        for j = 2:i
            R(i,j) = R(i,j-1) + (R(i,j-1) - R(i-1,j-1))/(4^(j-1) - 1);
        end      
        % 检查收敛性
        if i>1 && abs(R(i,i)-R(i-1,i-1))<tol
            break;
        end
    end
    %Romberg表
    fprintf('\nRomberg表:\n');
    for row = 1:i
        for col = 1:row
            if abs(R(row,col))>1e-3
                fprintf('%10.6f ',R(row,col));
            else
                fprintf('%10.2e ',R(row,col));
            end
        end
        fprintf('\n');
    end
    fprintf('\n最终结果: %.10e\n',R(i,i));
    fprintf('迭代次数: %d\n',i);
    fprintf('估计误差: %.2e\n',abs(R(i,i)-R(i-1,i-1)));
end