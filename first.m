%第一题
%1.a)下三角方程组函数
function x = LowerTri(L, b)
    n = length(b);      
    x = zeros(n, 1);    
    %从第一行开始逐步向下求解
    for k = 1:n
        if L(k, k) == 0
            error('矩阵L不可逆，存在零对角元素');
        end
        x(k) = (b(k) - L(k, 1:k-1) * x(1:k-1)) / L(k, k);
    end
end


%%
%1.a)上三角方程组函数
function x = UpperTri(U, b)
    n = length(b);    
    x = zeros(n, 1);   
    %从最后一行开始逐步向上求解
    for k = n:-1:1
        if U(k, k) == 0
            error('矩阵U不可逆，存在零对角元素');
        end
        x(k) = (b(k) - U(k, k+1:n) * x(k+1:n)) / U(k, k);
    end
end


%%
%1.a)P111例1
clear
clc
L = [1, 0, 0;
     2, 1, 0;
    -1, 2, 1];
b = [3; 1; -7];

x = LowerTri(L, b);
disp(['1a)例1解：']);
disp(x);


%%
%1.a)P111例2
clear
clc
U = [2, 2, 3;
     0, 3, 1;
     0, 0, 6];
b = [3; -5; 6];

x = UpperTri(U, b);
disp('1a)例2解：');
disp(x);


%%
%1.b)追赶法
function x = tridiagonal(a, b, c, d)
    n = length(b);  
    %前向消元，分解A为LU
    for i = 2:n
        factor = a(i-1) / b(i-1);
        b(i) = b(i) - factor * c(i-1);
        d(i) = d(i) - factor * d(i-1);
    end
    %回代求解，先解Ux=y
    x = zeros(n, 1);
    x(n) = d(n) / b(n); 
    for i = n-1:-1:1
        x(i) = (d(i) - c(i) * x(i+1)) / b(i);
    end
end


%%
%1.b)追赶法检验
clear
clc
n = 5;

a = [1; 2; 3; 4];         
b = [4; 5; 6; 7; 8];      
c = [1; 2; 3; 4];         
d = [7; 8; 15; 22; 30];  

x = tridiagonal(a, b, c, d);

disp('解向量 x：');
disp(x);

%构造原矩阵 A
A = diag(b) + diag(a, -1) + diag(c, 1);
disp('验证 Ax = d:');
disp(A*x);  %结果接近 d


%%
%1.b)以前的样条计算程序
function [yp] = build_spline3(x, y, BcType)
% yy = build_spline3(x, y, BcType)
%
% 3次样条插值：
% x，y：数据点（要求按顺序排列）
% Bctype：边界条件
%         ‘'period’  : 周期边界条件
%         ‘'natural’ : 自然边界条件，即两个端点的二阶导数值为0
%         ‘'endslope’: 指定两个端点的斜率（一阶导），存放在y的最后两个分量
% 函数返回各个节点的导数值yp，用于构造每个区间上的3次Hermite插值
%

n = length(x); % 插值点的维数
m = n -1;      % 区间个数
h = zeros(1,m); dy = zeros(1,m); alpha = zeros(1,m-1);

h = x(2:m+1) - x(1:m); % 计算各个区间的长度
dy= y(2:m+1) - y(1:m); % y的差分

alpha = h(1:m-1)./(h(1:m-1)+h(2:m)); % 计算系数，见书p39， （2.58）

switch BcType
    case 'natural' % 自然
        A = zeros(n,n); beta = zeros(n,1); yp = zeros(n,1);
        for i = 2:n-1
            A(i,i-1:i+1) = [1-alpha(i-1), 2, alpha(i-1)];
        end
        A(1, 1:2)   = [2, 1];
        A(n, n-1:n) = [1, 2];
        beta(2:n-1) = 3*((1-alpha).*dy(1:m-1)./h(1:m-1) + alpha.*dy(2:m)./h(2:m));
        beta(1)     = 3*dy(1)/h(1);
        beta(n)     = 3*dy(n-1)/h(n-1);
        yp = A\beta;
    case 'not_a_knot' % 请上网搜索该边界条件的解释
        % 意指没有节点上的边界条件要求，但要求首尾分界节点上的三次导数一致
        A = zeros(n,n); beta = zeros(n,1); yp = zeros(n,1);
        for i = 2:n-1
            A(i,i-1:i+1) = [1-alpha(i-1), 2, alpha(i-1)];
        end
        A(1, 1:3)   = [1, 1-(h(1)/h(2))^2, -(h(1)/h(2))^2];
        A(n, n-2:n) = [1, 1-(h(n-2)/h(n-1))^2, -(h(n-2)/h(n-1))^2];
        beta(2:n-1) = 3*((1-alpha).*dy(1:m-1)./h(1:m-1) + alpha.*dy(2:m)./h(2:m));
        beta(1)     = -2*y(1)/h(1) + 2*y(2)*(1/h(1)+h(1)^2/h(2)^3) -2*y(3)*h(1)^2/h(2)^3;
        beta(n)     = -2*y(n-2)/h(n-2) + 2*y(n-1)*(1/h(n-2)+h(n-2)^2/h(n-1)^3) -2*y(n)*h(n-2)^2/h(n-1)^3;
        yp = A\beta;
    case 'endslope'
        A = zeros(n-2,n-2); beta = zeros(n-2,1); yp = zeros(n,1);
        for i = 2:n-3
            A(i,i-1:i+1) = [1-alpha(i), 2, alpha(i)];
        end
        A(1, 1:2)   = [2, alpha(1)];
        A(n-2, n-3:n-2) = [1-alpha(n-2), 2];
%         beta(2:n-3) = 3*((1-alpha(2:n-3)).*dy(2:n-3)./h(2:n-3) + alpha(2:n-3).*dy(3:n-2)./h(3:n-2));
        beta(1:n-2) = 3*((1-alpha).*dy(1:m-1)./h(1:m-1) + alpha.*dy(2:m)./h(2:m));
        beta(1)     = beta(1) - (1-alpha(1))*y(n+1);
        beta(n-2)   = beta(n-2)- alpha(n-2)*y(n+2);
        yp(2:n-1) = A\beta; yp(1) = y(n+1); yp(n) = y(n+2);
    otherwise
        disp('其他边界条件尚未实现！')
end
end


%%
%1.b)以前的样条计算程序的检验
clear
clc
x = [75, 76, 77, 78, 79, 80];
y = [2.768, 2.833, 2.903, 2.979, 3.062, 3.153];

BcType = 'natural';
yp = build_spline3(x, y, BcType);
%定义插值点
x_query = [75.5, 78.3];

%构造样条插值函数
spline_value = zeros(size(x_query));
for k = 1:length(x_query)
    for i = 1:length(x)-1
        if x_query(k) >= x(i) && x_query(k) <= x(i+1)
            h = x(i+1) - x(i);
            t = (x_query(k) - x(i)) / h; % 归一化变量
            
            %Hermite插值公式
            h00 = 2*t^3 - 3*t^2 + 1;
            h10 = t^3 - 2*t^2 + t;
            h01 = -2*t^3 + 3*t^2;
            h11 = t^3 - t^2;
            
            spline_value(k) = h00*y(i) + h10*h*yp(i) + ...
                              h01*y(i+1) + h11*h*yp(i+1);
            break;
        end
    end
end

disp('插值点的函数值：');
disp(table(x_query', spline_value', 'VariableNames', {'x_query', 'SplineValue'}));


%%
%1.b)把函数代入以前的样条计算程序中
function [yp] = build_spline3_new(x, y, BcType)
    n = length(x); 
    m = n - 1;     

    h = diff(x);                
    dy = diff(y);               
    alpha = h(1:end-1) ./ (h(1:end-1) + h(2:end)); 
    beta = zeros(n, 1);         

    switch BcType
        case 'natural' 
            %三对角矩阵的系数
            a = [0; (1 - alpha(:))];        %下对角线
            b = [2; 2 * ones(n - 2, 1); 2]; %主对角线 
            c = [alpha(:); 0];              %上对角线

            %右端项
            beta(2:n-1) = 3 * ((1 - alpha) .* (dy(1:end-1) ./ h(1:end-1)) + ...
                               alpha .* (dy(2:end) ./ h(2:end)));
            beta(1) = 3 * dy(1) / h(1);
            beta(n) = 3 * dy(end) / h(end);

            %调用tridiagonal函数求解
            yp = tridiagonal(a, b, c, beta);

        otherwise
            disp('其他边界条件尚未实现！');
            yp = [];
    end
end


%%
%1.b)再进行一次检验
clear
clc
x = [75, 76, 77, 78, 79, 80];
y = [2.768, 2.833, 2.903, 2.979, 3.062, 3.153];

BcType = 'natural';

yp = build_spline3_new(x, y, BcType);

x_query = [75.5, 78.3];

spline_value = zeros(size(x_query));
for k = 1:length(x_query)
    for i = 1:length(x)-1
        if x_query(k) >= x(i) && x_query(k) <= x(i+1)
            h = x(i+1) - x(i);
            t = (x_query(k) - x(i)) / h;
           
            %Hermite插值
            h00 = 2*t^3 - 3*t^2 + 1;
            h10 = t^3 - 2*t^2 + t;
            h01 = -2*t^3 + 3*t^2;
            h11 = t^3 - t^2;
            
            spline_value(k) = h00*y(i) + h10*h*yp(i) + ...
                              h01*y(i+1) + h11*h*yp(i+1);
            break;
        end
    end
end

disp('插值点的函数值：');
disp(table(x_query', spline_value', 'VariableNames', {'x_query', 'SplineValue'}));
