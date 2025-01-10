function [Sn, n, tol_c] = int_Simpson_subdivision(f, a, b, opts)
% 使用逐次分半法的 Simpson 公式计算积分，优化内存使用
%
% 输入：
%   f: 被积函数
%   a, b: 积分区间
%   opts.n: 初始区间数，可选，默认值=6
%   opts.tol: 容许误差，可选，默认值=1e-10
%
% 输出：
%   Sn: 积分结果
%   n: 最终的区间个数（偶数）
%   tol_c: 最终误差

% 参数验证和默认值
arguments
    f function_handle
    a double
    b double
    opts.n = 6; % 初始区间数
    opts.tol = 1e-10; % 容许误差
end

n = opts.n;
tol = opts.tol;

% 确保初始区间数为偶数
if mod(n, 2) ~= 0
    n = n + 1;
end

% 初始步长
h = (b - a) / n;

% 初始 Simpson 积分值
x1 = a;
x2 = a + h;
x3 = a + 2 * h;
Sn = 0;

% 分段计算初始积分
for i = 1:(n / 2)
    Sn = Sn + (h / 6) * (f(x1) + 4 * f(x2) + f(x3));
    x1 = x1 + 2 * h;
    x2 = x2 + 2 * h;
    x3 = x3 + 2 * h;
end

% 逐次分半法迭代
while true
    % 保存旧的积分值
    Sn_old = Sn;

    % 更新步长和区间数
    h = h / 2;
    n = n * 2;

    % 分段计算新增点的 Simpson 积分
    x1 = a;
    x2 = a + h;
    x3 = a + 2 * h;
    Sn = 0;

    for i = 1:(n / 2)
        Sn = Sn + (h / 6) * (f(x1) + 4 * f(x2) + f(x3));
        x1 = x1 + 2 * h;
        x2 = x2 + 2 * h;
        x3 = x3 + 2 * h;
    end

    % 计算误差
    tol_c = abs(Sn - Sn_old);

    % 检查误差是否满足容许要求
    if tol_c < tol
        break;
    end
end

end