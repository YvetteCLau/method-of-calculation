% 复化梯形公式函数
function [Tn] = int_Trapezoidal(f, a, b, n)
    h = (b-a)/n;
    x = a: h: b;
    y = f(x);
    Tn = h * (sum(y) - (y(1) + y(end))/2);
end