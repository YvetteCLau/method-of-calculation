function [I] = Simpson_Integration1(f,a,b,n)
    if mod(n,2) ~= 0
        n = n + 1;
    end
    h = (b-a)/n;
    x = a:h:b;
    y = f(x);
    I = h/3 * (y(1)+y(end)+4*sum(y(2:2:end-1))+2*sum(y(3:2:end-2)));
end
