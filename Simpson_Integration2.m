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

