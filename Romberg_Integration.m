function [I, n] = Romberg_Integration(f,a,b,tol)
    maxiter = 20;
    R = zeros(maxiter,maxiter);
    h = b - a;
    R(1,1) = h/2*(f(a)+f(b));
    for i = 2:maxiter
        h = h/2;
        x = a+h:2*h:b-h;
        R(i,1) = R(i-1,1)/2 + h*sum(f(x));
        for j = 2:i
            R(i,j) = R(i,j-1) + (R(i,j-1)-R(i-1,j-1))/(4^(j-1) - 1);
        end
        if i > 2 && abs(R(i,i)-R(i-1,i-1))<tol
            I = R(i,i);
            n = 2^(i-1)+1; 
            return;
        end
    end
    warning('达到最大迭代次数');
    I = R(maxiter,maxiter);
    n = 2^(maxiter-1) + 1;
end
