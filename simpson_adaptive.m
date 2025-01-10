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

