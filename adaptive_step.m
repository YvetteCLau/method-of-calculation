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
