%1.b)
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