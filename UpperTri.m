%1.a)
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