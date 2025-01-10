%1.a)
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