%1.b)
function [yp] = build_spline3_new(x, y, BcType)
    n = length(x); 
    m = n - 1;     

    h = diff(x);                
    dy = diff(y);               
    alpha = h(1:end-1) ./ (h(1:end-1) + h(2:end)); 
    beta = zeros(n, 1);         

    switch BcType
        case 'natural' 
            %三对角矩阵的系数
            a = [0; (1 - alpha(:))];        %下对角线
            b = [2; 2 * ones(n - 2, 1); 2]; %主对角线 
            c = [alpha(:); 0];              %上对角线

            %右端项
            beta(2:n-1) = 3 * ((1 - alpha) .* (dy(1:end-1) ./ h(1:end-1)) + ...
                               alpha .* (dy(2:end) ./ h(2:end)));
            beta(1) = 3 * dy(1) / h(1);
            beta(n) = 3 * dy(end) / h(end);

            %调用tridiagonal函数求解
            yp = tridiagonal(a, b, c, beta);

        otherwise
            disp('其他边界条件尚未实现！');
            yp = [];
    end
end