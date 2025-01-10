%1.b)
clear
clc
x = [75, 76, 77, 78, 79, 80];
y = [2.768, 2.833, 2.903, 2.979, 3.062, 3.153];

BcType = 'natural';

yp = build_spline3_new(x, y, BcType);

x_query = [75.5, 78.3];

spline_value = zeros(size(x_query));
for k = 1:length(x_query)
    for i = 1:length(x)-1
        if x_query(k) >= x(i) && x_query(k) <= x(i+1)
            h = x(i+1) - x(i);
            t = (x_query(k) - x(i)) / h;
           
            %Hermite插值
            h00 = 2*t^3 - 3*t^2 + 1;
            h10 = t^3 - 2*t^2 + t;
            h01 = -2*t^3 + 3*t^2;
            h11 = t^3 - t^2;
            
            spline_value(k) = h00*y(i) + h10*h*yp(i) + ...
                              h01*y(i+1) + h11*h*yp(i+1);
            break;
        end
    end
end

disp('插值点的函数值：');
disp(table(x_query', spline_value', 'VariableNames', {'x_query', 'SplineValue'}));