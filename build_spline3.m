%1.b)
function [yp] = build_spline3(x, y, BcType)
% yy = build_spline3(x, y, BcType)
%
% 3��������ֵ��
% x��y�����ݵ㣨Ҫ��˳�����У�
% Bctype���߽�����
%         ��'period��  : ���ڱ߽�����
%         ��'natural�� : ��Ȼ�߽��������������˵�Ķ��׵���ֵΪ0
%         ��'endslope��: ָ�������˵��б�ʣ�һ�׵����������y�������������
% �������ظ����ڵ�ĵ���ֵyp�����ڹ���ÿ�������ϵ�3��Hermite��ֵ
%

n = length(x); % ��ֵ���ά��
m = n -1;      % �������
h = zeros(1,m); dy = zeros(1,m); alpha = zeros(1,m-1);

h = x(2:m+1) - x(1:m); % �����������ĳ���
dy= y(2:m+1) - y(1:m); % y�Ĳ��

alpha = h(1:m-1)./(h(1:m-1)+h(2:m)); % ����ϵ��������p39�� ��2.58��

switch BcType
    case 'natural' % ��Ȼ
        A = zeros(n,n); beta = zeros(n,1); yp = zeros(n,1);
        for i = 2:n-1
            A(i,i-1:i+1) = [1-alpha(i-1), 2, alpha(i-1)];
        end
        A(1, 1:2)   = [2, 1];
        A(n, n-1:n) = [1, 2];
        beta(2:n-1) = 3*((1-alpha).*dy(1:m-1)./h(1:m-1) + alpha.*dy(2:m)./h(2:m));
        beta(1)     = 3*dy(1)/h(1);
        beta(n)     = 3*dy(n-1)/h(n-1);
        yp = A\beta;
    case 'not_a_knot' % �����������ñ߽������Ľ���
        % ��ָû�нڵ��ϵı߽�����Ҫ�󣬵�Ҫ����β�ֽ�ڵ��ϵ����ε���һ��
        A = zeros(n,n); beta = zeros(n,1); yp = zeros(n,1);
        for i = 2:n-1
            A(i,i-1:i+1) = [1-alpha(i-1), 2, alpha(i-1)];
        end
        A(1, 1:3)   = [1, 1-(h(1)/h(2))^2, -(h(1)/h(2))^2];
        A(n, n-2:n) = [1, 1-(h(n-2)/h(n-1))^2, -(h(n-2)/h(n-1))^2];
        beta(2:n-1) = 3*((1-alpha).*dy(1:m-1)./h(1:m-1) + alpha.*dy(2:m)./h(2:m));
        beta(1)     = -2*y(1)/h(1) + 2*y(2)*(1/h(1)+h(1)^2/h(2)^3) -2*y(3)*h(1)^2/h(2)^3;
        beta(n)     = -2*y(n-2)/h(n-2) + 2*y(n-1)*(1/h(n-2)+h(n-2)^2/h(n-1)^3) -2*y(n)*h(n-2)^2/h(n-1)^3;
        yp = A\beta;
    case 'endslope'
        A = zeros(n-2,n-2); beta = zeros(n-2,1); yp = zeros(n,1);
        for i = 2:n-3
            A(i,i-1:i+1) = [1-alpha(i), 2, alpha(i)];
        end
        A(1, 1:2)   = [2, alpha(1)];
        A(n-2, n-3:n-2) = [1-alpha(n-2), 2];
%         beta(2:n-3) = 3*((1-alpha(2:n-3)).*dy(2:n-3)./h(2:n-3) + alpha(2:n-3).*dy(3:n-2)./h(3:n-2));
        beta(1:n-2) = 3*((1-alpha).*dy(1:m-1)./h(1:m-1) + alpha.*dy(2:m)./h(2:m));
        beta(1)     = beta(1) - (1-alpha(1))*y(n+1);
        beta(n-2)   = beta(n-2)- alpha(n-2)*y(n+2);
        yp(2:n-1) = A\beta; yp(1) = y(n+1); yp(n) = y(n+2);
    otherwise
        disp('�����߽�������δʵ�֣�')
end
end

