function [x,y] = Sigmacircle(cx,cy,P,N,N2)
if nargin < 4
    N=1;
    N2 = 20;
elseif nargin < 5
    N2 = 20;            % 画椭圆的采样点
end
[a,b,phi] = calculateEllipara(P);
Rot = [cos(phi),-sin(phi);
    sin(phi),cos(phi)];

theta = 0:pi/N2:2*pi;

xy = [N*a*cos(theta); N*b*sin(theta)];           % 椭圆上的边界点
xy = Rot*xy +[cx*ones(1,length(theta));cy*ones(1,length(theta))];      % 相对于质心的边界点

x = xy(1,:)';
y = xy(2,:)';

function [a,b,phi] = calculateEllipara(A)
% a,b 长短轴
% phi 旋转角
E = eig(A);
a = max(sqrt(E));
b = min(sqrt(E));
if a-b < 1e-5
    phi =0;
else
    phi = real(asin(2*A(1,2)/(a^2-b^2)))/2+2*pi;
end