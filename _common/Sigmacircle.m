function [x,y] = Sigmacircle(cx,cy,P,N,N2)
if nargin < 4
    N=1;
    N2 = 20;
elseif nargin < 5
    N2 = 20;            % ����Բ�Ĳ�����
end
[a,b,phi] = calculateEllipara(P);
Rot = [cos(phi),-sin(phi);
    sin(phi),cos(phi)];

theta = 0:pi/N2:2*pi;

xy = [N*a*cos(theta); N*b*sin(theta)];           % ��Բ�ϵı߽��
xy = Rot*xy +[cx*ones(1,length(theta));cy*ones(1,length(theta))];      % ��������ĵı߽��

x = xy(1,:)';
y = xy(2,:)';

function [a,b,phi] = calculateEllipara(A)
% a,b ������
% phi ��ת��
E = eig(A);
a = max(sqrt(E));
b = min(sqrt(E));
if a-b < 1e-5
    phi =0;
else
    phi = real(asin(2*A(1,2)/(a^2-b^2)))/2+2*pi;
end