

theta = 0:0.01:2*pi;
xl = -800:1:800;
%circle
l1 = 350;
l2 = 300;
b = 200;
r = sqrt(l1^2 +l2^2);
x = r*cos(theta) + b;
y = r*sin(theta);

%line
th1 = 45 * pi/180;
th2 = -20 * pi/180;

SH = [b 0];
EL = SH + [sin(th1)*l2, cos(th1)*l2];
H = EL + [sin(th1 + th2)*l1, cos(th1 + th2)*l1]; 

figure(1)
plot([SH(1) EL(1)],[SH(2) EL(2)],'b'), hold on
plot([EL(1) H(1)],[EL(2) H(2)],'b'),
plot(H(1),H(2),'r','markersize',10),

plot([0 H(1)],[0 H(2)],'m')

% y = x

t = H(2) / H(1);

yl = t*xl;
% for i = 1:length(x)
% y(i) = sqrt(x(i)^2 - r^2)
% end

figure(1)
plot(xl,yl,'c--')
plot(x,y)
axis square

syms x y
% r = sqrt(l1^2 +l2^2);
% t = H(2) / H(1);
a = (x - b)^2 + (t*x)^2 - r^2 
% - 3*x - y

solve(a == 0)

