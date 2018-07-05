function[cir] = circle(r,x,y)
th = 0:pi/50:2*pi;
xc = r*cos(th) + x;
yc = r*sin(th) + y;
cir = [xc' yc'];
end