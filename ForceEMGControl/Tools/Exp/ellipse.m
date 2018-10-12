function[ellip] = ellipse(rx,ry,x,y)
th = 0:pi/50:2*pi;
xe = rx*cos(th) + x;
ye = ry*sin(th) + y;
ellip = [xe' ye'];
end