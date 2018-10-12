function[ellip] = ellipse(rx,ry,x,y,ang)
th = 0:pi/50:2*pi;
xe = rx*cos(th);
ye = ry*sin(th);
if ~isempty(ang)
    ellip = Rot([xe' ye'],ang);
else
    ellip = [xe' ye'];
end
ellip = [ellip(:,1)+x ellip(:,2)+y];
end