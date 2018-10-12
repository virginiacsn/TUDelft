function[xyRot] = Rot(xy,ang)
if size(xy,2)<size(xy,1)
    xytemp = xy';
end

xyRot = [cos(ang) sin(ang); -sin(ang) cos(ang)]*xytemp;

if size(xy,2)<size(xy,1)
    xyRot = xyRot';
end
end