function[xyRot] = Rot(xy,ang)
if size(xy,1)~=2
    xytemp = xy';
else
    xytemp = xy;
end

xyRot = [cos(ang) -sin(ang); sin(ang) cos(ang)]*xytemp;

if size(xy,1)~=2
    xyRot = xyRot';
end
end