function[xyRot] = Rot(xy,ang)
if size(xy,2)>size(xy,1)
    xy = xy';
end

xyRot = [cosd(ang) sind(ang); -sind(ang) cosd(ang)]*xy;

if size(xy,2)>size(xy,1)
    xyRot = xyRot';
end
end