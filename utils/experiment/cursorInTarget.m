function[cursorin] = cursorInTarget(cursorCir,targetCir)

rcursor = sqrt((cursorCir(1,1)-cursorCir(round(length(cursorCir)/2),1))^2+(cursorCir(1,2)-cursorCir(round(length(cursorCir)/2),2))^2)/2;

rtargetx = sqrt((targetCir(1,1)-(targetCir(round(length(targetCir)/2),1)))^2+(targetCir(1,2)-(targetCir(round(length(targetCir)/2),2)))^2)/2;
rtargety = sqrt((targetCir(round(3*length(targetCir)/4),1)-(targetCir(round(length(targetCir)/4)+1,1)))^2+(targetCir(round(3*length(targetCir)/4),2)-(targetCir(round(length(targetCir)/4)+1,2)))^2)/2;

if rtargetx == rtargety
    rtarget = sqrt((targetCir(1,1)-targetCir(round(length(targetCir)/2),1))^2+(targetCir(1,2)-targetCir(round(length(targetCir)/2),2))^2)/2;
    
    ccursor = [cursorCir(1,1)-rcursor,cursorCir(1,2)];
    ctarget = [targetCir(1,1)-rtarget,targetCir(1,2)];
else
    ccursor = [cursorCir(1,1)-rcursor,cursorCir(1,2)];
    ctarget = [(targetCir(1,1)+(targetCir(round(length(targetCir)/2),1)))/2,(targetCir(round(3*length(targetCir)/4),2)+(targetCir(round(length(targetCir)/4)+1,2)))/2];
    
    dv = [ccursor(1)-ctarget(1) ccursor(2)-ctarget(2)];
    ange = atan((targetCir(1,2)-ctarget(2))/(targetCir(1,1)-ctarget(1)));
    angct = atan(dv(2)/dv(1));
    
    rtarget = (rtargetx*rtargety)/sqrt((rtargety*cos(angct-ange))^2+(rtargetx*sin(angct-ange))^2);
end

d = sqrt((ccursor(1)-ctarget(1))^2+(ccursor(2)-ctarget(2))^2);

if d < rcursor+rtarget
    cursorin = true;
else
    cursorin = false;
end
end
