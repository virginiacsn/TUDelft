function[cursorin] = cursorInTarget(cursorCir,targetCir)
rcursor = sqrt((cursorCir(1,1)-cursorCir(round(length(cursorCir)/2),1))^2+(cursorCir(1,2)-cursorCir(round(length(cursorCir)/2),2))^2)/2;
rtarget = sqrt((targetCir(1,1)-targetCir(round(length(targetCir)/2),1))^2+(targetCir(1,2)-targetCir(round(length(targetCir)/2),2))^2)/2;

ccursor = [cursorCir(1,1)-rcursor,cursorCir(1,2)];
ctarget = [targetCir(1,1)-rtarget,targetCir(1,2)];


d = sqrt((ccursor(1)-ctarget(1))^2+(ccursor(2)-ctarget(2))^2);

if d < rcursor+rtarget
    cursorin = true;
else
    cursorin = false;
end
end
