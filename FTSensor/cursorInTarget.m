function[cursorin] = cursorInTarget(cursorCir,targetCir)
areaoverlap = polyarea(cursorCir(:,1)-targetCir(:,1),cursorCir(:,2)-targetCir(:,2));
areatarget = polyarea(targetCir(:,1),targetCir(:,2));
if areaoverlap/areatarget >= 0.9
    cursorin = true;
else
    cursorin = false;
end
end
