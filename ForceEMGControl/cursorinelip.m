
cursorCir = circle(1,2,0);
targetCir = ellipse(2,1,0,0,[]);

figure
plot(targetCir(:,1),targetCir(:,2),'r','Linewidth',3);
hold on; 
plot(cursorCir(:,1),cursorCir(:,2),'g','Linewidth',3); 
axis square;

rcursor = sqrt((cursorCir(1,1)-cursorCir(round(length(cursorCir)/2),1))^2+(cursorCir(1,2)-cursorCir(round(length(cursorCir)/2),2))^2)/2;
rt1 = sqrt((targetCir(1,1)-targetCir(round(length(targetCir)/2),1))^2)/2;
rt2 = sqrt((targetCir(1,2)-targetCir(round(length(targetCir)/2),1))^2)/2;

ccursor = [cursorCir(1,1)-rcursor,cursorCir(1,2)];
ctarget = [targetCir(1,1)-rt1,targetCir(1,2)];

dvec = [ccursor(1)-ctarget(1) ccursor(2)-ctarget(2)];
ang = atan(dvec(2)/dvec(1));
reli = rt1*rt2/sqrt((rt1*cos(ang))^2+(rt2*sin(ang))^2);

d = sqrt((ccursor(1)-ctarget(1))^2+(ccursor(2)-ctarget(2))^2);

if d < rcursor+reli
    cursorin = true;
else
    cursorin = false;
end