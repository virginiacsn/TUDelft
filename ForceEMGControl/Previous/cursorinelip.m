
cursorCir = circle(1,1,1);
targetCir = ellipse(2,1,0,0,0);

figure
plot(targetCir(:,1),targetCir(:,2),'r','Linewidth',3);
axis equal
hold on; grid on;
plot(cursorCir(:,1),cursorCir(:,2),'g','Linewidth',3); 
plot(targetCir(1,1),targetCir(1,2),'*')
plot(targetCir(round(length(targetCir)/2),1),targetCir(round(length(targetCir)/2),2),'co')
plot(targetCir(round(length(targetCir)/4)+1,1),targetCir(round(length(targetCir)/4)+1,2),'o')
xlim([-4 4]);ylim([-4 4]);

rcursor = sqrt((cursorCir(1,1)-cursorCir(round(length(cursorCir)/2),1))^2+(cursorCir(1,2)-cursorCir(round(length(cursorCir)/2),2))^2)/2;

rtx = sqrt((targetCir(1,1)-(targetCir(round(length(targetCir)/2),1)))^2+(targetCir(1,2)-(targetCir(round(length(targetCir)/2),2)))^2)/2;
rty = sqrt((targetCir(round(3*length(targetCir)/4),1)-(targetCir(round(length(targetCir)/4)+1,1)))^2+(targetCir(round(3*length(targetCir)/4),2)-(targetCir(round(length(targetCir)/4)+1,2)))^2)/2;

line([targetCir(1,1) targetCir(round(length(cursorCir)/2),1)],[targetCir(1,2) targetCir(round(length(cursorCir)/2),2)])
line([targetCir(round(3*length(targetCir)/4),1) targetCir(round(length(targetCir)/4)+1,1)],[targetCir(round(3*length(targetCir)/4),2) targetCir(round(length(targetCir)/4)+1,2)])

ccursor = [cursorCir(1,1)-rcursor,cursorCir(1,2)];
ctarget = [(targetCir(1,1)+(targetCir(round(length(targetCir)/2),1)))/2,(targetCir(round(3*length(targetCir)/4),2)+(targetCir(round(length(targetCir)/4)+1,2)))/2];

plot(ccursor(1),ccursor(2),'g*')
plot(ctarget(1),ctarget(2),'r*')

dvec = [ccursor(1)-ctarget(1) ccursor(2)-ctarget(2)];
ange = atan((targetCir(1,2)-ctarget(2))/(targetCir(1,1)-ctarget(1)));
angct = atan(dvec(2)/dvec(1));
angd = rad2deg(angct)+rad2deg(ange);

%  if (ange == 0) 
    reli = (rtx*rty)/sqrt((rty*cos(angct-ange))^2+(rtx*sin(angct-ange))^2);
% else
%   reli = (rtx*rty)/sqrt((rtx*cos(angct))^2+(rty*sin(angct))^2);
% end

d = sqrt((ccursor(1)-ctarget(1))^2+(ccursor(2)-ctarget(2))^2);

if d < rcursor+reli
    cursorin = true;
else
    cursorin = false;
end