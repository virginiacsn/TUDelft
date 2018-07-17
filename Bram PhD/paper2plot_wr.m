close all
clear all

datavec1 = {'PP1' 'PP2' 'PP3' 'PP4' 'PP5' 'PP6' 'PP7' 'PP8' 'PP9' 'PP10' 'PP11' 'PP12'};
datavec = {'PP1' 'PP2' 'PP3' 'PP4' 'PP5' 'PP6' 'PP7' 'PP8' 'PP9' 'PP10' 'PP11' 'PP12'};
cond = {'_p1' '_p2' '_p3' '_p4'};
fl = 20;
sc = 7.5;
plxp = [0:1:360] * pi/180;
xan = [1:1:360];
% for f = 1:3
% for f = 3
for q =1:4
    
for p = 1:12
% loading data

    txt = ['D:\matlab\Data\Exp2\data\' char(datavec(p)) char(cond(q))];
  
    load(txt)
   
for i = 1:128
    
    FX(i,:) = mean(dataout.data(1,i).FX);
    FY(i,:) = mean(dataout.data(1,i).FY);
    FZ(i,:) = mean(dataout.data(1,i).FZ);
    
end

%Differences in angle between target and reproduction
hoek = dataout.protocol.krachthoek *180/pi;

for i = 1:128
    
    if hoek(i) == 180 && atan2d(FY(i),FX(i)) < 0
        
        ang(i,1) = hoek(i) + atan2d(FY(i),FX(i));
        
    else
        
        ang(i,1) = atan2d(FY(i),FX(i)) - hoek(i);
    end
    
%     ang2(i,1) = atan2d(FY(i),FX(i)) - hoek(i);
end 

%All data in one vector
FA = [hoek FX FY FZ ang];
FAv = FA(1:2:end,:);
FAb = FA(2:2:end,:);

%Calculating force sqrt(x^2 + y^2) (vector length)
for e = 1:64
    
    FAv(e,6) = sqrt(FAv(e,2)^2 + FAv(e,3)^2);
    FAb(e,6) = sqrt(FAb(e,2)^2 + FAb(e,3)^2);
    
    %and in 3d, we dont use this in the analysis
    Fv3(e) = sqrt(sqrt(FAv(e,2)^2 + FAv(e,3)^2)^2 + FAv(e,4)^2);
    Fb3(e) = sqrt(sqrt(FAb(e,2)^2 + FAb(e,3)^2)^2 + FAb(e,4)^2);
end

%Finding the indices for the directions
I225 = find(FAv(:,1) == -135);
I270 = find(FAv(:,1) == -90);
I315= find(FAv(:,1) == -45);
I0 = find(FAv(:,1) == 0);
I45 = find(FAv(:,1) == 45);
I90 = find(FAv(:,1) == 90);
I135 = find(FAv(:,1) == 135);
I180 = find(FAv(:,1) == 180);

% Iall = [I225 I270 I315 I0 I45 I90 I135 I180];
Iall = [I0 I45 I90 I135 I180 I225 I270 I315];
%group average vectors (visual trials)
G225xv(:,p) = FAv(I225,2);
G270xv(:,p) = FAv(I270,2);
G315xv(:,p) = FAv(I315,2);
G0xv(:,p) = FAv(I0,2);
G45xv(:,p) = FAv(I45,2);
G90xv(:,p) = FAv(I90,2);
G135xv(:,p) = FAv(I135,2);
G180xv(:,p) = FAv(I180,2);

G225yv(:,p) = FAv(I225,3);
G270yv(:,p) = FAv(I270,3);
G315yv(:,p) = FAv(I315,3);
G0yv(:,p) = FAv(I0,3);
G45yv(:,p) = FAv(I45,3);
G90yv(:,p) = FAv(I90,3);
G135yv(:,p) = FAv(I135,3);
G180yv(:,p) = FAv(I180,3);

%group average vectors (reproduction trials)
G225x(:,p) = FAb(I225,2);
G270x(:,p) = FAb(I270,2);
G315x(:,p) = FAb(I315,2);
G0x(:,p) = FAb(I0,2);
G45x(:,p) = FAb(I45,2);
G90x(:,p) = FAb(I90,2);
G135x(:,p) = FAb(I135,2);
G180x(:,p) = FAb(I180,2);

Fbmx = [mean(G225x(:,p)); mean(G270x(:,p)); mean(G315x(:,p)); mean(G0x(:,p)); mean(G45x(:,p)); mean(G90x(:,p)); mean(G135x(:,p)); mean(G180x(:,p))];

G225y(:,p) = FAb(I225,3);
G270y(:,p) = FAb(I270,3);
G315y(:,p) = FAb(I315,3);
G0y(:,p) = FAb(I0,3);
G45y(:,p) = FAb(I45,3);
G90y(:,p) = FAb(I90,3);
G135y(:,p) = FAb(I135,3);
G180y(:,p) = FAb(I180,3);

Fbmy = [mean(G225y(:,p)); mean(G270y(:,p)); mean(G315y(:,p)); mean(G0y(:,p)); mean(G45y(:,p)); mean(G90y(:,p)); mean(G135y(:,p)); mean(G180y(:,p))];


FREr(p,:,q) = [mean(sqrt(FAb(I0,2).^2 + FAb(I0,3).^2)) - mean(sqrt(FAv(I0,2).^2 + FAv(I0,3).^2)), ...
    mean(sqrt(FAb(I45,2).^2 + FAb(I45,3).^2)) - mean(sqrt(FAv(I45,2).^2 + FAv(I45,3).^2)), ...
    mean(sqrt(FAb(I90,2).^2 + FAb(I90,3).^2)) - mean(sqrt(FAv(I90,2).^2 + FAv(I90,3).^2)), ...
    mean(sqrt(FAb(I135,2).^2 + FAb(I135,3).^2)) - mean(sqrt(FAv(I135,2).^2 + FAv(I135,3).^2)), ...
    mean(sqrt(FAb(I180,2).^2 + FAb(I180,3).^2)) - mean(sqrt(FAv(I180,2).^2 + FAv(I180,3).^2)), ...
    mean(sqrt(FAb(I225,2).^2 + FAb(I225,3).^2)) - mean(sqrt(FAv(I225,2).^2 + FAv(I225,3).^2)), ...
    mean(sqrt(FAb(I270,2).^2 + FAb(I270,3).^2)) - mean(sqrt(FAv(I270,2).^2 + FAv(I270,3).^2)), ...
    mean(sqrt(FAb(I315,2).^2 + FAb(I315,3).^2)) - mean(sqrt(FAv(I315,2).^2 + FAv(I315,3).^2))];

FRE_an(p,:,q) = [mean(atan2(FAv(I0,3),FAv(I0,2))*180/pi - atan2(FAb(I0,3),FAb(I0,2))*180/pi),...
    mean(atan2(FAv(I45,3),FAv(I45,2))*180/pi - atan2(FAb(I45,3),FAb(I45,2))*180/pi),...
    mean(atan2(FAv(I90,3),FAv(I90,2))*180/pi - atan2(FAb(I90,3),FAb(I90,2))*180/pi),...
    mean(atan2(FAv(I135,3),FAv(I135,2))*180/pi - atan2(FAb(I135,3),FAb(I135,2))*180/pi),...
    mean(abs(atan2(FAv(I180,3),FAv(I180,2))*180/pi) - abs(atan2(FAb(I180,3),FAb(I180,2))*180/pi)),...
    mean(atan2(FAv(I225,3),FAv(I225,2))*180/pi - atan2(FAb(I225,3),FAb(I225,2))*180/pi),...
    mean(atan2(FAv(I270,3),FAv(I270,2))*180/pi - atan2(FAb(I270,3),FAb(I270,2))*180/pi),...
    mean(atan2(FAv(I315,3),FAv(I315,2))*180/pi - atan2(FAb(I315,3),FAb(I315,2))*180/pi)];

%% Position ellipses individual
%calculating position ellipses per person and plot them
[EL_p, xm_p, ym_p, phi, Ew_p, Kel] = ellipse_fit(FAb(:,2),FAb(:,3));
[EL_pv, xm_pv, ym_pv, phiv, Ew_pv, Kelv] = ellipse_fit(FAv(:,2),FAv(:,3));
Opp(p,q) = pi*Ew_p(1)*Ew_p(2);
Ratio(p,q) = Ew_p(1)/Ew_p(2);
Orient(p,q) = phi*180/pi;
EL_p(1,:) = EL_p(1,:) + xm_p;
EL_p(2,:) = EL_p(2,:) + ym_p;
EL_pv(1,:) = EL_pv(1,:) + xm_pv;
EL_pv(2,:) = EL_pv(2,:) + ym_pv;

%2 vectors, (:,1) angle of force estimates and (:,2) magnitude of force
%estimates
% ANG_MAG = [atan2(FAb(:,3),FAb(:,2))*180/pi sqrt(FAb(:,2).^2 + FAb(:,3).^2)];
% C = find(ANG_MAG < 0);
% ANG_MAG(C,1) = ANG_MAG(C,1) + 360;
% 
% for i = 1:length(ANG_MAG)
%     if round(ANG_MAG(i,1) - phi) > 0
%         an(i,1) = round(ANG_MAG(i,1) - phi);
%     elseif round(ANG_MAG(i,1) - phi) <= 0
%         an(i,1) = round(ANG_MAG(i,1) + 360 - phi);
%     end
% end
% 
% %same for ellipse (at the same angles, rounded)
% EL_fit = [ANG_MAG(:,1) sqrt(EL_p(1,an)'.^2 + EL_p(2,an)'.^2)];
% 
% AM = sortrows(ANG_MAG,1);
% EF = sortrows(EL_fit,1);
% 
% rs_fit = fitlm(AM(:,1),AM(:,2));
% 
% figure(99 + q)
% subplot(4,3,p)
% plot(AM(:,1),AM(:,2),'b*'), hold on
% plot(EF(:,1),EF(:,2),'r')
% xlabel({num2str(rs_fit.Rsquared.Ordinary)})

figure(q)
subplot(4,3,p)
plot(EL_p(1,:),EL_p(2,:),'color', 'r','linewidth',2), hold on
plot(FAb(:,2),FAb(:,3),'r*')
plot(Fbmx,Fbmy,'g*')
plot(EL_pv(1,:) + xm_pv,EL_pv(2,:) + ym_pv ,'color', 'b','linewidth',2), hold on
plot(xm_p,ym_p,'marker','+','color','r','linewidth',2)
% xlabel('Force [N]','fontsize',14)
% ylabel('Force [N]','fontsize',14)
% txt1 = ['Position_' num2str(q)];
% title(txt1,'fontsize',14)
axis equal

% figure(50+q)
% subplot(4,3,p)
% plot(EL_p(1,:),EL_p(2,:),'color', 'r','linewidth',2), hold on
% plot(FAb(:,2),FAb(:,3),'r*')
% plot(Fbmx,Fbmy,'g*')
% plot(EL_pv(1,:),EL_pv(2,:) ,'color', 'b','linewidth',2), hold on
% plot(xm_p,ym_p,'marker','+','color','r','linewidth',2)
% xlabel('Force [N]','fontsize',14)
% ylabel('Force [N]','fontsize',14)
% txt1 = ['Position_' num2str(q)];
% title(txt1,'fontsize',14)
% axis equal
%% joint contributions
%data of the subjects

txt = ['D:\matlab\Data\Exp2\data\' char(datavec1(p)) '_info'];
load(txt)

l1 = info(2)/1000;
l2 = info(3)/1000;
l3 = 50/1000;
tsh1 = info(5)*pi/180;
tsh2 = info(6)*pi/180;
tel1 = info(7)*pi/180;
tel2 = info(8)*pi/180;
twr = 0;

Jinfo(:,p) = [l1; l2; tsh1; tsh2; tel1; tel2];
 %calculating the joint contributions by taking one point of the ellipse and
%transforming by using the jacobian

%The ellipses are rotated and start in the direction phi, here the the
%vectors are shifted
for d = 1:360
    
    kv = d + round(phiv*180/pi);
    k = d + round(phi*180/pi);
    if k > 360
        e = k - 360;
    else
        e = k;
    end
    if kv >= 361
        ev = kv - 360;
    else
        ev = kv;
    end
    
    FRE(p,e,q) = sqrt(EL_p(1,d).^2 + EL_p(2,d).^2) - sqrt(EL_pv(1,d).^2 + EL_pv(2,d).^2);
    EL_s(:,e) = EL_p(:,d);
    %calculating the joint contributions (the FRJp vectors are used for
    %group averages)
    if q == 1

%         J = [-sin(tsh1)*l1-sin(tsh1+tel1)*l2,     -sin(tsh1+tel1)*l2;
%             cos(tsh1)*l1+cos(tsh1+tel1)*l2,      cos(tsh1+tel1)*l2];
        
        J = [-sin(tsh1)*l1-sin(tsh1+tel1)*l2-sin(tsh1+tel1+twr)*l3,     -sin(tsh1+tel1)*l2-sin(tsh1+tel1+twr)*l3 -sin(tsh1+tel1+twr)*l3;
            cos(tsh1)*l1+cos(tsh1+tel1)*l2+cos(tsh1+tel1+twr)*l3,      cos(tsh1+tel1)*l2+cos(tsh1+tel1+twr)*l3 cos(tsh1+tel1+twr)*l3];
   
        FRJp1(:,e,p) = (J' * EL_p(:,d));
        VFRJp1(:,ev,p) = (J' * EL_pv(:,d));
        
%         dSS_p1(:,p) = regress(FRJp1(1,:,:)',[VFRJp1(1,:,:)' ones(360,1)])
%         dSS_p1(:,p) = regress(FRJp1(1,:,:)',[VFRJp1(1,:,:)' ones(360,1)])
    elseif q == 2
        
%         J = [-sin(tsh1)*l1-sin(tsh1+tel2)*l2,     -sin(tsh1+tel2)*l2;
%             cos(tsh1)*l1+cos(tsh1+tel2)*l2,      cos(tsh1+tel2)*l2];
        
        J = [-sin(tsh1)*l1-sin(tsh1+tel2)*l2-sin(tsh1+tel2+twr)*l3,     -sin(tsh1+tel2)*l2-sin(tsh1+tel2+twr)*l3 -sin(tsh1+tel2+twr)*l3;
            cos(tsh1)*l1+cos(tsh1+tel2)*l2+cos(tsh1+tel2+twr)*l3,      cos(tsh1+tel2)*l2+cos(tsh1+tel2+twr)*l3 cos(tsh1+tel2+twr)*l3];
        
        FRJp2(:,e,p) = (J' * EL_p(:,d));
        VFRJp2(:,ev,p) = (J' * EL_pv(:,d));
        
    elseif q == 3
       
%         J = [-sin(tsh2)*l1-sin(tsh2+tel1)*l2,     -sin(tsh2+tel1)*l2;
%             cos(tsh2)*l1+cos(tsh2+tel1)*l2,      cos(tsh2+tel1)*l2];
        
        J = [-sin(tsh2)*l1-sin(tsh2+tel1)*l2-sin(tsh2+tel1+twr)*l3,     -sin(tsh2+tel1)*l2-sin(tsh2+tel1+twr)*l3 -sin(tsh2+tel1+twr)*l3;
            cos(tsh2)*l1+cos(tsh2+tel1)*l2+cos(tsh2+tel1+twr)*l3,      cos(tsh2+tel1)*l2+cos(tsh2+tel1+twr)*l3 cos(tsh2+tel1+twr)*l3];
        
        FRJp3(:,e,p) = (J' * EL_p(:,d));
        VFRJp3(:,ev,p) = (J' * EL_pv(:,d));
        
    elseif q == 4
     
%         J = [-sin(tsh2)*l1-sin(tsh2+tel2)*l2,     -sin(tsh2+tel2)*l2;
%             cos(tsh2)*l1+cos(tsh2+tel2)*l2,      cos(tsh2+tel2)*l2];
        
        J = [-sin(tsh2)*l1-sin(tsh2+tel2)*l2-sin(tsh2+tel2+twr)*l3,     -sin(tsh2+tel2)*l2-sin(tsh2+tel2+twr)*l3 -sin(tsh2+tel2+twr)*l3;
            cos(tsh2)*l1+cos(tsh2+tel2)*l2+cos(tsh2+tel2+twr)*l3,      cos(tsh2+tel2)*l2+cos(tsh2+tel2+twr)*l3 cos(tsh2+tel2+twr)*l3];

        FRJp4(:,e,p) = (J' * EL_p(:,d));
        VFRJp4(:,ev,p) = (J' * EL_pv(:,d));
        
    end
end
%%
% % global Jinfo
% options=optimset('tolfun',1e-10,'diffminchange',1e-5,'maxfunevals',20000,'maxiter',1000,'display','iter','tolx',1e-5);
% a = 1.1;
% b = 0.1;
% c = 1.2;
% % Ps = [1.3 0.2 1.4];
% x0 = [1.1 0.1 1.4 0.01 0.01];
% lb = [1 -2 1 0 0]%./x0;
% ub = [2 2 2 1 1]%./x0;
% % x0 = [1.1 1.4 0.01 0.01];
% % lb = [1 1 0 0 ]./x0;
% % ub = [2 2 1 1]./x0;
% 
% % e = erfit(a,b,c,Jinfo,ELf)
% P = lsqnonlin('erfits',x0,[],[],options,EL_s,J);
% Ptot(:,p,q) = P;
% 
% %predict per subject
% n = 0:0.0175:2*pi;
% tar = 10 * [cos(n); sin(n)];
% 
% Joint = J' * tar;  
% 
% SS = P(1);
% EE = P(3);
% Bi = P(2);
% SO = P(4);
% EO = P(5);
% 
% %with bi-articular muscles
% Joint_reps = [SS   Bi; Bi EE ] * Joint;
% Joint_reps(1,:) = Joint_reps(1,:) + SO;
% Joint_reps(2,:) = Joint_reps(2,:) + EO;
% rep_el = inv(J') * Joint_reps;
% 
% %directions of mono- and bi-articular muscles
% sdir = inv(J') * [1 0; 0 0];
% edir = inv(J') * [0 0; 0 1];
% bdir = inv(J') * [1 1; 1 1];
% MS = sdir * [cos(n); sin(n)];
% ME = edir * [cos(n); sin(n)];
% BB = bdir * [cos(n); sin(n)];
% 
% figure(q)
% subplot(4,3,p)
% plot(rep_el(1,:),rep_el(2,:),'color', 'c','linewidth',2), hold on
% plot(MS(1,:),MS(2,:),'color', 'm'), hold on
% plot(ME(1,:),ME(2,:),'color', 'g','linewidth',2), hold on
% plot(BB(1,:),BB(2,:),'color', 'k','linewidth',2), hold on
% axis equal
%%
if q == 1
    dSS_p1(:,p) = regress(FRJp1(1,:,p)',[VFRJp1(1,:,p)' ones(360,1)]);
    dEE_p1(:,p) = regress(FRJp1(2,:,p)',[VFRJp1(2,:,p)' ones(360,1)]);
elseif q == 2
    dSS_p2(:,p) = regress(FRJp2(1,:,p)',[VFRJp2(1,:,p)' ones(360,1)]);
    dEE_p2(:,p) = regress(FRJp2(2,:,p)',[VFRJp2(2,:,p)' ones(360,1)]);
elseif q == 3
    dSS_p3(:,p) = regress(FRJp3(1,:,p)',[VFRJp3(1,:,p)' ones(360,1)]);
    dEE_p3(:,p) = regress(FRJp3(2,:,p)',[VFRJp3(2,:,p)' ones(360,1)]);
elseif q == 4
    dSS_p4(:,p) = regress(FRJp4(1,:,p)',[VFRJp4(1,:,p)' ones(360,1)]);
    dEE_p4(:,p) = regress(FRJp4(2,:,p)',[VFRJp4(2,:,p)' ones(360,1)]);
end



%plotting the joint contributions
% figure(6)
%     if q == 1
%     subplot(221)
%     polar(plxp, sc*ones(1,length(plxp)),'w'), hold on
%     polar(plxp, abs([VFRJp1(1,:,p) VFRJp1(1,1,p)]),'b'), hold all
%     polar(plxp, abs([VFRJp1(2,:,p) VFRJp1(2,1,p)]),'r')
%     polar(plxp, abs([FRJp1(1,:,p) FRJp1(1,1,p)]),'c')
%     polar(plxp, abs([FRJp1(2,:,p) FRJp1(2,1,p)]),'m')
%     title('Position 1')
%     
%     elseif q == 2
%     subplot(222)
%     polar(plxp, sc*ones(1,length(plxp)),'w'), hold on
%     polar(plxp, abs([VFRJp2(1,:,p) VFRJp2(1,1,p)]),'b'), hold all
%     polar(plxp, abs([VFRJp2(2,:,p) VFRJp2(2,1,p)]),'r')
%     polar(plxp, abs([FRJp2(1,:,p) FRJp2(1,1,p)]),'c')
%     polar(plxp, abs([FRJp2(2,:,p) FRJp2(2,1,p)]),'m')
%     title('Position 2')
%     
%     elseif q == 3
%     subplot(223)
%     polar(plxp, sc*ones(1,length(plxp)),'w'), hold on
%     polar(plxp, abs([VFRJp3(1,:,p) VFRJp3(1,1,p)]),'b'), hold all
%     polar(plxp, abs([VFRJp3(2,:,p) VFRJp3(2,1,p)]),'r')
%     polar(plxp, abs([FRJp3(1,:,p) FRJp3(1,1,p)]),'c')
%     polar(plxp, abs([FRJp3(2,:,p) FRJp3(2,1,p)]),'m')
%     title('Position 3')
%     
%     elseif q == 4
%     subplot(224)
%     polar(plxp, sc*ones(1,length(plxp)),'w'), hold on
%     polar(plxp, abs([VFRJp4(1,:,p) VFRJp4(1,1,p)]),'b'), hold all
%     polar(plxp, abs([VFRJp4(2,:,p) VFRJp4(2,1,p)]),'r')
%     polar(plxp, abs([FRJp4(1,:,p) FRJp4(1,1,p)]),'c')
%     polar(plxp, abs([FRJp4(2,:,p) FRJp4(2,1,p)]),'m')
%     title('Position 4')
%     end
%     

end




% Group averages
markv = {'r*' 'ro' 'r*' 'ro' 'r*' 'ro' 'r*' 'ro' 'r*' 'ro' 'r*' 'ro'};
%plot average angle and force reproduction
Xvtot = [mean(G225xv,1) mean(G270xv,1) mean(G315xv,1) mean(G0xv,1) mean(G45xv,1) mean(G90xv,1) mean(G135xv,1) mean(G180xv,1)]';
Yvtot = [mean(G225yv,1) mean(G270yv,1) mean(G315yv,1) mean(G0yv,1) mean(G45yv,1) mean(G90yv,1) mean(G135yv,1) mean(G180yv,1)]';
% Xbtot = [mean(mean(G225x,1)) mean(mean(G270x,1)) mean(mean(G315x,1)) mean(mean(G0x,1)) mean(mean(G45x,1)) mean(mean(G90x,1)) mean(mean(G135x,1)) mean(mean(G180x,1))]';
% Ybtot = [mean(mean(G225y,1)) mean(mean(G270y,1)) mean(mean(G315y,1)) mean(mean(G0y,1)) mean(mean(G45y,1)) mean(mean(G90y,1)) mean(mean(G135y,1)) mean(mean(G180y,1))]';
Xbtot = [mean(G225x,1) mean(G270x,1) mean(G315x,1) mean(G0x,1) mean(G45x,1) mean(G90x,1) mean(G135x,1) mean(G180x,1)]';
Ybtot = [mean(G225y,1) mean(G270y,1) mean(G315y,1) mean(G0y,1) mean(G45y,1) mean(G90y,1) mean(G135y,1) mean(G180y,1)]';

xd_g = [mean(G225x,1); mean(G270x,1); mean(G315x,1); mean(G0x,1); mean(G45x,1); mean(G90x,1); mean(G135x,1); mean(G180x,1)]';
yd_g = [mean(G225y,1); mean(G270y,1); mean(G315y,1); mean(G0y,1); mean(G45y,1); mean(G90y,1); mean(G135y,1); mean(G180y,1)]';


RFv = [sqrt(mean(G225xv,1).^2 + mean(G225yv,1).^2); sqrt(mean(G270xv,1).^2 + mean(G270yv,1).^2); sqrt(mean(G315xv,1).^2 + mean(G315yv,1).^2); sqrt(mean(G0xv,1).^2 + mean(G0yv,1).^2); sqrt(mean(G45xv,1).^2 + mean(G45yv,1).^2); sqrt(mean(G90xv,1).^2 + mean(G90yv,1).^2); sqrt(mean(G135xv,1).^2 + mean(G135yv,1).^2); sqrt(mean(G180xv,1).^2 + mean(G180yv,1).^2)];
RF = [sqrt(mean(G225x,1).^2 + mean(G225y,1).^2); sqrt(mean(G270x,1).^2 + mean(G270y,1).^2); sqrt(mean(G315x,1).^2 + mean(G315y,1).^2); sqrt(mean(G0x,1).^2 + mean(G0y,1).^2); sqrt(mean(G45x,1).^2 + mean(G45y,1).^2); sqrt(mean(G90x,1).^2 + mean(G90y,1).^2); sqrt(mean(G135x,1).^2 + mean(G135y,1).^2); sqrt(mean(G180x,1).^2 + mean(G180y,1).^2)];

% target = [0 45 90 135 180 -135 -90 -45];


figure(20)
subplot(2,2,q)
for i = 1:8
[EL_dg] = confell95(xd_g(:,i),yd_g(:,i));
% plot([0 Xvtot(i)],[0 Yvtot(i)],'color','b','linewidt',2), hold on
% plot([0 Xbtot(i)],[0 Ybtot(i)],'color','r','linewidt',2), hold on
plot(xd_g(:,i),yd_g(:,i),char(markv(i)),'color',[128 128 128]/255), hold on
% set(gca,'Facealpha',0.5)
% plot(EL_dg(1,:),EL_dg(2,:),'r')
end

% [EL_pgr, xm_pg, ym_pg, phig, Ew_pg, Kelg] = ellipse_fit(Xbtot,Ybtot);
% [EL_pgv, xm_pgv, ym_pgv, phigv, Ew_pgv, Kelgv] = ellipse_fit(Xvtot,Yvtot);
[EL_pg, xm_pg, ym_pg, phig, Ew_pg, Kelg] = ellipse_fit(Xbtot,Ybtot);
[EL_pgv, xm_pgv, ym_pgv, phigv, Ew_pgv, Kelgv] = ellipse_fit(Xvtot,Yvtot);

dx(q) = xm_pg;
dy(q) = ym_pg;
EL_pgv(1,:) = EL_pgv(1,:) + xm_pgv;
EL_pgv(2,:) = EL_pgv(2,:) + ym_pgv;

% for ii = 1:360
%     k = ii + round(phig*180/pi);
%     if k > 360
%         e = k - 360;
%     else
%         e = k;
%     end
%     
%     EL_pg(:,e) = EL_pgr(:,ii);
%     
% end

Xf(:,:,q) = [Xbtot Ybtot Xvtot Yvtot];
EL_pg(1,:) = EL_pg(1,:) + xm_pg;
EL_pg(2,:) = EL_pg(2,:) + ym_pg;
ELf(:,:,q) = EL_pg;
phiELf(q) = phig;

or_x = -30:30;
or_line = tan(phig)*or_x;

plot(EL_pg(1,:),EL_pg(2,:),'r','linewidt',3)
plot(EL_pgv(1,:),EL_pgv(2,:),'b','linewidt',3)
plot(xm_pg,ym_pg,'r+','markersize',10);
plot(xm_pgv,ym_pgv,'k+','markersize',40);
plot(or_x + xm_pg,or_line + ym_pg,'k-.')
% plot([0 Xbtot(4)],[0 Ybtot(4)],'color','k','linewidt',2), hold on
% plot([Xvtot; Xvtot(1)],[Yvtot; Yvtot(1)],'color','b','linewidt',2), hold on
% plot([Xbtot; Xbtot(1)],[Ybtot; Ybtot(1)],'color','r','linewidt',2), hold on
axis equal
axis([-fl fl -fl fl])
axis equal
box off
title(['Position ' num2str(q)])

% EL_pg(1,:) = EL_pg(1,:) + xm_pg;
% EL_pg(2,:) = EL_pg(2,:) + ym_pg;
% EL_pgv(1,:) = EL_pgv(1,:) + xm_pgv;
% EL_pgv(2,:) = EL_pgv(2,:) + ym_pgv;


% FRE_pol = sqrt(((EL_pg(1,:) + xm_pg) - (EL_pgv(1,:) + xm_pgv)).^2 + ((EL_pg(2,:) + ym_pg) - (EL_pgv(2,:) + ym_pgv)).^2);
% FRE_pol
% figure(31)
% for i = 1:4
% subplot(2,2,i)    
% polar(1:360,FRE_pol ,'g'); hold on
% end

figure(7)
    if q == 1
    subplot(221)
    polar(plxp, sc*ones(1,length(plxp)),'w'), hold on
    polar(plxp, abs([mean(VFRJp1(1,:,:),3) mean(VFRJp1(1,1,:),3)]),'b--'), hold all
%     polar(plxp, abs([mean(VFRJp1(1,:,:),3)+std(VFRJp1(1,:,:),0,3) mean(VFRJp1(1,1,:),3)+std(VFRJp1(1,1,:),0,3)]),'b--')
%     polar(plxp, abs([mean(VFRJp1(1,:,:),3)-std(VFRJp1(1,:,:),0,3) mean(VFRJp1(1,1,:),3)-std(VFRJp1(1,1,:),0,3)]),'b--')
    polar(plxp, abs([mean(VFRJp1(2,:,:),3) mean(VFRJp1(2,1,:),3)]),'b:')
    polar(plxp, abs([mean(VFRJp1(3,:,:),3) mean(VFRJp1(3,1,:),3)]),'b')
%     polar(plxp, abs([mean(VFRJp1(2,:,:),3)+std(VFRJp1(2,:,:),0,3) mean(VFRJp1(2,1,:),3)+std(VFRJp1(2,1,:),0,3)]),'r--')
%     polar(plxp, abs([mean(VFRJp1(2,:,:),3)-std(VFRJp1(2,:,:),0,3) mean(VFRJp1(2,1,:),3)-std(VFRJp1(2,1,:),0,3)]),'r--')
    polar(plxp, abs([mean(FRJp1(1,:,:),3) mean(FRJp1(1,1,:),3)]),'r--')
    %     polar(plxp, abs([mean(FRJp1(1,:,:),3)+std(FRJp1(1,:,:),0,3) mean(FRJp1(1,1,:),3)+std(FRJp1(1,1,:),0,3)]),'c--')
%     polar(plxp, abs([mean(FRJp1(1,:,:),3)-std(FRJp1(1,:,:),0,3) mean(FRJp1(1,1,:),3)-std(FRJp1(1,1,:),0,3)]),'c--')
    polar(plxp, abs([mean(FRJp1(2,:,:),3) mean(FRJp1(2,1,:),3)]),'r:')
    polar(plxp, abs([mean(FRJp1(3,:,:),3) mean(FRJp1(3,1,:),3)]),'r')
%     polar(plxp, abs([mean(FRJp1(2,:,:),3)+std(FRJp1(2,:,:),0,3) mean(FRJp1(2,1,:),3)+std(FRJp1(2,1,:),0,3)]),'m--')
%     polar(plxp, abs([mean(FRJp1(2,:,:),3)-std(FRJp1(2,:,:),0,3) mean(FRJp1(2,1,:),3)-std(FRJp1(2,1,:),0,3)]),'m--')
    title('Position 1')
    
    elseif q == 2
    subplot(222)
    polar(plxp, sc*ones(1,length(plxp)),'w'), hold on
    polar(plxp, abs([mean(VFRJp2(1,:,:),3) mean(VFRJp2(1,1,:),3)]),'b--'), hold all
    polar(plxp, abs([mean(VFRJp2(2,:,:),3) mean(VFRJp2(2,1,:),3)]),'b:')
    polar(plxp, abs([mean(VFRJp2(3,:,:),3) mean(VFRJp2(3,1,:),3)]),'b')
    polar(plxp, abs([mean(FRJp2(1,:,:),3) mean(FRJp2(1,1,:),3)]),'r--')
    polar(plxp, abs([mean(FRJp2(2,:,:),3) mean(FRJp2(2,1,:),3)]),'r:')
    polar(plxp, abs([mean(FRJp2(3,:,:),3) mean(FRJp2(3,1,:),3)]),'r')
    title('Position 2')
    legend('','SH vis', 'EL vis', 'SH rep', 'EL rep')
    
    elseif q == 3
    subplot(223)
     polar(plxp, sc*ones(1,length(plxp)),'w'), hold on
    polar(plxp, abs([mean(VFRJp3(1,:,:),3) mean(VFRJp3(1,1,:),3)]),'b--'), hold all
    polar(plxp, abs([mean(VFRJp3(2,:,:),3) mean(VFRJp3(2,1,:),3)]),'b:')
    polar(plxp, abs([mean(VFRJp3(3,:,:),3) mean(VFRJp3(3,1,:),3)]),'b')
    polar(plxp, abs([mean(FRJp3(1,:,:),3) mean(FRJp3(1,1,:),3)]),'r--')
    polar(plxp, abs([mean(FRJp3(2,:,:),3) mean(FRJp3(2,1,:),3)]),'r:')
    polar(plxp, abs([mean(FRJp3(3,:,:),3) mean(FRJp3(3,1,:),3)]),'r')
    title('Position 3')
    
    elseif q == 4
    subplot(224)
     polar(plxp, sc*ones(1,length(plxp)),'w'), hold on
    polar(plxp, abs([mean(VFRJp4(1,:,:),3) mean(VFRJp4(1,1,:),3)]),'b--'), hold all
    polar(plxp, abs([mean(VFRJp4(2,:,:),3) mean(VFRJp4(2,1,:),3)]),'b:')
    polar(plxp, abs([mean(VFRJp4(3,:,:),3) mean(VFRJp4(3,1,:),3)]),'b')
    polar(plxp, abs([mean(FRJp4(1,:,:),3) mean(FRJp4(1,1,:),3)]),'r--')
    polar(plxp, abs([mean(FRJp4(2,:,:),3) mean(FRJp4(2,1,:),3)]),'r:')
    polar(plxp, abs([mean(FRJp4(3,:,:),3) mean(FRJp4(3,1,:),3)]),'r')
    title('Position 4')
    end
    
    figure(100)
    if q == 1
    subplot(221)
     polar(plxp, ones(1,length(plxp)),'w'), hold on
    polar(plxp, abs([mean(VFRJp1(3,:,:),3) mean(VFRJp1(3,1,:),3)]),'b'), hold on
    polar(plxp, abs([mean(FRJp1(3,:,:),3) mean(FRJp1(3,1,:),3)]),'r')

    title('Position 1')
    
    elseif q == 2
    subplot(222)
      polar(plxp, ones(1,length(plxp)),'w'), hold on
    polar(plxp, abs([mean(VFRJp2(3,:,:),3) mean(VFRJp2(3,1,:),3)]),'b'), hold on
    polar(plxp, abs([mean(FRJp2(3,:,:),3) mean(FRJp2(3,1,:),3)]),'r')
    title('Position 2')
    legend('','SH vis', 'EL vis', 'SH rep', 'EL rep')
    
    elseif q == 3
    subplot(223)
      polar(plxp, ones(1,length(plxp)),'w'), hold on
    polar(plxp, abs([mean(VFRJp3(3,:,:),3) mean(VFRJp3(3,1,:),3)]),'b'), hold on
    polar(plxp, abs([mean(FRJp3(3,:,:),3) mean(FRJp3(3,1,:),3)]),'r')
    title('Position 3')
    
    elseif q == 4
    subplot(224)
      polar(plxp, ones(1,length(plxp)),'w'), hold on
    polar(plxp, abs([mean(VFRJp4(3,:,:),3) mean(VFRJp4(3,1,:),3)]),'b'), hold on
    polar(plxp, abs([mean(FRJp4(3,:,:),3) mean(FRJp4(3,1,:),3)]),'r')
    title('Position 4')
    end
end

for p = 1:12
        
    MA(p,:,1) = abs([FRJp1(1,round(Orient(p,1)),p) FRJp2(1,round(Orient(p,2)),p) FRJp3(1,round(Orient(p,3)),p) FRJp4(1,round(Orient(p,4)),p)]);
    MA(p,:,2) = abs([FRJp1(2,round(Orient(p,1)),p) FRJp2(2,round(Orient(p,2)),p) FRJp3(2,round(Orient(p,3)),p) FRJp4(2,round(Orient(p,4)),p)]);
    MI(p,:,1) = abs([FRJp1(1,round(Orient(p,1)+90),p) FRJp2(1,round(Orient(p,2)+90),p) FRJp3(1,round(Orient(p,3)+90),p) FRJp4(1,round(Orient(p,4)+90),p)]);
    MI(p,:,2) = abs([FRJp1(2,round(Orient(p,1)+90),p) FRJp2(2,round(Orient(p,2)+90),p) FRJp3(2,round(Orient(p,3)+90),p) FRJp4(2,round(Orient(p,4)+90),p)]);
end
    
NC = [MA(:,:,1)./MI(:,:,1) MA(:,:,2)./MI(:,:,2)];
NCm = [mean(MA(:,:,1)./MI(:,:,1),2) mean(MA(:,:,2)./MI(:,:,2),2)];
xlswrite('NCm',NCm)

target = [0 45 90 135 180 225 270 315];
% figure(100)
% for i = 1:4
%     subplot(2,2,i)
% errorbar(target,mean(SE_son(:,:,i),1),std(SE_son(:,:,i),0,1)/12,'r'), hold on
% errorbar(target,mean(SE_eon(:,:,i),1),std(SE_eon(:,:,i),0,1)/12,'b')
% errorbar(target,mean(SE_soff(:,:,i),1),std(SE_soff(:,:,i),0,1)/12,'r--'), hold on
% errorbar(target,mean(SE_eoff(:,:,i),1),std(SE_eoff(:,:,i),0,1)/12,'b--')
% end
FREc_p1 = FRE(:,1:45:end,1);
FREc_p2 = FRE(:,1:45:end,2);
FREc_p3 = FRE(:,1:45:end,3);
FREc_p4 = FRE(:,1:45:end,4);

% FREc_p1 = FRE(:,1:45:end,1);
% FREc_p2 = FRE(:,1:45:end,2);
% FREc_p3 = FRE(:,1:45:end,3);
% FREc_p4 = FRE(:,1:45:end,4);

posv = [1 2 5 6];
colv = ['r' 'b' 'g' 'k'];
figure(30)
for i = 1:4
subplot(2,2,i)    

fill([1:360, fliplr(1:360)],[(mean(FRE(:,:,i),1) + std(FRE(:,:,i),0,1)/sqrt(12))'; fliplr(mean(FRE(:,:,i),1) - std(FRE(:,:,i),0,1)/sqrt(12))'],'g','linewidt',2), hold on
% errorbar(target,mean(FRE(:,1:45:end,i),1),std(FRE(:,1:45:end,i),0,1)/sqrt(12),'m','linewidt',2), hold on
% errorbar(target,mean(FREr(:,:,i),1),std(FREr(:,:,i),0,1)/sqrt(12),'b','linewidt',2)
errorbarxy(target + mean(FRE_an(:,:,i)), mean(FREr(:,:,i)),std(FRE_an(:,:,i),0,1)/sqrt(12), std(FRE_an(:,:,i),0,1)/sqrt(12), std(FREr(:,:,i),0,1)/sqrt(12), std(FREr(:,:,i),0,1)/sqrt(12),{'k.' 'k' 'k'})
ylabel('FRE [N]')
xlabel('force directions [deg.]')
axis([-5 360 0 8])
box off
end

figure(31)
for i = 1:4
subplot(2,2,i)    
polar(plxp,abs([mean(FRE(:,:,i),1) mean(FRE(:,1,i),1)]),'g'); hold on
end

target = [0 45 90 135 180 225 270 315];
%     figure(10)
%     subplot(3,2,[1 2])
%     errorbar(target,mean(mean(SJE_p1(:,:,:),1),3),std(mean(SJE_p1(:,:,:),1),0,3)/sqrt(12),'r','linewidt',2), hold on
%     errorbar(target,mean(mean(EJE_p1(:,:,:),1),3),std(mean(EJE_p1(:,:,:),1),0,3)/sqrt(12),'b','linewidt',2)
%     errorbar(target,mean(mean(SJE_p2(:,:,:),1),3),std(mean(SJE_p2(:,:,:),1),0,3)/sqrt(12),'g','linewidt',2), hold on
%     errorbar(target,mean(mean(EJE_p2(:,:,:),1),3),std(mean(EJE_p2(:,:,:),1),0,3)/sqrt(12),'c','linewidt',2)
%     errorbar(target,mean(mean(SJE_p3(:,:,:),1),3),std(mean(SJE_p3(:,:,:),1),0,3)/sqrt(12),'r--','linewidt',2), hold on
%     errorbar(target,mean(mean(EJE_p3(:,:,:),1),3),std(mean(EJE_p3(:,:,:),1),0,3)/sqrt(12),'b--','linewidt',2)
%     errorbar(target,mean(mean(SJE_p4(:,:,:),1),3),std(mean(SJE_p4(:,:,:),1),0,3)/sqrt(12),'g--','linewidt',2), hold on
%     errorbar(target,mean(mean(EJE_p4(:,:,:),1),3),std(mean(EJE_p4(:,:,:),1),0,3)/sqrt(12),'c--','linewidt',2)
%     legend('Shoulder P1', 'Elbow P1', 'Shoulder P2', 'Elbow P2', 'Shoulder P3', 'Elbow P3', 'Shoulder P4', 'Elbow P4')
%     
    figure(40)
    subplot(2,2,1)
%     errorbar(target,mean(mean(SJE_p1(:,:,:),1),3),std(mean(SJE_p1(:,:,:),1),0,3)/sqrt(12),'m','linewidt',2), hold on
%     errorbar(target,mean(mean(EJE_p1(:,:,:),1),3),std(mean(EJE_p1(:,:,:),1),0,3)/sqrt(12),'c','linewidt',2)
    errorbar(target,mean((FRJp1(1,1:45:end,:) - VFRJp1(1,1:45:end,:)),3),std((FRJp1(1,1:45:end,:) - VFRJp1(1,1:45:end,:)),0,3)/sqrt(12),'m','linewidt',2), hold on
    errorbar(target,mean((FRJp1(2,1:45:end,:) - VFRJp1(2,1:45:end,:)),3),std((FRJp1(2,1:45:end,:) - VFRJp1(2,1:45:end,:)),0,3)/sqrt(12),'c','linewidt',2)
    ylabel('JCE [N]')
    xlabel('force directions [deg.]')
    
    plot(target,zeros(8,1),'k--')
    axis([-5 320 -2 2])
    subplot(2,2,2)
%     errorbar(target,mean(mean(SJE_p2(:,:,:),1),3),std(mean(SJE_p2(:,:,:),1),0,3)/sqrt(12),'m','linewidt',2), hold on
%     errorbar(target,mean(mean(EJE_p2(:,:,:),1),3),std(mean(EJE_p2(:,:,:),1),0,3)/sqrt(12),'c','linewidt',2)
    errorbar(target,mean((FRJp2(1,1:45:end,:) - VFRJp2(1,1:45:end,:)),3),std((FRJp2(1,1:45:end,:) - VFRJp2(1,1:45:end,:))/sqrt(12),0,3),'m','linewidt',2), hold on
    errorbar(target,mean((FRJp2(2,1:45:end,:) - VFRJp2(2,1:45:end,:)),3),std((FRJp2(2,1:45:end,:) - VFRJp2(2,1:45:end,:))/sqrt(12),0,3),'c','linewidt',2)
    ylabel('JCE [N]')
    xlabel('force directions [deg.]')
    plot(target,zeros(8,1),'k--')
    axis([-5 320 -2 2])
    subplot(2,2,3)
%     errorbar(target,mean(mean(SJE_p3(:,:,:),1),3),std(mean(SJE_p3(:,:,:),1),0,3)/sqrt(12),'m','linewidt',2), hold on
%     errorbar(target,mean(mean(EJE_p3(:,:,:),1),3),std(mean(EJE_p3(:,:,:),1),0,3)/sqrt(12),'c','linewidt',2)
    errorbar(target,mean((FRJp3(1,1:45:end,:) - VFRJp3(1,1:45:end,:)),3),std((FRJp3(1,1:45:end,:) - VFRJp3(1,1:45:end,:))/sqrt(12),0,3),'m','linewidt',2), hold on
    errorbar(target,mean((FRJp3(2,1:45:end,:) - VFRJp3(2,1:45:end,:)),3),std((FRJp3(2,1:45:end,:) - VFRJp3(2,1:45:end,:))/sqrt(12),0,3),'c','linewidt',2)
    ylabel('JCE [N]')
    xlabel('force directions [deg.]')
    plot(target,zeros(8,1),'k--')
    axis([-5 320 -2 2])
    subplot(2,2,4)
%     errorbar(target,mean(mean(SJE_p4(:,:,:),1),3),std(mean(SJE_p4(:,:,:),1),0,3)/sqrt(12),'m','linewidt',2), hold on
%     errorbar(target,mean(mean(EJE_p4(:,:,:),1),3),std(mean(EJE_p4(:,:,:),1),0,3)/sqrt(12),'c','linewidt',2)
    errorbar(target,mean((FRJp4(1,1:45:end,:) - VFRJp4(1,1:45:end,:)),3),std((FRJp4(1,1:45:end,:) - VFRJp4(1,1:45:end,:))/sqrt(12),0,3),'m','linewidt',2), hold on
    errorbar(target,mean((FRJp4(2,1:45:end,:) - VFRJp4(2,1:45:end,:)),3),std((FRJp4(2,1:45:end,:) - VFRJp4(2,1:45:end,:))/sqrt(12),0,3),'c','linewidt',2)
    ylabel('JCE [N]')
    xlabel('force directions [deg.]')
    plot(target,zeros(8,1),'k--')
    axis([-5 320 -2 2])
   
    x_des = [-8:8];
    figure(5)
    subplot(2,2,1)
    plot(abs(mean(VFRJp1(1,:,:),3)),abs(mean(FRJp1(1,:,:),3)),'m','linewidt',2), hold on
    plot(abs(mean(VFRJp1(2,:,:),3)),abs(mean(FRJp1(2,:,:),3)),'c','linewidt',2), hold on
    plot(abs(mean(VFRJp1(1,1,:),3)),abs(mean(FRJp1(1,1,:),3)),'m*','linewidt',2), hold on
    plot(abs(mean(VFRJp1(2,1,:),3)),abs(mean(FRJp1(2,1,:),3)),'c*','linewidt',2), hold on
    plot(x_des,x_des,'k--')
     axis([0 8 0 8])
    ylabel('Reproduced torque')
    axis equal
%     box off
    
    subplot(2,2,2)
    plot(abs(mean(VFRJp2(1,:,:),3)),abs(mean(FRJp2(1,:,:),3)),'m','linewidt',2), hold on
    plot(abs(mean(VFRJp2(2,:,:),3)),abs(mean(FRJp2(2,:,:),3)),'c','linewidt',2), hold on
    plot(abs(mean(VFRJp2(1,1,:),3)),abs(mean(FRJp2(1,1,:),3)),'m*','linewidt',2), hold on
    plot(abs(mean(VFRJp2(2,1,:),3)),abs(mean(FRJp2(2,1,:),3)),'c*','linewidt',2), hold on
    plot(x_des,x_des,'k--')
    axis([0 8 0 8])
    axis equal
%     box off
    
    subplot(2,2,3)
    plot(abs(mean(VFRJp3(1,:,:),3)),abs(mean(FRJp3(1,:,:),3)),'m','linewidt',2), hold on
    plot(abs(mean(VFRJp3(2,:,:),3)),abs(mean(FRJp3(2,:,:),3)),'c','linewidt',2), hold on
    plot(abs(mean(VFRJp3(1,1,:),3)),abs(mean(FRJp3(1,1,:),3)),'m*','linewidt',2), hold on
    plot(abs(mean(VFRJp3(2,1,:),3)),abs(mean(FRJp3(2,1,:),3)),'c*','linewidt',2), hold on
    plot(x_des,x_des,'k--')
    axis([0 8 0 8])
    xlabel('Desired torque')
    ylabel('Reproduced torque')
    axis equal
%     box off
    
    subplot(2,2,4)
    plot(abs(mean(VFRJp4(1,:,:),3)),abs(mean(FRJp4(1,:,:),3)),'m','linewidt',2), hold on
    plot(abs(mean(VFRJp4(2,:,:),3)),abs(mean(FRJp4(2,:,:),3)),'c','linewidt',2), hold on
    plot(abs(mean(VFRJp4(1,1,:),3)),abs(mean(FRJp4(1,1,:),3)),'m*','linewidt',2), hold on
    plot(abs(mean(VFRJp4(2,1,:),3)),abs(mean(FRJp4(2,1,:),3)),'c*','linewidt',2), hold on
    plot(x_des,x_des,'k--')
     axis([0 8 0 8])
    xlabel('Desired torque')
    axis equal
%     box off
    
    
     x_des = [-8:8];
     figure(41)
     plot(mean(VFRJp1(1,:,:),3),mean(FRJp1(1,:,:),3),'k','linewidt',2), hold on
     plot(mean(VFRJp1(1,1,:),3),mean(FRJp1(1,1,:),3),'k*','linewidt',2), hold on
     plot(mean(VFRJp2(1,:,:),3),mean(FRJp2(1,:,:),3),'color',[64 64 64]/255,'linewidt',2), hold on
     plot(mean(VFRJp2(1,1,:),3),mean(FRJp2(1,1,:),3),'k*','color',[64 64 64]/255,'linewidt',2), hold on
     plot(mean(VFRJp3(1,:,:),3),mean(FRJp3(1,:,:),3),'color',[128 128 128]/255,'linewidt',2), hold on
     plot(mean(VFRJp3(1,1,:),3),mean(FRJp3(1,1,:),3),'b*','color',[128 128 128]/255,'linewidt',2), hold on
     plot(mean(VFRJp4(1,:,:),3),mean(FRJp4(1,:,:),3),'color',[192 192 192]/255,'linewidt',2), hold on
     plot(mean(VFRJp4(1,1,:),3),mean(FRJp4(1,1,:),3),'g*','color',[192 192 192]/255,'linewidt',2), hold on
     plot(x_des,x_des,'k--')
     axis([-8 8 -8 8])
     axis equal
     

     figure(42)
     plot(mean(VFRJp1(2,:,:),3),mean(FRJp1(2,:,:),3),'k','linewidt',2), hold on
     plot(mean(VFRJp1(2,1,:),3),mean(FRJp1(2,1,:),3),'k*','linewidt',2), hold on
     plot(mean(VFRJp2(2,:,:),3),mean(FRJp2(2,:,:),3),'color',[64 64 64]/255,'linewidt',2), hold on
     plot(mean(VFRJp2(2,1,:),3),mean(FRJp2(2,1,:),3),'r*','color',[64 64 64]/255,'linewidt',2), hold on
     plot(mean(VFRJp3(2,:,:),3),mean(FRJp3(2,:,:),3),'color',[128 128 128]/255,'linewidt',2), hold on
     plot(mean(VFRJp3(2,1,:),3),mean(FRJp3(2,1,:),3),'b*','color',[128 128 128]/255,'linewidt',2), hold on
     plot(mean(VFRJp4(2,:,:),3),mean(FRJp4(2,:,:),3),'color',[192 192 192]/255,'linewidt',2), hold on
     plot(mean(VFRJp4(2,1,:),3),mean(FRJp4(2,1,:),3),'g*','color',[192 192 192]/255,'linewidt',2), hold on
     plot(x_des,x_des,'k--')
     axis([-5 5 -5 5])
     axis equal



figure(8)
subplot(221)
plot(xan,mean(VFRJp1(1,:,:),3),'b','linewidt',2), hold on
plot(xan,mean(VFRJp1(2,:,:),3),'r','linewidt',2)
plot(xan,mean(FRJp1(1,:,:),3),'c','linewidt',2)
plot(xan,mean(FRJp1(2,:,:),3),'m','linewidt',2)
legend('SH vis', 'EL vis', 'SH rep', 'EL rep')
title('Position 1')
xlabel('Force direction')
ylabel('Torque')
axis([0 360 -7 7])

subplot(222)
plot(xan,mean(VFRJp2(1,:,:),3),'b','linewidt',2), hold on
plot(xan,mean(VFRJp2(2,:,:),3),'r','linewidt',2)
plot(xan,mean(FRJp2(1,:,:),3),'c','linewidt',2)
plot(xan,mean(FRJp2(2,:,:),3),'m','linewidt',2)
title('Position 2')
xlabel('Force direction')
ylabel('Torque')
axis([0 360 -7 7])

subplot(223)
plot(xan,mean(VFRJp3(1,:,:),3),'b','linewidt',2), hold on
plot(xan,mean(VFRJp3(2,:,:),3),'r','linewidt',2)
plot(xan,mean(FRJp3(1,:,:),3),'c','linewidt',2)
plot(xan,mean(FRJp3(2,:,:),3),'m','linewidt',2)
title('Position 3')
xlabel('Force direction')
ylabel('Torque')
axis([0 360 -7 7])

subplot(224)
plot(xan,mean(VFRJp4(1,:,:),3),'b','linewidt',2), hold on
plot(xan,mean(VFRJp4(2,:,:),3),'r','linewidt',2)
plot(xan,mean(FRJp4(1,:,:),3),'c','linewidt',2)
plot(xan,mean(FRJp4(2,:,:),3),'m','linewidt',2)
title('Position 4')
xlabel('Force direction')
ylabel('Torque')
axis([0 360 -7 7])
break
%%
% global Jinfo
options=optimset('tolfun',1e-10,'diffminchange',1e-5,'maxfunevals',20000,'maxiter',1000,'display','iter','tolx',1e-5);
a = 1.1;
b = 0.1;
c = 1.2;
% Ps = [1.3 0.2 1.4];
x0 = [1.1 0 1.4 0.01 0.01];
lb = [1 -2 1 0 0]%./x0;
ub = [2 2 2 1 1]%./x0;
% x0 = [1.1 1.4 0.01 0.01];
% lb = [1 1 0 0 ]./x0;
% ub = [2 2 1 1]./x0;

% e = erfit(a,b,c,Jinfo,ELf)
P = lsqnonlin('erfit',x0,[],[],options,ELf);
% P = lsqnonlin('erfitF',x0,[],[],options,Xf);
%% Model
% close figure 1
% close figure 2

for k = 1:4;
SH1 = mean(Jinfo(3,:));
SH2 = mean(Jinfo(4,:));
EL1 = mean(Jinfo(5,:));
EL2 = mean(Jinfo(6,:));
WR = 0;
L1 = mean(Jinfo(1,:));
L2 = mean(Jinfo(2,:));
L3 = 5;


 n1 = (0:0.0175:2*pi) + phiELf(k);
    ind = find(n1 > 2*pi);
    n = n1;
    n(ind) = n1(ind) - 2*pi;
    
    tar = 10 * [cos(n); sin(n)];
    
    
% n = 0:0.0175:2*pi;
% tar = 10 * [cos(n); sin(n)];

if k == 1
%     J = [-sin(SH1)*L1-sin(SH1+EL1)*L2,     -sin(SH1+EL1)*L2;
%         cos(SH1)*L1+cos(SH1+EL1)*L2,      cos(SH1+EL1)*L2];
    
    J = [-sin(SH1)*l1-sin(SH1+EL1)*l2-sin(SH1+EL1+WR)*l3,     -sin(SH1+EL1)*l2-sin(SH1+EL1+WR)*l3 -sin(SH1+EL1+WR)*l3;
          cos(SH1)*l1+cos(SH1+EL1)*l2+cos(SH1+EL1+WR)*l3,      cos(SH1+EL1)*l2+cos(SH1+EL1+WR)*l3 cos(SH1+EL1+WR)*l3];
elseif k == 2
%     J = [-sin(SH1)*L1-sin(SH1+EL2)*L2,     -sin(SH1+EL2)*L2;
%         cos(SH1)*L1+cos(SH1+EL2)*L2,      cos(SH1+EL2)*L2];
    
      J = [-sin(SH1)*l1-sin(SH1+EL2)*l2-sin(SH1+EL2+WR)*l3,     -sin(SH1+EL2)*l2-sin(SH1+EL2+WR)*l3 -sin(SH1+EL2+WR)*l3;
          cos(SH1)*l1+cos(SH1+EL2)*l2+cos(SH1+EL2+WR)*l3,      cos(SH1+EL2)*l2+cos(SH1+EL2+WR)*l3 cos(SH1+EL2+WR)*l3];
elseif k == 3
%     J = [-sin(SH2)*L1-sin(SH2+EL1)*L2,     -sin(SH2+EL1)*L2;
%         cos(SH2)*L1+cos(SH2+EL1)*L2,      cos(SH2+EL1)*L2];
    
    J = [-sin(SH2)*l1-sin(SH2+EL1)*l2-sin(SH2+EL1+WR)*l3,     -sin(SH2+EL1)*l2-sin(SH2+EL1+WR)*l3 -sin(SH2+EL1+WR)*l3;
          cos(SH2)*l1+cos(SH2+EL1)*l2+cos(SH2+EL1+WR)*l3,      cos(SH2+EL1)*l2+cos(SH2+EL1+WR)*l3 cos(SH2+EL1+WR)*l3];
elseif k == 4
%     J = [-sin(SH2)*L1-sin(SH2+EL2)*L2,     -sin(SH2+EL2)*L2;
%         cos(SH2)*L1+cos(SH2+EL2)*L2,      cos(SH2+EL2)*L2];
    
    J = [-sin(SH2)*l1-sin(SH2+EL2)*l2-sin(SH2+EL2+WR)*l3,     -sin(SH2+EL2)*l2-sin(SH2+EL2+WR)*l3 -sin(SH2+EL2+WR)*l3;
          cos(SH2)*l1+cos(SH2+EL2)*l2+cos(SH2+EL2+WR)*l3,      cos(SH2+EL2)*l2+cos(SH2+EL2+WR)*l3 cos(SH2+EL2+WR)*l3];
end
  

Joint = J' * tar;  

% SS = mean([mean(dSS_p1(1,:)) mean(dSS_p2(1,:)) mean(dSS_p4(1,:))]);
% EE = mean([mean(dEE_p1(1,:)) mean(dEE_p2(1,:)) mean(dEE_p4(1,:))]);
% SO = mean([mean(dSS_p1(2,:)) mean(dSS_p2(2,:)) mean(dSS_p4(2,:))]);
% EO = mean([mean(dEE_p1(2,:)) mean(dEE_p2(2,:)) mean(dEE_p4(2,:))]);
% SS = mean(dSS_p1(1,:));
% EE = mean(dEE_p1(1,:));
% SO = mean(dSS_p1(2,:));
% EO = mean(dEE_p1(2,:));
% % % SS = P(1);
% % % EE = P(3);
% % % Bi = P(2);
% % % SO = P(4);
% % % EO = P(5);
% SS = P(1);
% EE = P(4);
% Bi1 = 0;% P(2);
% Bi2 = 0%;P(3);
% SO = P(5);
% EO = P(6);
% SO = P(3);
% EO = P(4);
SS = 1.22;
EE = 1.40;
Bi = -0.15;
SO = 0.13;
EO = 0.12;

%estimate reproduction ellipse
Joint_rep(1,:) = SS * Joint(1,:) + SO;
Joint_rep(2,:) = EE * Joint(2,:) + EO;

figure(101)
subplot(121)
plot(Joint(1,:),Joint_rep(1,:)),hold on
subplot(122)
plot(Joint(2,:),Joint_rep(2,:)), hold on
%with bi-articular muscles
Joint_rep2 = [SS   Bi; Bi EE ] * Joint;
Joint_rep2(1,:) = Joint_rep2(1,:) + SO;
Joint_rep2(2,:) = Joint_rep2(2,:) + EO;
rep_el = inv(J') * Joint_rep2;

VAFx(k,:) = VAFs(ELf(1,:,k),rep_el(1,:));
VAFy(k,:) = VAFs(ELf(2,:,k),rep_el(2,:));


magEL = sqrt(ELf(1,:,k).^2 + ELf(2,:,k).^2);
magSIM = sqrt(rep_el(1,:).^2 + rep_el(2,:).^2);


VAFug(k,:) = (1-abs(magEL-magSIM) / magEL)*100;

JrSS = [SS 0; 0 0] * Joint;
JrEE = [0 0; 0 EE] * Joint;
JrBB = [Bi Bi; Bi Bi] * Joint;

figure(101)
subplot(121)
plot(Joint(1,:),Joint_rep2(1,:),'r')
subplot(122)
plot(Joint(2,:),Joint_rep2(2,:),'r')

%directions of mono- and bi-articular muscles
sdir = inv(J') * [1 0; 0 0];
edir = inv(J') * [0 0; 0 1];
bdir = inv(J') * [1 1; 1 1];
MS = sdir * [cos(n); sin(n)];
ME = edir * [cos(n); sin(n)];
BB = bdir * [cos(n); sin(n)];

figure(20)
subplot(2,2,k)
plot(tar(1,:),tar(2,:)), hold on
plot(rep_el(1,:),rep_el(2,:),'c','linewidt',2), hold on
% plot(MS(1,:), MS(2,:),'c')
% plot(ME(1,:), ME(2,:),'m')
% plot(BB(1,:), BB(2,:),'g')
axis([-20 20 -20 20])
axis equal

figure(7)
subplot(2,2,k)
polar(n,abs(Joint(1,:))),hold on
polar(n,abs(Joint(2,:)))
polar(n,abs(Joint_rep(1,:)),'r')
polar(n,abs(Joint_rep(2,:)),'r')
polar(n,abs(JrSS(1,:)),'c')
polar(n,abs(JrEE(2,:)),'m')
polar(n,abs(JrBB(1,:)),'k')
end

