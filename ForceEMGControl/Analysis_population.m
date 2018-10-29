%% Data Analysis for population
% clear all;
addpath(genpath('Tools'));

switch computer
    case 'PCWIN64'
        filepath =  ['D:\Student_experiments\Virginia\Data\TD\'];
    case 'MACI64'
        filepath =  ['/Users/virginia/Documents/MATLAB/Thesis/Data/TD/'];
end

Aparams_pp = [];
trial_pp_force = [];
trial_pp_EMG = [];

dirs = dir(filepath);
files = dirs(~[dirs.isdir]);
for i = 1:length(files)
    load([filepath,files(i).name]);
    Aparams_pp = [Aparams_pp, trial_pp.Aparams];
    trial_pp_force = [trial_pp_force, trial_pp.forceCO];
    trial_pp_EMG = [trial_pp_EMG, trial_pp.EMGCO];
end

%% FIGURES
fields_pp = {'force','EMG.rect'};
trial_avg_force = ppAngleAvg(trial_pp_force, fields_pp, Aparams_pp);
trial_avg_EMG = ppAngleAvg(trial_pp_EMG, fields_pp, Aparams_pp);

targetAnglesForce = unique([trial_avg_force.angle]);
targetAnglesEMG = unique([trial_avg_EMG.angle]);
angComp = targetAnglesForce(ismember(targetAnglesForce,targetAnglesEMG));


%% Force mean - Plot
figure('Name','Force Magnitude Mean');
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
h1 = plot(rad2deg([trial_avg_force.angle]),[trial_avg_force.force_mean],'bo');
hold on;
errorbar(rad2deg([trial_avg_force.angle]),[trial_avg_force.force_mean],[trial_avg_force.force_std],'b.');
h2 = plot(rad2deg([trial_avg_EMG.angle]),[trial_avg_EMG.force_mean],'ro');
errorbar(rad2deg([trial_avg_EMG.angle]),[trial_avg_EMG.force_mean],[trial_avg_EMG.force_std],'r.');
set(gca,'XTick',rad2deg(angComp))
ylabel('Force [N]'); xlabel('Target');
title('Mean Force Magnitude');
legend([h1,h2],{'ForceCO','EMGCO'});

%% Force mean - Polar plot
th_fill = targetAnglesForce;
rf_fill_up = [trial_avg_force.force_mean]+[trial_avg_force.force_std];
rf_fill_do = [trial_avg_force.force_mean]-[trial_avg_force.force_std];
[xf_fill_up,yf_fill_up] = pol2cart(th_fill,rf_fill_up);
[xf_fill_do,yf_fill_do] = pol2cart(th_fill,rf_fill_do);

th_fill = targetAnglesEMG;
re_fill_up = [trial_avg_EMG.force_mean]+[trial_avg_EMG.force_std];
re_fill_do = [trial_avg_EMG.force_mean]-[trial_avg_EMG.force_std];
[xe_fill_up,ye_fill_up] = pol2cart(th_fill,re_fill_up);
[xe_fill_do,ye_fill_do] = pol2cart(th_fill,re_fill_do);

figure('Name','Force Magnitude Mean');
polar(0,max([yf_fill_up,ye_fill_up]),'w'); hold on;
f = polar(targetAnglesForce,[trial_avg_force.force_mean],'b-.o');
f.MarkerFaceColor = 'b';
hold on;
e = polar(targetAnglesEMG,[trial_avg_EMG.force_mean],'r-.o');
e.MarkerFaceColor = 'r';
fill([xf_fill_up fliplr(xf_fill_do)],[yf_fill_up fliplr(yf_fill_do)],'b','edgecolor','none','facealpha',0.3);
fill([xe_fill_up fliplr(xe_fill_do)],[ye_fill_up fliplr(ye_fill_do)],'r','edgecolor','none','facealpha',0.3);
% xticks(rad2deg(angComp)); %thetaticklabels(Aparams.muscComp);
% xticklabels(rad2deg(angComp)); 
% legend([f,e],{'ForceCO','EMGCO'},'Location','bestoutside');
set(gca,'FontSize',12);

%% Mean force polarplot
figure('Name','Force Magnitude Mean');
f = polarplot(targetAnglesForce,[trial_avg_force.force_mean],'b-.o');
f.MarkerFaceColor = 'b'; f.MarkerSize = 8;
hold on;
e = polarplot(targetAnglesEMG,[trial_avg_EMG.force_mean],'r-.o');
e.MarkerFaceColor = 'r'; e.MarkerSize = 8;
polarplot(targetAnglesForce,rf_fill_up,'b-','Color',[0 0 1 0.4]);
polarplot(targetAnglesForce,rf_fill_do,'b-','Color',[0 0 1 0.4]);
polarplot(targetAnglesEMG,re_fill_up,'r-','Color',[1 0 0 0.4]);
polarplot(targetAnglesEMG,re_fill_do,'r-','Color',[1 0 0 0.4]);
thetaticks(rad2deg(angComp));
thetaticklabels(rad2deg(angComp));
% xticklabels(rad2deg(angComp)); 
% legend([f,e],{'ForceCO','EMGCO'},'Location','bestoutside');
set(gca,'FontSize',12);

%%
    h=polar(repmat(pdData.(move_cor),2,1),maxRadius*[0;1],linspec);
    set(h,'linewidth',2,'color',color)
    th_fill = [pdData.([move_cor, 'CI'])(2) pdData.(move_cor) pdData.([move_cor, 'CI'])(1) 0];
    r_fill = [maxRadius maxRadius maxRadius 0];
    [x_fill,y_fill] = pol2cart(th_fill,r_fill);
    patch(x_fill,y_fill,color,'edgecolor','none','facealpha',0.3);

    
yu = y+.1;
yl = y-.1;
fill([x;flipud(x)],[y-dy;flipud(y+dy)]

fill([x fliplr(x)], [yu fliplr(yl)], [.9 .9 .9], 'linestyle', 'none')


