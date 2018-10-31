%% Data Analysis for population
clear all;
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
fields_pp = {'force.filtmag_mean','EMG.rect'};
trial_avg_force = ppAngleAvg(trial_pp_force, fields_pp, Aparams_pp);
trial_avg_EMG = ppAngleAvg(trial_pp_EMG, fields_pp, Aparams_pp);

targetAnglesForce = unique([trial_avg_force.angle]);
targetAnglesEMG = unique([trial_avg_EMG.angle]);
angComp = targetAnglesForce(ismember(targetAnglesForce,targetAnglesEMG));

avgforce_force = [trial_avg_force.force];
avgEMG_force = [trial_avg_EMG.force];

avgforce_EMG = [trial_avg_force.EMG];
avgEMG_EMG = [trial_avg_EMG.EMG];

avgCV = CV(trial_pp_force,trial_pp_EMG,angComp);

[xf,yf] = getBound(targetAnglesForce,[avgforce_force.filtmag_mean_mean],[avgforce_force.filtmag_mean_std],'cart');
[xE,yE] = getBound(targetAnglesEMG,[avgEMG_force.filtmag_mean_mean],[avgEMG_force.filtmag_mean_std],'cart');

[thf,rf] = getBound(targetAnglesForce,[avgforce_force.filtmag_mean_mean],[avgforce_force.filtmag_mean_std],'polar');
[thE,rE] = getBound(targetAnglesEMG,[avgEMG_force.filtmag_mean_mean],[avgEMG_force.filtmag_mean_std],'polar');


%% Force mean - Plot
figure('Name','Force Magnitude Mean');
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
h1 = plot(rad2deg([trial_avg_force.angle]),[avgforce_force.filtmag_mean_mean],'bo');
hold on;
errorbar(rad2deg([trial_avg_force.angle]),[avgforce_force.filtmag_mean_mean],[avgforce_force.filtmag_mean_std],'b.');
h2 = plot(rad2deg([trial_avg_EMG.angle]),[avgEMG_force.filtmag_mean_mean],'ro');
errorbar(rad2deg([trial_avg_EMG.angle]),[avgEMG_force.filtmag_mean_mean],[avgEMG_force.filtmag_mean_std],'r.');
set(gca,'XTick',rad2deg(angComp));
ylabel('Force [N]'); xlabel('Target');
title('Mean Force Magnitude');
legend([h1,h2],{'ForceCO','EMGCO'});

%% Force mean - Polar plot
figure('Name','Force Magnitude Mean');
polar(0,max([yf,yE]),'w'); hold on;
f = polar(targetAnglesForce,[avgforce_force.filtmag_mean_mean],'b-.o');
f.MarkerFaceColor = 'b';
hold on;
e = polar(targetAnglesEMG,[avgEMG_force.filtmag_mean_mean],'r-.o');
e.MarkerFaceColor = 'r';
fill(xf,yf,'b','edgecolor','none','facealpha',0.3);
fill(xE,yE,'r','edgecolor','none','facealpha',0.3);
% xticks(rad2deg(angComp)); %thetaticklabels(Aparams.muscComp);
% xticklabels(rad2deg(angComp)); 
legend([f,e],{'ForceCO','EMGCO'},'Location','bestoutside');
set(gca,'FontSize',12);

%% Mean force polarplot
figure('Name','Force Magnitude Mean');
f = polarplot(targetAnglesForce,[avgforce_force.filtmag_mean_mean],'b-.o');
f.MarkerFaceColor = 'b'; f.MarkerSize = 8;
hold on;
e = polarplot(targetAnglesEMG,[avgEMG_force.filtmag_mean_mean],'r-.o');
e.MarkerFaceColor = 'r'; e.MarkerSize = 8;
polarplot(targetAnglesForce,rf(1,:),'b-','Color',[0 0 1 0.4]);
polarplot(targetAnglesForce,rf(2,:),'b-','Color',[0 0 1 0.4]);
polarplot(targetAnglesEMG,rE(1,:),'r-','Color',[1 0 0 0.4]);
polarplot(targetAnglesEMG,rE(2,:),'r-','Color',[1 0 0 0.4]);
thetaticks(rad2deg(angComp));
thetaticklabels(rad2deg(angComp));
legend([f,e],{'ForceCO','EMGCO'},'Location','bestoutside');
set(gca,'FontSize',12);

%% Plot mean force magnitude for all participants
cols = hsv(length(trial_pp_force));
figure('Name','Force Magnitude Mean');
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(trial_pp_force)
    scatter(rad2deg(trial_pp_force(i).angles),trial_pp_force(i).force.filtmag_mean,60,cols(i,:));
    hold on;
    errorbar(rad2deg(trial_pp_force(i).angles),trial_pp_force(i).force.filtmag_mean,trial_pp_force(i).force.filtmag_pstd,'.','Color',cols(i,:));
    scatter(rad2deg(trial_pp_EMG(i).angles),trial_pp_EMG(i).force.filtmag_mean,60,cols(i,:),'filled');
    errorbar(rad2deg(trial_pp_EMG(i).angles),trial_pp_EMG(i).force.filtmag_mean,trial_pp_EMG(i).force.filtmag_pstd,'.','Color',cols(i,:));
end
plot(rad2deg([trial_avg_force.angle]),[avgforce_force.filtmag_mean_mean],'k');
plot(rad2deg([trial_avg_EMG.angle]),[avgEMG_force.filtmag_mean_mean],'k-.');
set(gca,'XTick',rad2deg(angComp))
ylabel('Force [N]'); xlabel('Target [deg]');
title('Mean Force Magnitude');
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko-');
h(2) = plot(NaN,NaN,'ko-.','MarkerFaceColor','k');
legend(h,{'ForceCO','EMGCO'});
set(gca,'FontSize',12);

%% Force CV - Plot
figure('Name','Force Magnitude Mean');
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
h1 = plot(rad2deg(avgCV.angle),avgCV.force.mean,'bo');
hold on;
errorbar(rad2deg(avgCV.angle),avgCV.force.mean,avgCV.force.sem,'bo');
h2 = plot(rad2deg(avgCV.angle),avgCV.EMG.mean,'ro');
errorbar(rad2deg(avgCV.angle),avgCV.EMG.mean,avgCV.EMG.sem,'r.');
set(gca,'XTick',rad2deg(angComp));
ylabel('CV [%]'); xlabel('Target');
title('Mean Force Magnitude');
legend([h1,h2],{'ForceCO','EMGCO'});


