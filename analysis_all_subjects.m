%% Data Analysis for all subjects/population
clear all;
addpath(genpath('Tools'));

switch computer
    case 'PCWIN64'
        filepath =  ['D:\Student_experiments\Virginia\Data\TD\'];
    case 'MACI64'
        filepath =  ['/Users/virginia/Documents/MATLAB/Thesis/Data/TD/UB/'];
end

Aparams_pp = [];
trial_pp_force = [];
trial_pp_EMG = [];

dirs = dir(filepath);
files = dirs(~[dirs.isdir]);
for i = 1:length(files)
    if ~isempty(strfind(files(i).name,'_s'))
        load([filepath,files(i).name]);
        if isempty(strfind(files(i).name,'s04'))
            Aparams_pp = [Aparams_pp, trial_pp.Aparams];
            trial_pp_force = [trial_pp_force, trial_pp.forceCO];
            trial_pp_EMG = [trial_pp_EMG, trial_pp.EMGCO];
        end
    end
end

%% Analysis pre-figures
fields_avg = {'force.filt_mean','force.mag_mean','EMG.rect'};
trial_avg_force = ppAngleAvg(trial_pp_force, fields_avg, Aparams_pp);
trial_avg_EMG = ppAngleAvg(trial_pp_EMG, fields_avg, Aparams_pp);

targetAnglesForce = unique([trial_avg_force.angle]);
targetAnglesEMG = unique([trial_avg_EMG.angle]);
angComp = targetAnglesForce(ismember(targetAnglesForce,targetAnglesEMG));
degAngComp = rad2deg(angComp);

avgforce_force = [trial_avg_force.force];
avgEMG_force = [trial_avg_EMG.force];

avgforce_EMG = [trial_avg_force.EMG];
avgEMG_EMG = [trial_avg_EMG.EMG];

avgCV = meanppCV(trial_pp_force,trial_pp_EMG,angComp);
avgEMGerr = meanppEMGerr(trial_pp_force,trial_pp_EMG,Aparams_pp,angComp,'abs');
avgCoh = meanppCoh(trial_pp_force,trial_pp_EMG,Aparams_pp,angComp);

[xf,yf] = getBound(targetAnglesForce,[avgforce_force.mag_mean_mean],[avgforce_force.mag_mean_std],'cart');
[xE,yE] = getBound(targetAnglesEMG,[avgEMG_force.mag_mean_mean],[avgEMG_force.mag_mean_std],'cart');

[thf,rf] = getBound(targetAnglesForce,[avgforce_force.mag_mean_mean],[avgforce_force.mag_mean_std],'polar');
[thE,rE] = getBound(targetAnglesEMG,[avgEMG_force.mag_mean_mean],[avgEMG_force.mag_mean_std],'polar');

[~,fit_force_fmag_mean] = fitCurve(rad2deg([trial_avg_force.angle]),[avgforce_force.mag_mean_mean],3);
[~,fit_EMG_fmag_mean] = fitCurve(rad2deg([trial_avg_EMG.angle]),[avgEMG_force.mag_mean_mean],3);

[~,fit_force_fCV] = fitCurve(rad2deg(avgCV.angle),avgCV.force.mean,3);
[fit_angle,fit_EMG_fCV] = fitCurve(rad2deg(avgCV.angle),avgCV.EMG.mean,3);

force_EMG_mean = reshape([avgforce_EMG.rect_mean],length(avgforce_EMG(1).rect_mean),length([trial_avg_force.angle]))';
EMG_EMG_mean = reshape([avgEMG_EMG.rect_mean],length(avgEMG_EMG(1).rect_mean),length([trial_avg_EMG.angle]))';
force_EMG_sem = reshape([avgforce_EMG.rect_sem],length(avgforce_EMG(1).rect_sem),length([trial_avg_force.angle]))';
EMG_EMG_sem = reshape([avgEMG_EMG.rect_sem],length(avgEMG_EMG(1).rect_sem),length([trial_avg_EMG.angle]))';
fit_force_EMG_mean = [];
fit_EMG_EMG_mean = [];
for i = 1:length(avgforce_EMG(1).rect_mean)
    [~,fit_force_EMG_mean(:,i)] = fitCurve(rad2deg([trial_avg_force.angle]),force_EMG_mean(:,i)',3);
    [~,fit_EMG_EMG_mean(:,i)] = fitCurve(rad2deg([trial_avg_EMG.angle]),EMG_EMG_mean(:,i)',3);
end

force_EMG_smean = reshape([avgforce_EMG.rect_scale_mean],length(avgforce_EMG(1).rect_scale_mean),length([trial_avg_force.angle]))';
EMG_EMG_smean = reshape([avgEMG_EMG.rect_scale_mean],length(avgEMG_EMG(1).rect_scale_mean),length([trial_avg_EMG.angle]))';
force_EMG_ssem = reshape([avgforce_EMG.rect_scale_sem],length(avgforce_EMG(1).rect_scale_sem),length([trial_avg_force.angle]))';
EMG_EMG_ssem = reshape([avgEMG_EMG.rect_scale_sem],length(avgEMG_EMG(1).rect_scale_sem),length([trial_avg_EMG.angle]))';
fit_force_EMG_smean = [];
fit_EMG_EMG_smean = [];
for i = 1:length(avgforce_EMG(1).rect_scale_mean)
    [~,fit_force_EMG_smean(:,i)] = fitCurve(degAngComp,force_EMG_smean(:,i)',3);
    [~,fit_EMG_EMG_smean(:,i)] = fitCurve(degAngComp,EMG_EMG_smean(:,i)',3);
end

fc = 80;

%% Statistical analysis

fields_stat = {'force.mag_mean','EMG.rect','trial_coh.rect.asig_coh','trial_coh.rect.nasig_coh',...
    'trial_coh.rect.asig_z','trial_coh.rect.nasig_z'};
ppStats = statAnalysis(trial_pp_force,trial_pp_EMG,Aparams_pp,angComp,fields_stat);

%% Force
%% Force mean - Plot
figure('Name','Force Magnitude Mean');
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
h1 = plot(degAngComp,[avgforce_force.mag_mean_mean],'bo','MarkerFaceColor','b');
hold on;
errorbar(degAngComp,[avgforce_force.mag_mean_mean],[avgforce_force.mag_mean_std],'b.');
h2 = plot(degAngComp,[avgEMG_force.mag_mean_mean],'ro','MarkerFaceColor','r');
errorbar(degAngComp,[avgEMG_force.mag_mean_mean],[avgEMG_force.mag_mean_std],'r.');
plot(fit_angle,fit_force_fmag_mean,'b--');
plot(fit_angle,fit_EMG_fmag_mean,'r--')
set(gca,'XTick',rad2deg(angComp));
ylabel('Force [N]'); xlabel('Target [deg]');
title('Mean Force Magnitude');
legend([h1,h2],{'ForceCO','EMGCO'});
set(gca,'FontSize',12);

%% Force trajectory mean - Plot
figure('Name','Force Trajectory Mean');
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(angComp)
    h1 = plot(avgforce_force(i).filt_mean_mean(1),avgforce_force(i).filt_mean_mean(2),'bo','MarkerFaceColor','b');
    hold on;
    %errorbar(degAngComp,[avgforce_force.filt_mean_mean],[avgforce_force.mag_mean_std],'b.');
    h2 = plot(avgEMG_force(i).filt_mean_mean(1),avgEMG_force(i).filt_mean_mean(2),'ro','MarkerFaceColor','r');
    %errorbar(degAngComp,[avgEMG_force.mag_mean_mean],[avgEMG_force.mag_mean_std],'r.');
    text(avgforce_force(i).filt_mean_mean(1),avgforce_force(i).filt_mean_mean(2),num2str(degAngComp(i)),'VerticalAlignment','top','HorizontalAlignment','left');
    text(avgEMG_force(i).filt_mean_mean(1),avgEMG_force(i).filt_mean_mean(2),num2str(degAngComp(i)),'VerticalAlignment','top','HorizontalAlignment','left');
end
axis square;
grid on;
xlim([-30 30]); ylim([-30 30]); 
ylabel('Fy [N]'); xlabel('Fx [N]');
title('Mean Force Magnitude');
legend([h1,h2],{'ForceCO','EMGCO'});
set(gca,'FontSize',12);

%% Force mean - Polar plot
figure('Name','Force Magnitude Mean');
g = polar(0,max([yf,yE]),'w'); 
polarticks(8,g); hold on;
f = polar(angComp,[avgforce_force.mag_mean_mean],'b-.o');
f.MarkerFaceColor = 'b';
hold on;
e = polar(angComp,[avgEMG_force.mag_mean_mean],'r-.o');
e.MarkerFaceColor = 'r';
fill(xf,yf,'b','edgecolor','none','facealpha',0.3);
fill(xE,yE,'r','edgecolor','none','facealpha',0.3);
title('Force Magnitude Mean')
% xticks(rad2deg(angComp)); %thetaticklabels(Aparams.muscComp);
% xticklabels(rad2deg(angComp));
legend([f,e],{'ForceCO','EMGCO'},'Location','bestoutside');
set(gca,'FontSize',12);

%% Mean force polarplot
figure('Name','Force Magnitude Mean');
f = polarplot(angComp,[avgforce_force.mag_mean_mean],'b-.o');
f.MarkerFaceColor = 'b'; f.MarkerSize = 8;
hold on;
e = polarplot(angComp,[avgEMG_force.mag_mean_mean],'r-.o');
e.MarkerFaceColor = 'r'; e.MarkerSize = 8;
polarplot(angComp,rf(1,:),'b-','Color',[0 0 1 0.4]);
polarplot(angComp,rf(2,:),'b-','Color',[0 0 1 0.4]);
polarplot(angComp,rE(1,:),'r-','Color',[1 0 0 0.4]);
polarplot(angComp,rE(2,:),'r-','Color',[1 0 0 0.4]);
thetaticks(degAngComp);
thetaticklabels(degAngComp);
legend([f,e],{'ForceCO','EMGCO'},'Location','bestoutside');
set(gca,'FontSize',12);

%% Plot mean force magnitude for all participants
cols = hsv(length(trial_pp_force));
figure('Name','Force Magnitude Mean');
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(trial_pp_force)
    scatter(rad2deg(trial_pp_force(i).angle),trial_pp_force(i).force.mag_mean,60,cols(i,:));
    hold on;
    errorbar(rad2deg(trial_pp_force(i).angle),trial_pp_force(i).force.mag_mean,trial_pp_force(i).force.mag_pstd,'.','Color',cols(i,:));
    scatter(rad2deg(trial_pp_EMG(i).angle),trial_pp_EMG(i).force.mag_mean,60,cols(i,:),'filled');
    errorbar(rad2deg(trial_pp_EMG(i).angle),trial_pp_EMG(i).force.mag_mean,trial_pp_EMG(i).force.mag_pstd,'.','Color',cols(i,:));
end
plot(fit_angle,fit_force_fmag_mean,'k');
plot(fit_angle,fit_EMG_fmag_mean,'k-.');
set(gca,'XTick',rad2deg(angComp))
ylabel('Force [N]'); xlabel('Target [deg]');
title('Mean Force Magnitude');
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko-');
h(2) = plot(NaN,NaN,'ko-.','MarkerFaceColor','k');
legend(h,{'ForceCO','EMGCO'});
set(gca,'FontSize',12);

%% Force CV - Plot
figure('Name','Force CV');
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
h1 = plot(rad2deg(avgCV.angle),avgCV.force.mean,'bo','MarkerFaceColor','b');
hold on;
errorbar(rad2deg(avgCV.angle),avgCV.force.mean,avgCV.force.sem,'bo');
h2 = plot(rad2deg(avgCV.angle),avgCV.EMG.mean,'ro','MarkerFaceColor','r');
errorbar(rad2deg(avgCV.angle),avgCV.EMG.mean,avgCV.EMG.sem,'r.');
plot(fit_angle,fit_force_fCV,'b--');
plot(fit_angle,fit_EMG_fCV,'r--');
set(gca,'XTick',rad2deg(angComp));
ylabel('CV [%]'); xlabel('Target [deg]');
title('Mean Force CV');
legend([h1,h2],{'ForceCO','EMGCO'});
set(gca,'FontSize',12);

%% Force mean and CV - Polar
figure('Name','Force CV');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
h1 = polarscatter(angComp,[avgforce_force.mag_mean_mean],60,'filled','b');
hold on;
h2 = polarscatter(angComp,[avgEMG_force.mag_mean_mean],60,'filled','r');
errorpolar(angComp,[avgforce_force.mag_mean_mean],[avgforce_force.mag_mean_sem],'b')
errorpolar(angComp,[avgEMG_force.mag_mean_mean],[avgEMG_force.mag_mean_sem],'r')
thetaticks(degAngComp); thetaticklabels(degAngComp);
title({'Force Magnitude'; 'Mean'})
set(gca,'FontSize',18);

subplot(1,2,2)
h1 = polarscatter(avgCV.angle,avgCV.force.mean,60,'filled','b');
hold on;
h2 = polarscatter(avgCV.angle,avgCV.EMG.mean,60,'filled','r');
errorpolar(avgCV.angle,avgCV.force.mean,avgCV.force.sem,'b')
errorpolar(avgCV.angle,avgCV.EMG.mean,avgCV.EMG.sem,'r')
thetaticks(degAngComp); thetaticklabels(degAngComp);
title({'Force Magnitude'; 'CV'})
set(gca,'FontSize',18);
legend([h1,h2],'FC','MC','Location','bestoutside');

%% EMG
%% EMG mean
figure('Name','EMG Mean');
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for i = 1:size(force_EMG_mean,2)
    subplot(size(force_EMG_mean,2),1,i)
    h1 = plot(degAngComp,force_EMG_mean(:,i),'bo','MarkerFaceColor','b');
    hold on;
    h2 = plot(degAngComp,EMG_EMG_mean(:,i),'ro','MarkerFaceColor','r');
    errorbar(degAngComp,force_EMG_mean(:,i),force_EMG_sem(:,i),'b.')
    errorbar(degAngComp,EMG_EMG_mean(:,i),EMG_EMG_sem(:,i),'r.')
    plot(fit_angle,fit_force_EMG_mean(:,i),'b--');
    plot(fit_angle,fit_EMG_EMG_mean(:,i),'r--')
    set(gca,'XTick',degAngComp);
    if i == size(force_EMG_mean,2)
        xlabel('Target [deg]');
    end
    ylabel('EMG [-]');
    title(Aparams_pp(1).chanControlName{i});
    legend([h1,h2],{'ForceCO','EMGCO'});
    set(gca,'FontSize',12);
end

%% Scaled EMG mean
figure('Name','Scaled EMG Mean');
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for i = 1:size(force_EMG_smean,2)
    subplot(size(force_EMG_smean,2),1,i)
    h1 = plot(degAngComp,force_EMG_smean(:,i),'bo','MarkerFaceColor','b');
    hold on;
    h2 = plot(degAngComp,EMG_EMG_smean(:,i),'ro','MarkerFaceColor','r');
    errorbar(degAngComp,force_EMG_smean(:,i),force_EMG_ssem(:,i),'b.','Capsize',12)
    errorbar(degAngComp,EMG_EMG_smean(:,i),EMG_EMG_ssem(:,i),'r.','Capsize',12)
    ylim([0 1.2]); xlim([degAngComp(1)-10 degAngComp(end)+10])
    line(xlim,[1 1],'LineStyle',':','Color','k','Linewidth',1.7);
    line(xlim,[0.7 0.7],'LineStyle',':','Color','k','Linewidth',1.7);
    set(gca,'XTick',degAngComp);
    if i == size(force_EMG_smean,2)
        xlabel('Target [deg]');
    end
    ylabel('EMG [-]');
    title(Aparams_pp(1).chanControlName{i});
    if i == 1
    legend([h1,h2],{'FC','MC'});
    end
    set(gca,'FontSize',13);
end

%% EMG difference mean
figure('Name','Mean EMG Difference');
h = plot(rad2deg(avgEMGerr.angle),avgEMGerr.mean,'o-.');
hold on;
for i = 1:length(h)
    h(i).MarkerFaceColor = h(i).Color;
    errorbar(rad2deg(avgEMGerr.angle),avgEMGerr.mean(:,i),avgEMGerr.sem(:,i),'.','Color',h(i).Color)
end
xticks(rad2deg(Aparams_pp(1).targetAnglesForce));
ylim([0 max(avgEMGerr.mean(:)+avgEMGerr.sem(:))+0.1]);
xlabel('Target [deg]'); ylabel('Error [-]')
legend(Aparams_pp(1).chanControlName)
title('Mean EMG error Force-EMG tasks')
set(gca,'FontSize',12);

%% Coherence
%% Mean coh - Polar only angmusc
figure('Name','Significant Coherence Mean');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
k = 0; s = [];
bands = {'alp','beta','gam'};
bands_label = {'Alpha','Beta','Gamma'};
for i = 1:length(avgCoh.musc)
    iang = ismember(avgCoh.angle,avgCoh.angmusc{i});
    yf = []; ye = [];
    for j = 1:length(bands)
        k = k+1;
        s = subplot(length(avgCoh.musc),length(bands),k);
        [xf,yf] = getBound(avgCoh.angle(iang),avgCoh.force.msig_coh.(bands{j}).mean(iang,i),avgCoh.force.msig_coh.(bands{j}).sem(iang,i),'cart');
        [xe,ye] = getBound(avgCoh.angle(iang),avgCoh.EMG.msig_coh.(bands{j}).mean(iang,i),avgCoh.EMG.msig_coh.(bands{j}).sem(iang,i),'cart');       
        g = polar(0,max(abs([yf;ye;xf;xe]))+0.1*max(abs([yf;ye;xf;xe])),'w');
        polarticks(8,g);
        hold on;
        f = polar(avgCoh.angle(iang),avgCoh.force.msig_coh.(bands{j}).mean(iang,i),'b-.o');
        f.MarkerFaceColor = 'b';
        hold on;
        e = polar(avgCoh.angle(iang),avgCoh.EMG.msig_coh.(bands{j}).mean(iang,i),'r-.o');
        e.MarkerFaceColor = 'r';
        fill(xf,yf,'b','edgecolor','none','facealpha',0.3);
        fill(xe,ye,'r','edgecolor','none','facealpha',0.3);
        if j == 1
            ylabel([avgCoh.musc{i}{1},'-',avgCoh.musc{i}{2}],'FontWeight','bold');
        end
        if i == 1
            title(bands_label{j});
        end
        set(gca,'FontSize',14);
    end
end

%% Mean coh - Polar all angles
bands = {'alp','beta','gam'};
bands_label = {'Alpha','Beta','Gamma'};
for i = 1:length(avgCoh.musc)
    figure('Name','Significant Coherence Mean');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
yf = []; ye = [];
    for j = 1:length(bands)
        subplot(1,length(bands),j)
        [xf,yf] = getBound(avgCoh.angle,avgCoh.force.msig_coh.(bands{j}).mean(:,i),avgCoh.force.msig_coh.(bands{j}).sem(:,i),'cart');
        [xe,ye] = getBound(avgCoh.angle,avgCoh.EMG.msig_coh.(bands{j}).mean(:,i),avgCoh.EMG.msig_coh.(bands{j}).sem(:,i),'cart');
        g = polar(0,max(abs([yf;ye;xf;xe])),'w');
        polarticks(8,g);        
        hold on;
        f = polar(avgCoh.angle,avgCoh.force.msig_coh.(bands{j}).mean(:,i),'b-.o');
        f.MarkerFaceColor = 'b';
        hold on;
        e = polar(avgCoh.angle,avgCoh.EMG.msig_coh.(bands{j}).mean(:,i),'r-.o');
        e.MarkerFaceColor = 'r';
        fill(xf,yf,'b','edgecolor','none','facealpha',0.3);
        fill(xe,ye,'r','edgecolor','none','facealpha',0.3);
            title(bands_label{j})
        if j == 1
            ylabel([avgCoh.musc{i}{1},'-',avgCoh.musc{i}{2}],'FontWeight','bold');
        end
        set(gca,'FontSize',14);
    end
end

%% Mean coh - Plot
figure('Name','Mean Significant Coherence');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
k = 0;
bands = {'alp','beta','gam'};
bands_label = {'Alpha','Beta','Gamma'};
for i = 1:length(avgCoh.musc)
    iang = ismember(avgCoh.angle,avgCoh.angmusc{i});
    yf = []; ye = [];
    for j = 1:length(bands)
        k = k+1;
        subplot(length(avgCoh.musc),length(bands),k)
        f = plot(rad2deg(avgCoh.angle(iang)),avgCoh.force.msig_coh.(bands{j}).mean(iang,i),'b-.o');
        f.MarkerFaceColor = 'b';
        hold on;
        errorbar(rad2deg(avgCoh.angle(iang)),avgCoh.force.msig_coh.(bands{j}).mean(iang,i),avgCoh.force.msig_coh.(bands{j}).sem(iang,i),'b.');
        e = plot(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.msig_coh.(bands{j}).mean(iang,i),'r-.o');
        e.MarkerFaceColor = 'r';
        errorbar(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.msig_coh.(bands{j}).mean(iang,i),avgCoh.EMG.msig_coh.(bands{j}).sem(iang,i),'r.');
        if i == 1
            title(bands_label{j})
        end
        if j == 1
            ylabel([avgCoh.musc{i}{1},'-',avgCoh.musc{i}{2}],'FontWeight','bold');
        end
        if i == length(avgCoh.musc)
            xlabel('Target [deg]');
        end
        set(gca,'XTick',rad2deg(avgCoh.angle(iang)));
        set(gca,'FontSize',14);
    end
end

%% Area coh - Polar all angles
for i = 1:length(avgCoh.musc)
    figure('Name','Significant Coherence Mean');
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
yf = []; ye = [];
    for j = 1:length(bands)
        subplot(1,length(bands),j)
        [xf,yf] = getBound(avgCoh.angle,avgCoh.force.asig_coh.(bands{j}).mean(:,i),avgCoh.force.asig_coh.(bands{j}).sem(:,i),'cart');
        [xe,ye] = getBound(avgCoh.angle,avgCoh.EMG.asig_coh.(bands{j}).mean(:,i),avgCoh.EMG.asig_coh.(bands{j}).sem(:,i),'cart');
        g = polar(0,max(abs([yf;ye;xf;xe]))+0.3*max(abs([yf;ye;xf;xe])),'w');
        polarticks(8,g);        
        hold on;
        f = polar(avgCoh.angle,avgCoh.force.asig_coh.(bands{j}).mean(:,i),'b-.o');
        f.MarkerFaceColor = 'b';
        hold on;
        e = polar(avgCoh.angle,avgCoh.EMG.asig_coh.(bands{j}).mean(:,i),'r-.o');
        e.MarkerFaceColor = 'r';
        fill(xf,yf,'b','edgecolor','none','facealpha',0.3);
        fill(xe,ye,'r','edgecolor','none','facealpha',0.3);
            title(bands_label{j})
        if j == 1
            ylabel([avgCoh.musc{i}{1},'-',avgCoh.musc{i}{2}],'FontWeight','bold');
        end
        set(gca,'FontSize',14);
    end
end

%% Area coh - Polar all musccomb
figure('Name','Area Coherence Mean');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
k = 0;
bands = {'alp','beta','gam'};
bands_label = {'Alpha','Beta','Gamma'};
for i = 1:length(avgCoh.musc)
    iang = ismember(avgCoh.angle,avgCoh.angmusc{i});
    yf = []; ye = [];
    for j = 1:length(bands)
        k = k+1;        
        subplot(length(avgCoh.musc),length(bands),k)
        [xf,yf] = getBound(avgCoh.angle(iang),avgCoh.force.asig_coh.(bands{j}).mean(iang,i),avgCoh.force.asig_coh.(bands{j}).sem(iang,i),'cart');
        [xe,ye] = getBound(avgCoh.angle(iang),avgCoh.EMG.asig_coh.(bands{j}).mean(iang,i),avgCoh.EMG.asig_coh.(bands{j}).sem(iang,i),'cart');
        g = polar(0,max(abs([yf;ye;xf;xe])),'w');
        polarticks(8,g);
        hold on;
        f = polar(avgCoh.angle(iang),avgCoh.force.asig_coh.(bands{j}).mean(iang,i),'b-.o');
        f.MarkerFaceColor = 'b';
        hold on;
        e = polar(avgCoh.angle(iang),avgCoh.EMG.asig_coh.(bands{j}).mean(iang,i),'r-.o');
        e.MarkerFaceColor = 'r';        
        
        ff = fill(xf,yf,'b','edgecolor','none','facealpha',0.3);
        fe = fill(xe,ye,'r','edgecolor','none','facealpha',0.3);
        if i == 1
            title(bands_label{j})
        end
        if j == 1
            ylabel([avgCoh.musc{i}{1},'-',avgCoh.musc{i}{2}],'FontWeight','bold');
        end   
        set(gca,'FontSize',14);
    end
end

%% Area coh - Polar, one per musccomb
bands = {'alp','beta','gam'};
bands_label = {'Alpha','Beta','Gamma'};
for i = 1:length(avgCoh.musc) 
    figure('Name','Area Coherence Mean');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    iang = ismember(avgCoh.angle,avgCoh.angmusc{i});
    yf = []; ye = [];
    for j = 1:length(bands)
        subplot(1,length(bands),j)
        [xf,yf] = getBound(avgCoh.angle(iang),avgCoh.force.asig_coh.(bands{j}).mean(iang,i),avgCoh.force.asig_coh.(bands{j}).sem(iang,i),'cart');
        [xe,ye] = getBound(avgCoh.angle(iang),avgCoh.EMG.asig_coh.(bands{j}).mean(iang,i),avgCoh.EMG.asig_coh.(bands{j}).sem(iang,i),'cart');
        g = polar(0,max(abs([yf;ye;xf;xe]))+0.1*max(abs([yf;ye;xf;xe])),'w');
        polarticks(8,g);
        hold on;
        f = polar(avgCoh.angle(iang),avgCoh.force.asig_coh.(bands{j}).mean(iang,i),'b-.o');
        f.MarkerFaceColor = 'b';
        hold on;
        e = polar(avgCoh.angle(iang),avgCoh.EMG.asig_coh.(bands{j}).mean(iang,i),'r-.o');
        e.MarkerFaceColor = 'r';        
        
        ff = fill(xf,yf,'b','edgecolor','none','facealpha',0.3);
        fe = fill(xe,ye,'r','edgecolor','none','facealpha',0.3);
        if i == 1
            title(bands_label{j})
        end
        if j == 1
            ylabel([avgCoh.musc{i}{1},'-',avgCoh.musc{i}{2}],'FontWeight','bold');
        end   
        set(gca,'FontSize',14);
    end
end

%% Area coh - Plot
figure('Name','Significant Coherence Area');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
k = 0;
bands = {'alp','beta','gam'};
bands_label = {'Alpha','Beta','Gamma'};
for i = 1:length(avgCoh.musc)
    iang = ismember(avgCoh.angle,avgCoh.angmusc{i});
    yf = []; ye = [];
    for j = 1:length(bands)
        k = k+1;
        subplot(length(avgCoh.musc),length(bands),k)
        f = plot(rad2deg(avgCoh.angle(iang)),avgCoh.force.asig_coh.(bands{j}).mean(iang,i),'b-.o');
        f.MarkerFaceColor = 'b';
        hold on;
        errorbar(rad2deg(avgCoh.angle(iang)),avgCoh.force.asig_coh.(bands{j}).mean(iang,i),avgCoh.force.asig_coh.(bands{j}).sem(iang,i),'b.');
        e = plot(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.asig_coh.(bands{j}).mean(iang,i),'r-.o');
        e.MarkerFaceColor = 'r';
        errorbar(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.asig_coh.(bands{j}).mean(iang,i),avgCoh.EMG.asig_coh.(bands{j}).sem(iang,i),'r.');
        if i == 1
            title(bands_label{j})
        end
        if j == 1
            ylabel([avgCoh.musc{i}{1},'-',avgCoh.musc{i}{2}],'FontWeight','bold');
        end
        if i == length(avgCoh.musc)
            xlabel('Target [deg]');
        end
        set(gca,'XTick',rad2deg(avgCoh.angle(iang)));
        set(gca,'FontSize',14);
    end
end

%% Area z - Plot
figure('Name','Significant Z-score Area');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
k = 0;
bands = {'alp','beta','gam'};
bands_label = {'Alpha','Beta','Gamma'};
for i = 1:length(avgCoh.musc)
    iang = ismember(avgCoh.angle,avgCoh.angmusc{i});
    yf = []; ye = [];
    for j = 1:length(bands)
        k = k+1;
        subplot(length(avgCoh.musc),length(bands),k)
        f = plot(rad2deg(avgCoh.angle(iang)),avgCoh.force.asig_z.(bands{j}).mean(iang,i),'b-.o');
        f.MarkerFaceColor = 'b';
        hold on;
        errorbar(rad2deg(avgCoh.angle(iang)),avgCoh.force.asig_z.(bands{j}).mean(iang,i),avgCoh.force.asig_z.(bands{j}).sem(iang,i),'b.');
        e = plot(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.asig_z.(bands{j}).mean(iang,i),'r-.o');
        e.MarkerFaceColor = 'r';
        errorbar(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.asig_z.(bands{j}).mean(iang,i),avgCoh.EMG.asig_z.(bands{j}).sem(iang,i),'r.');
        if i == 1
            title(bands_label{j})
        end
        if j == 1
            ylabel([avgCoh.musc{i}{1},'-',avgCoh.musc{i}{2}],'FontWeight','bold');
        end
        if i == length(avgCoh.musc)
            xlabel('Target [deg]');
        end
        set(gca,'XTick',rad2deg(avgCoh.angle(iang)));
        set(gca,'FontSize',14);
    end
end

%% Norm area coh - Plot
figure('Name','Significant Norm Coherence Area');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
k = 0;
bands = {'alp','beta','gam'};
bands_label = {'Alpha','Beta','Gamma'};
for i = 1:length(avgCoh.musc)
    iang = ismember(avgCoh.angle,avgCoh.angmusc{i});
    for j = 1:length(bands)
        k = k+1;
        subplot(length(avgCoh.musc),length(bands),k)
        f = plot(rad2deg(avgCoh.angle(iang)),avgCoh.force.nasig_coh.(bands{j}).mean(iang,i),'b-.o');
        f.MarkerFaceColor = 'b';
        hold on;
        errorbar(rad2deg(avgCoh.angle(iang)),avgCoh.force.nasig_coh.(bands{j}).mean(iang,i),avgCoh.force.nasig_coh.(bands{j}).sem(iang,i),'b.');
        e = plot(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.nasig_coh.(bands{j}).mean(iang,i),'r-.o');
        e.MarkerFaceColor = 'r';
        errorbar(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.nasig_coh.(bands{j}).mean(iang,i),avgCoh.EMG.nasig_coh.(bands{j}).sem(iang,i),'r.');
        if i == 1
            title(bands_label{j})
        end
        if j == 1
            ylabel([avgCoh.musc{i}{1},'-',avgCoh.musc{i}{2}],'FontWeight','bold');
        end
        if i == length(avgCoh.musc)
            xlabel('Target [deg]');
        end
        set(gca,'XTick',rad2deg(avgCoh.angle(iang)));
        set(gca,'FontSize',13);      
    end
end

%% Norm area z - Plot
srt = [3,2,1];
sangpair = {avgCoh.angmusc{3},avgCoh.angmusc{2},avgCoh.angmusc{1}};
smuscpair = {avgCoh.muscComp{3},avgCoh.muscComp{2}(srt),avgCoh.muscComp{1}};

figure('Name','Significant Norm Z-score Area');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
k = 0;
bands = {'alp','beta','gam'};
bands_label = {'Alpha','Beta','Gamma'};
for j = 1:length(bands)
    for i = 1:length(avgCoh.musc)
        iang = ismember(avgCoh.angle,sangpair{i});        
        k = k+1;
        subplot(length(avgCoh.musc),length(bands),k)
        f = plot(rad2deg(avgCoh.angle(iang)),avgCoh.force.nasig_z.(bands{j}).mean(iang,srt(i)),'b-.o','MarkerFaceColor','b');
        hold on;
        errorbar(rad2deg(avgCoh.angle(iang)),avgCoh.force.nasig_z.(bands{j}).mean(iang,srt(i)),avgCoh.force.nasig_z.(bands{j}).sem(iang,srt(i)),'b.','CapSize',12);
        e = plot(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.nasig_z.(bands{j}).mean(iang,srt(i)),'r-.o','MarkerFaceColor','r');
        errorbar(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.nasig_z.(bands{j}).mean(iang,srt(i)),avgCoh.EMG.nasig_z.(bands{j}).sem(iang,srt(i)),'r.','CapSize',12);
        xlim(rad2deg([sangpair{i}(1)-0.2 sangpair{i}(3)+0.2]));
        if i == 1
            ylabel({['\bf',bands_label{j}];'\rmArea Z [-]'},'interpreter','tex');
        end
        if j == 1
            title([smuscpair{i}{1},'-',smuscpair{i}{3}],'FontWeight','bold');
        end
        if j == length(avgCoh.musc)
            xlabel('Target [deg]');
        end
        if j == 1 && i == 3
            legend([f,e],'FC','MC');
        end
        set(gca,'XTick',rad2deg(avgCoh.angle(iang)));
        set(gca,'FontSize',13);
    end
end

%% Norm area z - Plot only joint ang
srt = [3,2,1];
sangpair = {avgCoh.angmusc{3},avgCoh.angmusc{2},avgCoh.angmusc{1}};
smuscpair = {avgCoh.musc{3},avgCoh.musc{2}([2,1]),avgCoh.musc{1}};

figure('Name','Significant Norm Z-score Area');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
k = 0;
bands = {'alp','beta','gam'};
bands_label = {'Alpha','Beta','Gamma'};
for i = 1:length(avgCoh.musc)
    iang = ismember(avgCoh.angle,sangpair{i}(2));
    k = k+1;
    subplot(1,length(avgCoh.musc),k)
    for j = 1:length(bands)
        f = plot(j,avgCoh.force.nasig_z.(bands{j}).mean(iang,srt(i)),'bo','MarkerFaceColor','b');
        hold on;
        errorbar(j,avgCoh.force.nasig_z.(bands{j}).mean(iang,srt(i)),avgCoh.force.nasig_z.(bands{j}).sem(iang,srt(i)),'b.','CapSize',12);
        e = plot(j,avgCoh.EMG.nasig_z.(bands{j}).mean(iang,srt(i)),'ro','MarkerFaceColor','r');
        errorbar(j,avgCoh.EMG.nasig_z.(bands{j}).mean(iang,srt(i)),avgCoh.EMG.nasig_z.(bands{j}).sem(iang,srt(i)),'r.','CapSize',12);
    end
    xlim([0.5 3.5]); ylim([0 1.8])
    xlabel('Frequency Band');
    title([smuscpair{i}{1},'-',smuscpair{i}{2} ' (',num2str(rad2deg(sangpair{i}(2))),')'],'FontWeight','bold');
    if i == 1   
        ylabel('Area Z [-]');        
    end
    if i == 3
        legend([f,e],'FC','MC');
    end
    set(gca,'XTick',1:3,'XTickLabels',bands_label);
    set(gca,'FontSize',13);
end

%% Grand average z
srt = [3,2,1];
sangpair = {avgCoh.angmusc{3},avgCoh.angmusc{2},avgCoh.angmusc{1}};
smuscpair = {avgCoh.muscComp{3},avgCoh.muscComp{2}(srt),avgCoh.muscComp{1}};

h = 0;
figure('Name','Grand Average Z-score');
for k = 1:length(smuscpair)
    for i = 1:length(sangpair)
        
        iang = find([avgCoh.angle] == sangpair{i}(k));
        h = h+1;
        subplot(length(sangpair),length(sangpair),h);
        
        h1 = plot(avgCoh.force.my_fcoh,avgCoh.force.z(iang).mean(:,srt(i)),'b','Linewidth',2);
        hold on;
        h2 = plot(avgCoh.force.my_fcoh,avgCoh.EMG.z(iang).mean(:,srt(i)),'r','Linewidth',2);
        line(xlim,[0.62 0.62],'Color','k','LineStyle','-.');
        xlim([2 fc]); %ylim([0 coh_lim.(field)(Aparams.angCompUni==sangpair{i}(k),srt(i))]);
        if k == 1
            title({['\fontsize{15}',smuscpair{i}{1},'-',smuscpair{i}{3}];...
                ['Target: ',num2str(rad2deg(sangpair{i}(k))),' (',smuscpair{i}{k},')']},'interpreter','tex');
        else
            title(['Target: ',num2str(rad2deg(sangpair{i}(k))),' (',smuscpair{i}{k},')']);
        end
        if k == 3
            xlabel('Frequency [Hz]');
        end
        if i == 1
            ylabel('Z [-]');
        end
        if i == 3 && k == 1
            legend([h1,h2],'FC','EC')
        end
        set(gca,'FontSize',13);
    end
end

%% Z participant counts
sangpair = {avgCoh.angmusc{3},avgCoh.angmusc{2},avgCoh.angmusc{1}};
smuscpair = {avgCoh.musc{3},avgCoh.musc{2},avgCoh.musc{1}};
srt = [3,2,1];
fc = 80; h = 0;
figure('Name','Grand Average Z-score');
for k = 1:length(smuscpair)
    for i = 1:length(sangpair)
        
        iang = find([avgCoh.angle] == sangpair{i}(k));
        h = h+1;
        subplot(length(sangpair),length(sangpair),h);
        
        h1 = plot(avgCoh.force.my_fcoh,avgCoh.force.my_coh(iang).scount(:,srt(i)),'b','Linewidth',2);
        hold on;
        h2 = plot(avgCoh.force.my_fcoh,avgCoh.EMG.my_coh(iang).scount(:,srt(i)),'r','Linewidth',2);
        line(xlim,[1 1],'Color','k','LineStyle','-.');
        xlim([2 fc]); ylim([0 7.5]);
        title(['Target: ',num2str(rad2deg(sangpair{i}(k))),' deg (',avgCoh.muscComp{srt(i)}{k},')']);
        if k == 3
            xlabel('Frequency [Hz]');
        end
        if i == 1
            ylabel('Subjects');
        end
        if i == 3 && k == 1
            legend([h1,h2],'FC','EC')
        end
        set(gca,'FontSize',13);
    end
end