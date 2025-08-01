%% Final figs for all subjects/population
%% Force magnitude and CV mean - Polar
figure('Name','Force CV');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
h1 = polarscatter(angComp+pi/4,[avgforce_force.mag_mean_mean],80,'filled','b');
hold on;
h2 = polarscatter(angComp+pi/4,[avgEMG_force.mag_mean_mean],80,'filled','r');
errorpolar(angComp+pi/4,[avgforce_force.mag_mean_mean],[avgforce_force.mag_mean_sem],'b')
errorpolar(angComp+pi/4,[avgEMG_force.mag_mean_mean],[avgEMG_force.mag_mean_sem],'r')
thetaticks(degAngComp+45); thetaticklabels(degAngComp);
title({'Force Magnitude'; 'Mean [N]'})
set(gca,'FontSize',20);

subplot(1,2,2)
h1 = polarscatter(avgCV.angle+pi/4,avgCV.force.mean,80,'filled','b');
hold on;
h2 = polarscatter(avgCV.angle+pi/4,avgCV.EMG.mean,80,'filled','r');
errorpolar(avgCV.angle+pi/4,avgCV.force.mean,avgCV.force.sem,'b');
errorpolar(avgCV.angle+pi/4,avgCV.EMG.mean,avgCV.EMG.sem,'r');
thetaticks(degAngComp+45); thetaticklabels(degAngComp);
title({'Force Magnitude'; 'CV [%]'})
set(gca,'FontSize',20);
legend([h1,h2],'FC','MC','Location','bestoutside');

%% Scaled EMG mean - Plot. Subplot per control muscle.
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

% %% Normalized area significant coherence - Plot
% figure('Name','Significant Norm Coherence Area');
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% k = 0;
% bands = {'alp','beta','gam'};
% bands_label = {'Alpha','Beta','Gamma'};
% for i = 1:length(avgCoh.musc)
%     iang = ismember(avgCoh.angle,avgCoh.angmusc{i});
%     for j = 1:length(bands)
%         k = k+1;
%         subplot(length(avgCoh.musc),length(bands),k)
%         f = plot(rad2deg(avgCoh.angle(iang)),avgCoh.force.nasig_coh.(bands{j}).mean(iang,i),'b-.o');
%         f.MarkerFaceColor = 'b';
%         hold on;
%         errorbar(rad2deg(avgCoh.anglclosee(iang)),avgCoh.force.nasig_coh.(bands{j}).mean(iang,i),avgCoh.force.nasig_coh.(bands{j}).sem(iang,i),'b.');
%         e = plot(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.nasig_coh.(bands{j}).mean(iang,i),'r-.o');
%         e.MarkerFaceColor = 'r';
%         errorbar(rad2deg(avgCoh.angle(iang)),avgCoh.EMG.nasig_coh.(bands{j}).mean(iang,i),avgCoh.EMG.nasig_coh.(bands{j}).sem(iang,i),'r.');
%         plotang = avgCoh.angle(iang);
%         xlim(rad2deg([plotang(1)-0.2 plotang(end)+0.2]));
%         if i == 1
%             title(bands_label{j})
%         end
%         if j == 1
%             ylabel([avgCoh.musc{i}{1},'-',avgCoh.musc{i}{2}],'FontWeight','bold');
%         end
%         if i == length(avgCoh.musc)
%             xlabel('Target [deg]');
%         end
%         set(gca,'XTick',rad2deg(avgCoh.angle(iang)));
%         set(gca,'FontSize',13);      
%     end
% end

%% Normalized area significant Z - Plot
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
        xlim(rad2deg([sangpair{i}(1)-0.2 sangpair{i}(3)+0.2])); ylim([0, max(ylim)]);
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

%% Normalized area significant Z - Plot. Only joint-muscle targets.
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
    xlim([0.5 3.5]); ylim([0, max(ylim)]);
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

%% Grand average Z
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
        
        plot(avgCoh.force.my_fcoh,avgCoh.force.z(iang).mean(:,srt(i)),'b','Linewidth',2);
        hold on;
        plot(avgCoh.force.my_fcoh,avgCoh.EMG.z(iang).mean(:,srt(i)),'r','Linewidth',2);
        xlim([2 fc]);
        yl = ylim;
        area([8 12],[yl(2)-0.02 yl(2)-0.02],'FaceColor',[0.8,0.8,0.8],'FaceAlpha',0.8,'EdgeColor','none');
        area([15 30],[yl(2)-0.02 yl(2)-0.02],'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.8,'EdgeColor','none');
        h1 = plot(avgCoh.force.my_fcoh,avgCoh.force.z(iang).mean(:,srt(i)),'b','Linewidth',2);
        h2 = plot(avgCoh.force.my_fcoh,avgCoh.EMG.z(iang).mean(:,srt(i)),'r','Linewidth',2);
        line(xlim,[0.65 0.65],'Color','k','LineStyle','-.');

        
%         plot([15 15],ylim,'Color',[0.3,0.3,0.3],'Linewidth',1);
%         plot([30 30],ylim,'Color',[0.3,0.3,0.3],'Linewidth',1);
        ylim([0, max(ylim)]);
        if k == 1
            title({['\fontsize{15}',smuscpair{i}{1},'-',smuscpair{i}{3}];...
                ['Target: ',num2str(rad2deg(sangpair{i}(k))),' deg (',smuscpair{i}{k},')']},'interpreter','tex');
        else
            title(['Target: ',num2str(rad2deg(sangpair{i}(k))),' deg (',smuscpair{i}{k},')']);
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