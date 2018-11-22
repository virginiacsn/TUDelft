%% Final figs for single-subject
%% Force trajectory. Target angles in same subplot for both tasks.
figure('Name','Force trajectory');
%subplot(2,1,1)
for i = 1:length(Aparams.targetAnglesForce)
    h1 = plot(trial_avg_force(i).force.filt(:,1),trial_avg_force(i).force.filt(:,2),'b');
    hold on;
    plot(trial_avg_force(i).force.filt(1,1),trial_avg_force(i).force.filt(1,2),'go');
    plot(trial_avg_force(i).force.filt(end,1),trial_avg_force(i).force.filt(end,2),'mo');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]+1);ylim(taskparams.targetForce*[-Flim Flim]+1);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    text(trial_avg_force(i).force.filt(end,1),trial_avg_force(i).force.filt(end,2),num2str(rad2deg(Aparams.targetAnglesForce(i))),'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12);
    %legend([h1 h2],'Start','End');
end
set(gca,'FontSize',12);

%figure('Name','Force trajectory');
%subplot(2,1,2)
for i = 1:length(Aparams.targetAnglesEMG)
    h2 = plot(trial_avg_EMG(i).force.filt(:,1),trial_avg_EMG(i).force.filt(:,2),'r');
    hold on;
    plot(trial_avg_EMG(i).force.filt(1,1),trial_avg_EMG(i).force.filt(1,2),'go');
    plot(trial_avg_EMG(i).force.filt(end,1),trial_avg_EMG(i).force.filt(end,2),'mo');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]+1);ylim(taskparams.targetForce*[-Flim Flim]+1);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    text(trial_avg_EMG(i).force.filt(end,1),trial_avg_EMG(i).force.filt(end,2),num2str(rad2deg(Aparams.targetAnglesEMG(i))),'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12);
    %legend([h1 h2],'Start','End');
end
legend([h1 h2],'FC','MC');
set(gca,'FontSize',13);

%% Force magnitude and CV - Polar
figure('Name','Force Mean and CV');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

subplot(1,2,1)
fm = zeros(length(length(Aparams.angCompUni)),3);
Em = zeros(length(length(Aparams.angCompUni)),3);
fstd = zeros(length(length(Aparams.angCompUni)),3);
Estd = zeros(length(length(Aparams.angCompUni)),3);
for k = 1:length(Aparams.angCompUni)
    iangf = find([trial_avg_force.angle] == Aparams.angCompUni(k));
    iangE = find([trial_avg_EMG.angle] == Aparams.angCompUni(k));
    
    fm(k) = trial_avg_force(iangf).force.filtmag_mean;
    Em(k) = trial_avg_EMG(iangE).force.filtmag_mean;
    fstd(k) = trial_avg_force(iangf).force.filtmag_std;
    Estd(k) = trial_avg_EMG(iangE).force.filtmag_std;
end
h1 = polarscatter(Aparams.angCompUni,fm,60,'filled','b');
hold on
errorpolar(Aparams.angCompUni,fm,fstd,'b')
h2 = polarscatter(Aparams.angCompUni,Em,60,'filled','r');
errorpolar(Aparams.angCompUni,Em,Estd,'r')
thetaticks([rad2deg(Aparams.angCompUni)]); thetaticklabels(rad2deg(Aparams.angCompUni));
title({'Force Magnitude'; 'Mean [N]'})
set(gca,'FontSize',18);

subplot(1,2,2)
fm = zeros(length(length(Aparams.angCompUni)),3);
Em = zeros(length(length(Aparams.angCompUni)),3);
for k = 1:length(Aparams.angCompUni)
    iangf = find([trial_avg_force.angle] == Aparams.angCompUni(k));
    iangE = find([trial_avg_EMG.angle] == Aparams.angCompUni(k));
    
    fm(k) = trial_avg_force(iangf).force.filtmag_CV;
    Em(k) = trial_avg_EMG(iangE).force.filtmag_CV;
end
h1 = polarscatter(Aparams.angCompUni,fm,60,'filled','b');
hold on
h2 = polarscatter(Aparams.angCompUni,Em,60,'filled','r');
thetaticks([rad2deg(Aparams.angCompUni)]); thetaticklabels(rad2deg(Aparams.angCompUni));
title({'Force Magnitude'; 'CV [%]'})
set(gca,'FontSize',18);

legend([h1,h2],'FC','MC','Location','bestoutside');

%% EMG in time. Fig per task. Subplot per muscle (control) (col) and target (row). 
% Compare control muscles and targets for each task.
h = 0;
figure('Name','ForceCO');
for j = 1:length(Aparams.targetAnglesForce)
    for i = 1:length(Aparams.chanControl)
        h = h+1;
        subplot(length(Aparams.targetAnglesForce),length(Aparams.chanControl),h);
        plot(trial_avg_force(j).ts,trial_avg_force(j).EMG.rectScale(:,Aparams.chanControl(i)));
        hold on;
        plot(trial_avg_force(j).ts,trial_avg_force(j).EMG.avgScale(:,Aparams.chanControl(i)),'LineWidth',2);
        ylim([0 2]);%ylim([0 EMG_lim(j,Aparams.chanControl(i))]); %ylim([0 max(trial_data_avg_force(j).EMG.rect(:))+10]);%
        xlim([0 trial_avg_force(j).ts(end)]);
        if j == length(Aparams.targetAnglesForce)
            xlabel('Time [s]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.targetAnglesForce(j)))],'FontWeight','bold');
        end
        if j == 1
            title([Aparams.chanControlName{i}]);
        end
        set(gca,'FontSize',13);
    end
end

h = 0;
figure('Name','EMGCO');
for j = 1:length(Aparams.targetAnglesEMG)
    for i = 1:length(Aparams.chanControl)
        h = h+1;
        subplot(length(Aparams.targetAnglesEMG),length(Aparams.chanControl),h);
        
        plot(trial_avg_EMG(j).ts,trial_avg_EMG(j).EMG.rectScale(:,Aparams.chanControl(i)));
        hold on;
        plot(trial_avg_EMG(j).ts,trial_avg_EMG(j).EMG.avgScale(:,Aparams.chanControl(i)),'LineWidth',2);
        ylim([0 2]);%ylim([0 EMG_lim(j,Aparams.chanControl(i))]);%ylim([0 max(trial_data_avg_EMG(j).EMG.rect(:))+50]);
        xlim([0 trial_avg_EMG(j).ts(end)]);
        if j == length(Aparams.targetAnglesEMG)
            xlabel('Time [s]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.targetAnglesEMG(j)))],'FontWeight','bold');
        end
        if j == 1
            title([Aparams.chanControlName{i}]);
        end
        set(gca,'FontSize',13);
    end
end

%% Scaled EMG mean - Plot. Subplot per control muscle.
figure('Name','Scaled EMG mean')
for i = 1:length(Aparams.chanControl)
    subplot(length(Aparams.chanControl),1,i)
    for j = 1:length(Aparams.angComp)
        iangf = find([trial_avg_force.angle] == Aparams.angComp{j}(1));
        
        h1 = plot(Aparams.angComp{j}(1),mean(trial_avg_force(iangf).EMG.rectScale(:,Aparams.chanControl(i))),'bo','MarkerFaceColor','b');
        hold on
        errorbar(Aparams.angComp{j}(1),mean(trial_avg_force(iangf).EMG.rectScale(:,Aparams.chanControl(i))),...
            std(trial_avg_force(iangf).EMG.rectScale(:,Aparams.chanControl(i))),'b.','Capsize',12)
        iangE = find([trial_avg_EMG.angle] == Aparams.angComp{j}(2));
        h2 = plot(Aparams.angComp{j}(2),mean(trial_avg_EMG(iangE).EMG.rectScale(:,Aparams.chanControl(i))),'ro','MarkerFaceColor','r');
        errorbar(Aparams.angComp{j}(2),mean(trial_avg_EMG(iangE).EMG.rectScale(:,Aparams.chanControl(i))),...
            std(trial_avg_EMG(iangE).EMG.rectScale(:,Aparams.chanControl(i))),'r.','Capsize',12)
    end
    ylim([0 1.5]);xlim([Aparams.angCompUni(1)-0.2 Aparams.angCompUni(end)+0.2])
    line(xlim,[1 1],'LineStyle',':','Color','k','Linewidth',1.7);
    line(xlim,[0.7 0.7],'LineStyle',':','Color','k','Linewidth',1.7);
    ylabel('EMG [-]');
    xticklabels(rad2deg(Aparams.angCompUni));
    set(gca,'XTick',[Aparams.angCompUni])
    title(Aparams.chanControlName(i));
    if i == 1
        legend([h1,h2],'FC','MC')
    elseif i == length(Aparams.chanControl)
        xlabel('Target [deg]');
    end
    set(gca,'FontSize',13);
end

%% Scaled EMG mean - Polar. Subplot per control muscle.
figure('Name','Scaled mean EMG');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(Aparams.chanControl)
    subplot(2,2,i)
    fm = zeros(length(length(Aparams.angCompUni)),3);
    Em = zeros(length(length(Aparams.angCompUni)),3);
    fstd = zeros(length(length(Aparams.angCompUni)),3);
    Estd = zeros(length(length(Aparams.angCompUni)),3);
    for k = 1:length(Aparams.angCompUni)
        iangf = find([trial_avg_force.angle] == Aparams.angCompUni(k));
        iangE = find([trial_avg_EMG.angle] == Aparams.angCompUni(k));
        
        fm(k) = mean(trial_avg_force(iangf).EMG.rectScale(:,Aparams.chanControl(i)));
        Em(k) = mean(trial_avg_EMG(iangE).EMG.rectScale(:,Aparams.chanControl(i)));
        fstd(k) = std(trial_avg_force(iangf).EMG.rectScale(:,Aparams.chanControl(i)));
        Estd(k) = std(trial_avg_EMG(iangE).EMG.rectScale(:,Aparams.chanControl(i)));
    end
    h1 = polarscatter(Aparams.angCompUni,fm,60,'filled','b');
    hold on
    errorpolar(Aparams.angCompUni,fm,fstd,'b')
    h2 = polarscatter(Aparams.angCompUni,Em,60,'filled','r');
    errorpolar(Aparams.angCompUni,Em,Estd,'r')
    thetaticks([rad2deg(Aparams.angCompUni)]); thetaticklabels(rad2deg(Aparams.angCompUni));
    rlim([0 1.2])
    title(Aparams.chanControlName(i))
    set(gca,'FontSize',18);
end

legend([h1,h2],'FC','MC','Location','bestoutside');

%% EMG FFT.
h = 0;
figure('Name','FFT EMG');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for j = 1:length(Aparams.angComp)
    for i = 1:length(Aparams.chanControl)
        h = h+1;
        subplot(length(Aparams.targetAnglesEMG),length(Aparams.chanControl),h);
        
        iangf = find([trial_avg_force.angle] == Aparams.angComp{j}(1));
        iangE = find([trial_avg_EMG.angle] == Aparams.angComp{j}(2));

        h1 = plot(trial_avg_force(iangf).fv,abs(trial_avg_force(iangf).EMG.fftrect(:,Aparams.chanControl(i))),'b');
        hold on;
        h2 = plot(trial_avg_EMG(iangE).fv,abs(trial_avg_EMG(iangE).EMG.fftrect(:,Aparams.chanControl(i))),'r');
        xlim([2 fc]); %ylim([0 max(max(fft_lim.rect(:,Aparams.chanControl(i))))]);
        if j == length(Aparams.angComp)
            xlabel('Frequency [Hz]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.angComp{j}(1)))],'FontWeight','bold');
        end
        if j == 1
            title([Aparams.chanControlName{i}]);
            if i == length(Aparams.chanControl)
                legend([h1,h2],{'FC','MC'})
            end
        end
        set(gca,'Fontsize',13);
    end
end

%% Coherence. Subplot per muscle pair (plotmuscpair) (col) and per target from angCompPair (row). 
figure('Name','Coherence');
h = 0;
sangpair = {Aparams.angCompPair{2},Aparams.angCompPair{3},Aparams.angCompPair{1}};
smuscpair = {Aparams.muscCompPair{2},Aparams.muscCompPair{3},Aparams.muscCompPair{1}};

for k = 1:length(Aparams.angCompPair)
    for i = 1:length(Aparams.angCompPair{k})
        
        plotmuscpair = find(sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],smuscpair{i}(1)),...
            [2,nmusccomb]))&sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],smuscpair{i}(3)),[2,nmusccomb])));
        
        iangf = find([trial_coh_force.angle] == sangpair{i}(k));
        iangE = find([trial_coh_EMG.angle] == sangpair{i}(k));
        h = h+1;
        subplot(length(Aparams.angCompPair),length(Aparams.angCompPair),h);
        
        h1 = plot(trial_coh_force(iangf).(field).my_fcoh(:,plotmuscpair),trial_coh_force(iangf).(field).my_coh(:,plotmuscpair),'b','Linewidth',2);
        hold on;
        line(xlim,trial_coh_force(iangf).(field).my_CL(plotmuscpair)*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_coh_EMG(iangE).(field).my_fcoh(:,plotmuscpair),trial_coh_EMG(iangE).(field).my_coh(:,plotmuscpair),'r','Linewidth',2);
        hold on;
        line(xlim,trial_coh_EMG(iangE).(field).my_CL(plotmuscpair)*[1 1],'Color','k','LineStyle','-.');
        
        xlim([2 fc]); ylim([0 coh_lim.(field)(Aparams.angCompUni==sangpair{i}(k),plotmuscpair)]);
        if k == 1
            title({['\fontsize{15}',smuscpair{i}{1},'-',smuscpair{i}{3}];...
                ['Target: ',num2str(rad2deg(sangpair{i}(k))),' (',smuscpair{i}{k},')']},'interpreter','tex');            
        else
            title(['Target: ',num2str(rad2deg(sangpair{i}(k))),' (',smuscpair{i}{k},')']);
        end
        if  k == length(Aparams.angCompPair{i})
            xlabel('Frequency [Hz]');
        end
        if i == 1
            ylabel('Coh [-]');
        end
        if i == length(Aparams.angCompPair{i}) && k == 1
            legend([h1,h2],'FC','EC')
        end
        set(gca,'FontSize',13);
    end
end

%% Normalized area of significant coherence for muscle pair (plotmuscpair) - Plot
sangpair = {Aparams.angCompPair{2},Aparams.angCompPair{3},Aparams.angCompPair{1}};
smuscpair = {Aparams.muscCompPair{2},Aparams.muscCompPair{3},Aparams.muscCompPair{1}};

figure('Name','Normalized area of significant coherence');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
h = 0;
for j = 1:length(freqBand)
    for i = 1:length(Aparams.muscCompPair)
        plotmuscpair = find(sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],smuscpair{i}(1)),...
            [2,nmusccomb]))&sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],smuscpair{i}(3)),[2,nmusccomb])));        
        fsig = [];
        Esig = [];        
        for k = 1:length(Aparams.angCompPair{i})
            iangf = find([trial_coh_force.angle] == sangpair{i}(k));
            iangE = find([trial_coh_EMG.angle] == sangpair{i}(k));
            
            fsig(k,:) = trial_coh_force(iangf).(field).nasig_coh(1:3,plotmuscpair)';
            Esig(k,:) = trial_coh_EMG(iangE).(field).nasig_coh(1:3,plotmuscpair)';
        end        
        h = h+1;
        subplot(length(Aparams.angCompPair),length(freqBand),h)
        f = plot(rad2deg(sangpair{i}),fsig(:,j),'b-.o','MarkerFaceColor','b');
        hold on;
        e = plot(rad2deg(sangpair{i}),Esig(:,j),'r-.o','MarkerFaceColor','r');
        if j == 1
            title([smuscpair{i}{1},'-',smuscpair{i}{3}],'FontWeight','bold');          
        elseif j == length(Aparams.muscCompPair)
            xlabel('Target [deg]');
        end
        if i == 1
            ylabel({['\bf',freqBand{j}];'\rmArea Coh [-]'},'interpreter','tex');
        end
        if j == 1 && i == 3
            legend([f,e],'FC','MC');
        end
        xlim(rad2deg([sangpair{i}(1)-0.2 sangpair{i}(3)+0.2])); ylim([0, max(ylim)]);
        set(gca,'XTick',rad2deg(sangpair{i}));
        set(gca,'FontSize',13);
    end
end

% %% No overlap. Coherence. Subplot per muscle pair (plotmuscpair) (col) and per target from angCompPair (row). 
% figure('Name','Coherence');
% h = 0;
% sangpair = {Aparams.angCompPair{2},Aparams.angCompPair{3},Aparams.angCompPair{1}};
% smuscpair = {Aparams.muscCompPair{2},Aparams.muscCompPair{3},Aparams.muscCompPair{1}};
% 
% for k = 1:length(Aparams.angCompPair)
%     for i = 1:length(Aparams.angCompPair{k})
%         
%         plotmuscpair = find(sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],smuscpair{i}(1)),...
%             [2,nmusccomb]))&sum(reshape(contains([trial_coh_force(1).(field).muscles{:}],smuscpair{i}(3)),[2,nmusccomb])));
%         
%         iangf = find([trial_coh_force.angle] == sangpair{i}(k));
%         iangE = find([trial_coh_EMG.angle] == sangpair{i}(k));
%         h = h+1;
%         subplot(length(Aparams.angCompPair),length(Aparams.angCompPair),h);
%         
%         h1 = plot(trial_coh_force_nov(iangf).(field).my_fcoh(:,plotmuscpair),trial_coh_force_nov(iangf).(field).my_coh(:,plotmuscpair),'b','Linewidth',2);
%         hold on;
%         line(xlim,trial_coh_force_nov(iangf).(field).my_CL(plotmuscpair)*[1 1],'Color','k','LineStyle','--');
%         
%         h2 = plot(trial_coh_EMG_nov(iangE).(field).my_fcoh(:,plotmuscpair),trial_coh_EMG_nov(iangE).(field).my_coh(:,plotmuscpair),'r','Linewidth',2);
%         hold on;
%         line(xlim,trial_coh_EMG_nov(iangE).(field).my_CL(plotmuscpair)*[1 1],'Color','k','LineStyle','-.');
%         
%         xlim([2 fc]); ylim([0 coh_lim.(field)(Aparams.angCompUni==sangpair{i}(k),plotmuscpair)]);
%         if k == 1
%             title({['\fontsize{15}',smuscpair{i}{1},'-',smuscpair{i}{3}];...
%                 ['Target: ',num2str(rad2deg(sangpair{i}(k))),' (',smuscpair{i}{k},')']},'interpreter','tex');            
%         else
%             title(['Target: ',num2str(rad2deg(sangpair{i}(k))),' (',smuscpair{i}{k},')']);
%         end
%         if  k == length(Aparams.angCompPair{i})
%             xlabel('Frequency [Hz]');
%         end
%         if i == 1
%             ylabel('Coh [-]');
%         end
%         if i == length(Aparams.angCompPair{i}) && k == 1
%             legend([h1,h2],'FC','EC')
%         end
%         set(gca,'FontSize',13);
%     end
% end
