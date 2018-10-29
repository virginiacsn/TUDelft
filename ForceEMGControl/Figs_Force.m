%% FORCE FIGURES (low-pass filtered)
%% Force in time. Fig per task. Subplot per target angle.
figure('Name','ForceCO; Force in time');
for i = 1:length(Aparams.targetAnglesForce)
    if rem(length(Aparams.targetAnglesForce),2) == 0
        subplot(2,length(Aparams.targetAnglesForce)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesForce),i);
    end
    for j = 1:2
        plot(trial_avg_force(i).ts,trial_avg_force(i).force.filt(:,j));
        hold on;
    end
    ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Time [s]'); ylabel('Force [N]');
    title(['Target: ',num2str(rad2deg(trial_avg_force(i).angle)),' deg']);
    legend('Fx','Fy')
end

figure('Name','EMGCO; Force in time');
for i = 1:length(Aparams.targetAnglesEMG)
    if rem(length(Aparams.targetAnglesEMG),2) == 0
        subplot(2,length(Aparams.targetAnglesEMG)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesEMG),i);
    end
    for j = 1:2
        plot(trial_avg_EMG(i).ts,trial_avg_EMG(i).force.filt(:,j));
        hold on;
    end
    ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Time [s]'); ylabel('Force [N]');
    title(['Target: ',num2str(rad2deg(trial_avg_EMG(i).angle)),' deg']);
    legend('Fx','Fy')
end

%% Force trajectory. Fig per task. Subplot per target angle.
figure('Name','ForceCO; Force trajectory');
for i = 1:length(Aparams.targetAnglesForce)
    if rem(length(Aparams.targetAnglesForce),2) == 0
        subplot(2,length(Aparams.targetAnglesForce)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesForce),i);
    end
    plot(trial_avg_force(i).force.filt(:,1),trial_avg_force(i).force.filt(:,2));
    hold on;
    h1 = plot(trial_avg_force(i).force.filt(1,1),trial_avg_force(i).force.filt(1,2),'go');
    h2 = plot(trial_avg_force(i).force.filt(end,1),trial_avg_force(i).force.filt(end,2),'ro');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]);ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    title(['Target: ',num2str(rad2deg(trial_avg_force(i).angle)),' deg']);
    legend([h1,h2],'Start','End');
end

figure('Name','EMGCO; Force trajectory');
for i = 1:length(Aparams.targetAnglesEMG)
    if rem(length(Aparams.targetAnglesEMG),2) == 0
        subplot(2,length(Aparams.targetAnglesEMG)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesEMG),i);
    end
    plot(trial_avg_EMG(i).force.filt(:,1),trial_avg_EMG(i).force.filt(:,2));
    hold on;
    h1 = plot(trial_avg_EMG(i).force.filt(1,1),trial_avg_EMG(i).force.filt(1,2),'go');
    h2 = plot(trial_avg_EMG(i).force.filt(end,1),trial_avg_EMG(i).force.filt(end,2),'ro');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]);ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    title(['Target: ',num2str(rad2deg(trial_avg_EMG(i).angle)),' deg']);
    legend([h1,h2],'Start','End');
end

%% Force trajectory. Subplot per task. Target angles in same subplot.
figure('Name','Force trajectory');
subplot(1,2,1)
for i = 1:length(Aparams.targetAnglesForce)
    plot(trial_avg_force(i).force.filt(:,1),trial_avg_force(i).force.filt(:,2));
    hold on;
    h1 = plot(trial_avg_force(i).force.filt(1,1),trial_avg_force(i).force.filt(1,2),'go');
    h2 = plot(trial_avg_force(i).force.filt(end,1),trial_avg_force(i).force.filt(end,2),'ro');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]);ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    text(trial_avg_force(i).force.filt(end,1),trial_avg_force(i).force.filt(end,2),num2str(rad2deg(Aparams.targetAnglesForce(i))),'VerticalAlignment','top','HorizontalAlignment','left');
    legend([h1 h2],'Start','End');
end
title('ForceCO');

subplot(1,2,2)
for i = 1:length(Aparams.targetAnglesEMG)
    plot(trial_avg_EMG(i).force.filt(:,1),trial_avg_EMG(i).force.filt(:,2));
    hold on;
    h1 = plot(trial_avg_EMG(i).force.filt(1,1),trial_avg_EMG(i).force.filt(1,2),'go');
    h2 = plot(trial_avg_EMG(i).force.filt(end,1),trial_avg_EMG(i).force.filt(end,2),'ro');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]);ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    text(trial_avg_EMG(i).force.filt(end,1),trial_avg_EMG(i).force.filt(end,2),num2str(rad2deg(Aparams.targetAnglesEMG(i))),'VerticalAlignment','top','HorizontalAlignment','left');
    legend([h1 h2],'Start','End');
end
title('EMGCO');

%% Force Mean Bar
figure('Name','Force Mean');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(Aparams.angComp)
    iangf = find([trial_avg_force.angle] == Aparams.angComp{i}(1));
    iangE = find([trial_avg_EMG.angle] == Aparams.angComp{i}(2));
    if rem(length(Aparams.angComp),2) == 0
        subplot(length(Aparams.angComp)/2,2,i);
    else
        subplot(1,length(Aparams.angComp),i);
    end
    bar(1,trial_avg_force(iangf).force.filtmag_mean,'barwidth',0.9,'FaceColor','b');
    hold on;
    bar(2,trial_avg_EMG(iangE).force.filtmag_mean,'barwidth',0.9,'FaceColor','r');
    errorbar([1 2],[trial_avg_force(iangf).force.filtmag_mean,...
        trial_avg_EMG(iangE).force.filtmag_mean],...
        [trial_avg_force(iangf).force.filtmag_std,trial_avg_EMG(iangE).force.filtmag_std],'k.');
    ylim([0 taskparams.targetForce*Flim+5])
    ylabel('Mean [-]');
    title(['Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),' deg (',Aparams.muscComp{i},')']);
    set(gca,'XTick',1:2,'XTickLabel',{'ForceCO','EMGCO'});
end

%% Force CV Bar
figure('Name','Force CV');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(Aparams.angComp)
    iangf = find([trial_avg_force.angle] == Aparams.angComp{i}(1));
    iangE = find([trial_avg_EMG.angle] == Aparams.angComp{i}(2));
    if rem(length(Aparams.angComp),2) == 0
        subplot(length(Aparams.angComp)/2,2,i);
    else
        subplot(1,length(Aparams.angComp),i);
    end
    bar(1,trial_avg_force(iangf).force.filtmag_pstd/trial_avg_force(iangf).force.filtmag_mean,'barwidth',0.9,'FaceColor','b');
    hold on;
    bar(2,trial_avg_EMG(iangE).force.filtmag_pstd/trial_avg_EMG(iangE).force.filtmag_mean,'barwidth',0.9,'FaceColor','r');
    ylabel('CV [-]');
    title(['Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),' deg (',Aparams.muscComp{i},')']);
    set(gca,'XTick',1:2,'XTickLabel',{'ForceCO','EMGCO'});
end

%% Force CV Polar plot
figure('Name','Force CV');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
fCV = zeros(length(length(Aparams.angCompUni)),3);
ECV = zeros(length(length(Aparams.angCompUni)),3);
for k = 1:length(Aparams.angCompUni)
    iangf = find([trial_avg_force.angle] == Aparams.angCompUni(k));
    iangE = find([trial_avg_EMG.angle] == Aparams.angCompUni(k));
    
    fCV(k) = 100*trial_avg_force(iangf).force.filtmag_pstd/trial_avg_force(iangf).force.filtmag_mean;
    ECV(k) = 100*trial_avg_EMG(iangE).force.filtmag_pstd/trial_avg_EMG(iangE).force.filtmag_mean;
end
polarscatter(Aparams.angCompUni,fCV,60,'filled')
hold on
polarscatter(Aparams.angCompUni,ECV,60,'filled','r')
thetaticks(rad2deg(Aparams.angCompUni)); %thetaticklabels(Aparams.muscComp);
set(gca,'FontSize',12);
legend('ForceCO','EMGCO','Location','bestoutside');
