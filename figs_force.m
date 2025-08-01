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
%subplot(2,1,1)
for i = 1:length(Aparams.targetAnglesForce)
    plot(trial_avg_force(i).force.filt(:,1),trial_avg_force(i).force.filt(:,2),'b');
    hold on;
    h1 = plot(trial_avg_force(i).force.filt(1,1),trial_avg_force(i).force.filt(1,2),'go');
    h2 = plot(trial_avg_force(i).force.filt(end,1),trial_avg_force(i).force.filt(end,2),'mo');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]);ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    text(trial_avg_force(i).force.filt(end,1),trial_avg_force(i).force.filt(end,2),num2str(rad2deg(Aparams.targetAnglesForce(i))),'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12);
    legend([h1 h2],'Start','End');
end
title('Force-Control');
set(gca,'FontSize',12);

figure('Name','Force trajectory');
%subplot(2,1,2)
for i = 1:length(Aparams.targetAnglesEMG)
    plot(trial_avg_EMG(i).force.filt(:,1),trial_avg_EMG(i).force.filt(:,2),'r');
    hold on;
    h1 = plot(trial_avg_EMG(i).force.filt(1,1),trial_avg_EMG(i).force.filt(1,2),'go');
    h2 = plot(trial_avg_EMG(i).force.filt(end,1),trial_avg_EMG(i).force.filt(end,2),'mo');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]);ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    text(trial_avg_EMG(i).force.filt(end,1),trial_avg_EMG(i).force.filt(end,2),num2str(rad2deg(Aparams.targetAnglesEMG(i))),'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12);
    legend([h1 h2],'Start','End');
end
title('EMG-Control');
set(gca,'FontSize',12);

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
%title('Force-Control');
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

%% Force mean and CV Polar plot
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
title({'Force Magnitude'; 'Mean'})
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
title({'Force Magnitude'; 'CV'})
set(gca,'FontSize',18);

legend([h1,h2],'FC','MC','Location','bestoutside');

%% Force CV plot
figure('Name','Force CV');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for k = 1:length(Aparams.angCompUni)
    iangf = find([trial_avg_force.angle] == Aparams.angCompUni(k));
    iangE = find([trial_avg_EMG.angle] == Aparams.angCompUni(k));
    
    h1 = plot(rad2deg(trial_avg_force(iangf).angle),trial_avg_force(iangf).force.filtmag_CV,'bo','MarkerFaceColor','b');
    hold on;
    %errorbar(trial_avg_force(iangf).angle,trial_avg_force(iangf).force.filtmag_CV_mean,trial_avg_force(iangf).force.filtmag_CV_std,'b.');
    h2 = plot(rad2deg(trial_avg_EMG(iangE).angle),trial_avg_EMG(iangE).force.filtmag_CV,'ro','MarkerFaceColor','r');
    %errorbar(trial_avg_EMG(iangE).angle,trial_avg_EMG(iangE).force.filtmag_CV,trial_avg_EMG(iangE).force.filtmag_CV_std,'r.');
end
xticks(rad2deg(Aparams.angCompUni))
title('Force CV')
legend([h1,h2],'ForceCO','EMGCO');
set(gca,'FontSize',13);


