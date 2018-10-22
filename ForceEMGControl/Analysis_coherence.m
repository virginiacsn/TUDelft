%% COHERENCE ANALYSIS
Aparams.cohparams.data = 'app';
Aparams.cohparams.tseg = 1;
Aparams.cohparams.nseg = 10;
Aparams.cohparams.my_nseg = 10;
Aparams.cohparams.window = @(N) hanning(N);
Aparams.cohparams.CLoverlap = 1;
fields_coh = {'filt','rect'};
fc = 100;

if strcmp(Aparams.cohparams.data,'avg')
    trial_data_coh_force = cohStruct(trial_data_avg_force,EMGparams.channelName,fields_coh,Aparams.cohparams);
    trial_data_coh_EMG = cohStruct(trial_data_avg_EMG,EMGparams.channelName,fields_coh,Aparams.cohparams);
else
    trial_data_coh_force = cohStruct(trial_data_app_force,EMGparams.channelName,fields_coh,Aparams.cohparams);
    trial_data_coh_EMG = cohStruct(trial_data_app_EMG,EMGparams.channelName,fields_coh,Aparams.cohparams);
end

% Number of muscles
nmusc = length(EMGparams.channelName)-1;
% Number of muscle pair combinations
nmusccomb = (length(EMGparams.channelName)-1)*(length(EMGparams.channelName)-2)/2; % n*(n-1)/2

%% Field and channels to plot
field = 'rect';

% Find muscle in all muscle pair combinations
plotBB = sum(reshape(contains([trial_data_coh_force(1).(field).muscles{:}],'BB'),[2,nmusccomb]));
plotTLH = sum(reshape(contains([trial_data_coh_force(1).(field).muscles{:}],'TLH'),[2,nmusccomb]));
plotDA = sum(reshape(contains([trial_data_coh_force(1).(field).muscles{:}],'DA'),[2,nmusccomb]));
plotDP = sum(reshape(contains([trial_data_coh_force(1).(field).muscles{:}],'DP'),[2,nmusccomb]));

% Find indexes of muscle pair combinations to plot. Combinations of control
% muscles.
plotmusc = find((plotBB&plotDA)|(plotDP&plotBB)|(plotDA&plotTLH)|(plotDA&plotDP)|(plotBB&plotTLH)|(plotTLH&plotDP));

% Coherence (my_coh) limits for each field in stuct. coh_lim(i,j) will be max coh 
% level between force-control and EMG-control where rows are target angles 
% to compare and cols muscle combinations (all)
for h = 1:length(fields_coh)
    for i = 1:length(Aparams.angComp)
        for j = 1:nmusccomb
            coh_lim.(fields_coh{h})(i,j) = max([max(trial_data_coh_force(Aparams.targetAnglesForce == Aparams.angComp{i}(1)).(fields_coh{h}).my_coh(3:end-3,j)),...
                max(trial_data_coh_EMG(Aparams.targetAnglesEMG == Aparams.angComp{i}(2)).(fields_coh{h}).my_coh(3:end-3,j))]);
        end
    end
end

%% Coherence 
%% Fig per muscle pair in plotmusc. Subplot per target from angComp. Line per task.
% Compare coherence of each muscle combination (plotmusc) between tasks and targets from angComp
for j = 1:length(plotmusc)
    figure('Name',['My Coherence: Musc: ',trial_data_coh_force(1).(field).muscles{plotmusc(j)}{1},...
        ',',trial_data_coh_force(1).(field).muscles{plotmusc(j)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(length(Aparams.angComp),1,i);
        end
        h1 = plot(trial_data_coh_force(iangf).(field).my_fcoh(:,plotmusc(j)),trial_data_coh_force(iangf).(field).my_coh(:,plotmusc(j)));
        hold on;
        line(xlim,trial_data_coh_force(iangf).(field).my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_data_coh_EMG(iangE).(field).my_fcoh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).my_coh(:,plotmusc(j)),'r');
        hold on;
        line(xlim,trial_data_coh_EMG(iangE).(field).my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','-.');
        
        xlim([2 fc]); ylim([0 coh_lim.(field)(i,plotmusc(j))]);
        if i == length(Aparams.angComp)
            xlabel('Frequency [Hz]');
        end
        ylabel('Coh [-]');
        title(['Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),' deg (',Aparams.muscComp{i},')']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end

%% Fig per muscle pair (all). Subplot per target from angComp. Line per task.
% Compare coherence of each muscle combination (all) between tasks and targets from angComp
for j = 1:nmusccomb
    figure('Name','My Coherence');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(length(Aparams.angComp),1,i);
        end
        h1 = plot(trial_data_coh_force(iangf).(field).my_fcoh(:,j),trial_data_coh_force(iangf).(field).my_coh(:,j));
        hold on;
        line(xlim,trial_data_coh_force(iangf).(field).my_CL(j)*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_data_coh_EMG(iangE).(field).my_fcoh(:,j),trial_data_coh_EMG(iangE).(field).my_coh(:,j),'r');
        hold on;
        line(xlim,trial_data_coh_EMG(iangE).(field).my_CL(j)*[1 1],'Color','k','LineStyle','-.');
                
        xlim([trial_data_coh_force(iangf).(field).my_fcoh(2,j) fc]);
        %ylim([0 1]);
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_data_coh_force(iangf).(field).muscles{j}{1},...
            ',',trial_data_coh_force(iangf).(field).muscles{j}{2},...
            '; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end

%% Fig per muscle pair (plotmusc). Subplot per task (col) per target from angComp (row). Line per processing method.
% Compare processing method (filt/rect) for each muscle combination between
% task and target from angComp
for j = 1:length(plotmusc)
    figure('Name',['My Coherence - Comparison Proc; Musc: ',trial_data_coh_force(1).filt.muscles{plotmusc(j)}{1},...
        ',',trial_data_coh_force(1).filt.muscles{plotmusc(j)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    h = 1;
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        subplot(length(Aparams.angComp),2,h);
        h1 = plot(trial_data_coh_force(iangf).filt.my_fcoh(:,plotmusc(j)),trial_data_coh_force(iangf).filt.my_coh(:,plotmusc(j)),'b');
        hold on;
        line(xlim,trial_data_coh_force(iangf).filt.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','--');
        h2 = plot(trial_data_coh_force(iangf).rect.my_fcoh(:,plotmusc(j)),trial_data_coh_force(iangf).rect.my_coh(:,plotmusc(j)),'c');
        line(xlim,trial_data_coh_force(iangf).rect.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','-.');
        xlim([2 fc]); ylim([0 max([coh_lim.filt(i,plotmusc(j)),coh_lim.rect(i,plotmusc(j))])]);        
        if i == length(Aparams.angComp)
            xlabel('Frequency [Hz]'); 
        end
        ylabel('Coh [-]');
        title(['ForceCO; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'Filt','Rect')
        
        subplot(length(Aparams.angComp),2,h+1);
        h1 = plot(trial_data_coh_EMG(iangE).filt.my_fcoh(:,plotmusc(j)),trial_data_coh_EMG(iangE).filt.my_coh(:,plotmusc(j)),'r');
        hold on;
        line(xlim,trial_data_coh_EMG(iangE).filt.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','--');
        h2 = plot(trial_data_coh_EMG(iangE).rect.my_fcoh(:,plotmusc(j)),trial_data_coh_EMG(iangE).rect.my_coh(:,plotmusc(j)),'m');
        line(xlim,trial_data_coh_EMG(iangE).rect.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','-.');
        xlim([2 fc]); ylim([0 max([coh_lim.filt(i,plotmusc(j)),coh_lim.rect(i,plotmusc(j))])]);        
        if i == length(Aparams.angComp)
            xlabel('Frequency [Hz]'); 
        end
        ylabel('Coh [-]');
        title(['EMGCO; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'Filt','Rect')
        
        h = h+2;
    end
end

%% Significant coherence
%% Fig of mean significant coherence per muscle pair (plotmusc). Subplot per target from angComp. Bar per frequency band and task.
% Compare mean significant coherence (SEM for CI) for each muscle
% combination (plotmusc) between task and frequency bands 
bandCol = {'c','g','r'};
freqBand = {'Alpha','Beta','Gamma'};

for j = 1:length(plotmusc)
    figure('Name',['Significant coherence: Musc: ',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{1},...
        ',',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));

        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(1,length(Aparams.angComp),i);
        end
        bar([1 2 3],[trial_data_coh_force(iangf).(field).my_mean_sig_coh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).my_mean_sig_coh(:,plotmusc(j))],'barwidth',0.9);
        hold on;
        errorbar([1-0.15 1+0.15; 2-0.15 2+0.15; 3-0.15 3+0.15],[trial_data_coh_force(iangf).(field).my_mean_sig_coh(:,plotmusc(j)),...
            trial_data_coh_EMG(iangE).(field).my_mean_sig_coh(:,plotmusc(j))],...
            [trial_data_coh_force(iangf).(field).my_SEM_sig_coh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).my_SEM_sig_coh(:,plotmusc(j))],'k.');%,'facecolor',bandCol{i},'barwidth',0.9);
        ylim([0 1]);
        xlabel('Frequency [Hz]'); ylabel('Sig Coh [-]');
        title(['Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),' deg (',Aparams.muscComp{i},')']);
        set(gca,'XTick',1:3,'XTickLabel',freqBand);
        legend('ForceCO','EMGCO');
    end
end

%% Fig of percent of significant coherence per muscle pair (plotmusc). Subplot per target from angComp. Bar per frequency band and task.
% Compare percent of significant coherence of each frequency band for each 
% muscle pair (plotmusc) between task and frequency band
freqBand = {'Alpha','Beta','Gamma'};

for j = 1:length(plotmusc)
    figure('Name',['Significant coherence: Musc: ',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{1},...
        ',',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        %mean_coh_alp(i,:) = [trial_data_coh_force(iangf).(field).sig_coh(1,plotmusc(j))
        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(1,length(Aparams.angComp),i);
        end
        bar([1 2 3],[trial_data_coh_force(iangf).(field).my_perc_sig_coh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).my_perc_sig_coh(:,plotmusc(j))],'barwidth',0.9);
        hold on;
        errorbar([1-0.15 1+0.15; 2-0.15 2+0.15; 3-0.15 3+0.15],[trial_data_coh_force(iangf).(field).my_perc_sig_coh(:,plotmusc(j)),...
            trial_data_coh_EMG(iangE).(field).my_perc_sig_coh(:,plotmusc(j))],...
            [trial_data_coh_force(iangf).(field).my_perc_CI_sig_coh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).my_perc_CI_sig_coh(:,plotmusc(j))],'k.');%,'facecolor',bandCol{i},'barwidth',0.9);
        %ylim([0 1]);
        xlabel('Frequency [Hz]'); ylabel('Perc Sig Coh [-]');
        title(['Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),' deg (',Aparams.muscComp{i},')']);
        set(gca,'XTick',1:3,'XTickLabel',freqBand);
        legend('ForceCO','EMGCO');
    end
end

%% Correlation
%% Fig per muscle pair in plotmusc. Subplot per target from angComp. Line per task.
% Compare cross-correlation of each muscle combination (plotmusc) between 
% tasks and targets from angComp
for j = 1:length(plotmusc)
    figure('Name',['Cross-Correlation: Musc: ',trial_data_coh_force(1).(field).muscles{plotmusc(j)}{1},...
        ',',trial_data_coh_force(1).(field).muscles{plotmusc(j)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        if rem(length(Aparams.angComp),2) == 0
            subplot(2,length(Aparams.angComp)/2,i);
        else
            subplot(1,length(Aparams.angComp),i);
        end
        
        plot(trial_data_coh_force(iangf).(field).lags(:,plotmusc(j))*trial_data_avg_force(iangf).ts(2),trial_data_coh_force(iangf).(field).corr(:,plotmusc(j)));
        hold on;        
        plot(trial_data_coh_EMG(iangE).(field).lags(:,plotmusc(j))*trial_data_avg_EMG(iangE).ts(2),trial_data_coh_EMG(iangE).(field).corr(:,plotmusc(j)),'r');
        
        maxf = max(abs(trial_data_coh_force(iangf).(field).corr(:,plotmusc(j))));
        maxE = max(abs(trial_data_coh_EMG(iangE).(field).corr(:,plotmusc(j))));

        xlim([-0.1 0.1]);%ylim([-max(maxf,maxE) max(maxf,maxE)]);
        xlabel('Time [s]');
        if i == 1
            ylabel('Corr [-]');
        end
        title(['Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),' deg (',Aparams.muscComp{i},')']);
        legend('ForceCO','EMGCO')
    end
end

%% Fig per muscle pair (plotmusc) and task. Subplot per target from angComp.
% Compare correlation between targets from angComp for each task and muscle
% combination (plotmusc)
for i = 1:length(plotmusc)
    figure('Name',['ForceCO; Musc: ',trial_data_coh_force(j).(field).muscles{plotmusc(i)}{1},...
            ',',trial_data_coh_force(j).(field).muscles{plotmusc(i)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for j = 1:length(Aparams.targetAnglesForce)
        if rem(length(Aparams.targetAnglesForce),2) == 0
            subplot(2,length(Aparams.targetAnglesForce)/2,j);
        else
            subplot(1,length(Aparams.targetAnglesForce),j);
        end
        plot(trial_data_coh_force(j).(field).lags(:,plotmusc(i))*trial_data_avg_force(j).ts(2),trial_data_coh_force(j).(field).corr(:,plotmusc(i)));      
        maxf = max(abs(trial_data_coh_force(j).(field).corr(:,plotmusc(i))));
        xlim([-0.1 0.1]); ylim([-maxf maxf]);
        xlabel('Time lag [s]'); ylabel('Corr [-]');
        title(['Target: ',num2str(rad2deg(Aparams.targetAnglesForce(j))),' deg']);
    end
end

for i = 1:length(plotmusc)
    figure('Name',['EMGCO; Musc: ',trial_data_coh_EMG(j).(field).muscles{plotmusc(i)}{1},...
            ',',trial_data_coh_EMG(j).(field).muscles{plotmusc(i)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for j = 1:length(Aparams.targetAnglesEMG)
        if rem(length(Aparams.targetAnglesEMG),2) == 0
            subplot(2,length(Aparams.targetAnglesEMG)/2,j);
        else
            subplot(1,length(Aparams.targetAnglesEMG),j);
        end
        plot(trial_data_coh_EMG(j).(field).lags(:,plotmusc(i))*trial_data_avg_EMG(j).ts(2),trial_data_coh_EMG(j).(field).corr(:,plotmusc(i)),'r');
        maxE = max(abs(trial_data_coh_EMG(j).(field).corr(:,plotmusc(i))));
        xlim([-0.1 0.1]); ylim([-maxE maxE]);
        xlabel('Time lag [s]'); ylabel('Corr [-]');
        title(['Target: ',num2str(rad2deg(Aparams.targetAnglesEMG(j))),' deg']);
    end
end

%% MATLAB coherence
%% Fig per muscle pair (plotmusc). Subplot per target from angComp. Line per task.
% Compare coherence of each muscle combination (plotmusc) between tasks and 
% targets from angComp
for j = 1:length(plotmusc)
    figure('Name','MATLAB Coherence');
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(length(Aparams.angComp),1,i);
        end
        h1 = plot(trial_data_coh_force(iangf).(field).fcoh(:,plotmusc(j)),trial_data_coh_force(iangf).(field).coh(:,plotmusc(j)));
        hold on;
        line(xlim,trial_data_coh_force(iangf).(field).CL(plotmusc(j))*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_data_coh_EMG(iangE).(field).fcoh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).coh(:,plotmusc(j)),'r');
        hold on;
        line(xlim,trial_data_coh_EMG(iangE).(field).CL(plotmusc(j))*[1 1],'Color','k','LineStyle','-.');
        
        xlim([trial_data_coh_force(iangf).(field).fcoh(2,plotmusc(j)) fc]);
        ylim([0 0.8]);
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{1},...
            ',',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{2},...
            '; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end

%% Fig per muscle pair (all). Subplot per target from angComp. Line per task.
% Compare coherence of each muscle combination (all) between tasks and 
% targets from angComp
for j = 1:nmusccomb
    figure('Name','Coherence');
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));

        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(1,length(Aparams.angComp),i);
        end
        h1 = plot(trial_data_coh_force(iangf).(field).fcoh(:,j),trial_data_coh_force(iangf).(field).coh(:,j));
        hold on;
        line(xlim,trial_data_coh_force(iangf).(field).CL(j)*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_data_coh_EMG(iangE).(field).fcoh(:,j),trial_data_coh_EMG(iangE).(field).coh(:,j),'r');
        hold on;
        line(xlim,trial_data_coh_EMG(iangE).(field).CL(j)*[1 1],'Color','k','LineStyle','-.');
        
        xlim([trial_data_coh_force(iangf).(field).fcoh(2,j) fc])
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_data_coh_force(iangf).(field).muscles{j}{1},...
            ',',trial_data_coh_force(iangf).(field).muscles{j}{2},...
            '; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end

%% Scrap
% %% Save trial data structs
% if ~exist([filepath,'TrialData/'],'dir')
%     mkdir([filepath,'TrialData/']);
% end
% if exist('trial_data_EMG_calib','var')
%     save([filepath,'TrialData/trial_data_calib'],'trial_data_EMG_calib','trial_data_avg_calib');
% end
% save([filepath,'TrialData/trial_data'],'trial_data_force','trial_data_EMG','Aparams');
% save([filepath,'TrialData/trial_data_avg'],'trial_data_avg_force','trial_data_avg_EMG');
% save([filepath,'TrialData/trial_data_app'],'trial_data_app_force','trial_data_app_EMG');

% for i = 1:6
%     for kk = 1:4
%         subplot(2,3,i)
%         h(kk)=bar(kk,mpt(kk,i),'facecolor',cols{kk},'barwidth',0.9); hold on;
%         errorbar(kk,mpt(kk,i),spt(kk,i),'k.')
%         title(strcat(header{i+4},' (\mu,\sigma)'));ylabel('Score');xlim([0 5]);
%         set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
%     end
% end