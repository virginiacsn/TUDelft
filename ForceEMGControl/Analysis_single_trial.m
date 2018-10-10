%% SINGLE-TRIAL
%% Single trial from trial_data
itrial = 20;
trial_data = trial_data;
figure('Name',['Target: ', num2str(rad2deg(trial_data(itrial).angle)), ' deg']);
for i = 1:length(EMGparams.channelName)-1
    subplot(1,length(EMGparams.channelName)-1,i);
    plot(trial_data(itrial).ts,trial_data(itrial).EMG.rect(:,i));
    hold on;
    plot(trial_data(itrial).ts,trial_data(itrial).EMG.avg(:,i),'r');
    xlim([0 trial_data(itrial).ts(end)]);
    %ylim([0 max(trial_data(itrial).EMG.rect(:,i))]);
    xlabel('Time [s]'); ylabel('EMG [-]');
    %title(['Musc: ',Aparams.channelNameEMG{i}]);
end

%% Compare tasks
load([filepath,'Parameters/',[date,'_s',subject,'_params_','001','.mat']])

forceAngles = [trial_data_force.angle];
EMGAngles = [trial_data_EMG.angle];

compAngles = Aparams.targetAnglesForce(ismember(Aparams.targetAnglesForce,Aparams.targetAnglesEMG));

for i = 1:length(compAngles)
    trialsAngleForce{i} = find(forceAngles == compAngles(i));
    trialsAngleEMG{i}  = find(EMGAngles == compAngles(i));
end

%% EMG
% Rectified and smoothed EMG
iAngle = 5;
for h = 8:10
    itrialForce = trialsAngleForce{iAngle}(h);
    itrialEMG = trialsAngleEMG{iAngle}(h);
    
    figure('Name',[' Target: ',num2str(rad2deg(trial_data_force(itrialForce).angle)),' deg']);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for j = 1:length(Aparams.channelNameEMG)-1
        subplot(2,length(Aparams.channelNameEMG)-1,j);
        plot(trial_data_force(itrialForce).ts,trial_data_force(itrialForce).EMG.rect(:,j));
        hold on;
        plot(trial_data_force(itrialForce).ts,trial_data_force(itrialForce).EMG.avg(:,j));
        xlabel('Time [s]'); ylabel('EMG [-]');
        xlim([0 trial_data_force(itrialForce).ts(end)]);
        ylim([0 max(trial_data_force(itrialForce).EMG.rect(:,j))]);
        title(['Musc: ',Aparams.channelNameEMG{j}]);
        
    end
    for j = 1:length(Aparams.channelNameEMG)-1
        subplot(2,length(Aparams.channelNameEMG)-1,j+length(Aparams.channelNameEMG)-1);
        plot(trial_data_EMG(itrialEMG).ts,trial_data_EMG(itrialEMG).EMG.rect(:,j));
        hold on;
        plot(trial_data_EMG(itrialEMG).ts,trial_data_EMG(itrialEMG).EMG.avg(:,j));
        xlabel('Time [s]'); ylabel('EMG [-]');
        xlim([0 trial_data_EMG(itrialEMG).ts(end)]);
        ylim([0 max(trial_data_EMG(itrialEMG).EMG.rect(:,j))]);
        title(['Musc: ',Aparams.channelNameEMG{j}]);
    end
end
% %% EMG
% % Rectified and smoothed EMG
% figure;
% set(gcf,'Name','Rectified EMG');
% for j = 1:nmusc
%     subplot(1,nmusc,j);
%     plot(trial_data(itrial).ts,trial_data(itrial).EMG.rect(:,j));
%     hold on;
%     plot(trial_data(itrial).ts,trial_data(itrial).EMG.avg(:,j));
%     xlabel('Time [s]'); ylabel('EMG [-]');
%     ylim([0 max(trial_data(itrial).EMG.rect(:))+50]);
%     title(['Musc: ',EMGparams.channelName{j},'; Target: ',num2str(rad2deg(trial_data(itrial).angle)),' deg']);
% end
% 
% % EMG trajectory
% figure;
% set(gcf,'Name','Trajectory Rectified EMG');
% plot(trial_data(itrial).EMG.avg(:,musccont(1)),trial_data(itrial).EMG.avg(:,musccont(2)));
% hold on;
% h1 = plot(trial_data(itrial).EMG.avg(1,musccont(1)),trial_data(itrial).EMG.avg(1,musccont(2)),'go');
% h2 = plot(trial_data(itrial).EMG.avg(end,musccont(1)),trial_data(itrial).EMG.avg(end,musccont(2)),'ro');
% grid on;
% axis square;
% xlabel(EMGparams.channelName{musccont(1)}); ylabel(EMGparams.channelName{musccont(2)});
% title(['Musc: ',EMGparams.channelName{musccont(1)},',',EMGparams.channelName{musccont(2)},'; Target: ',num2str(rad2deg(trial_data(itrial).angle)),' deg']);
% legend([h1,h2],'Start','End');
% 
% % FFT of rectified EMG
% figure;
% set(gcf,'Name','FFT Rectified EMG');
% for j = 1:nmusc
%     plot(trial_data(itrial).fv,abs(trial_data(itrial).EMG.rect_fft(:,j))/length(trial_data(itrial).EMG.rect_fft(:,j)));
%     hold on;
% end
% xlabel('Frequency [Hz]'); ylabel('FFT EMG [-]');
% xlim([trial_data(itrial).fv(2) trial_data(itrial).fv(round(end/2))]);
% title(['Target: ',num2str(rad2deg(trial_data(itrial).angle)),' deg']);
% legend(EMGparams.channelName{1:end-1});
% 
% %% Force
% % LPF force in time
% figure;
% set(gcf,'Name','Force in time');
% for j = 1:2
%     plot(trial_data(itrial).ts,trial_data(itrial).force.filt(:,j));
%     hold on;
% end
% %ylim(taskparams.targetForce*[-1.5 1.5]);
% xlabel('Time [s]'); ylabel('Force [N]');
% title(['Trial: ',num2str(itrial),'; Target: ',num2str(rad2deg(trial_data(itrial).angle)),' deg']);
% legend('Fx','Fy')
% 
% % LPF force trajectory
% figure;
% set(gcf,'Name','Force trajectory');
% plot(trial_data(itrial).force.filt(:,1),trial_data(itrial).force.filt(:,2));
% hold on;
% h1 = plot(trial_data(itrial).force.filt(1,1),trial_data(itrial).force.filt(1,2),'go');
% h2 = plot(trial_data(itrial).force.filt(end,1),trial_data(itrial).force.filt(end,2),'ro');
% xlim(taskparams.targetForce*[-2 2]);ylim(taskparams.targetForce*[-2 2]);
% grid on;
% axis square;
% xlabel('Fx [N]'); ylabel('Fy [N]');
% title(['Trial: ',num2str(itrial),'; Target: ',num2str(rad2deg(trial_data(itrial).angle)),' deg']);
% legend([h1,h2],'Start','End');
