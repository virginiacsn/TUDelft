
figure;
plot(trial_data_force(1).force.raw(:,2))
hold on;
plot(trial_data_force(1).force.filt(:,2))
%%
itrial = 5;
sampleRateEMG = 1024;
fchEMG = 30;
fclEMG = [];
window = 800;
overlap = 100;

EMGmean = [];
forcemean = [];
for i = 1:length(trial_data_avg_force)
    EMGmean(i,:) = mean(trial_data_avg_force(i).EMG.rect,1);
    forcemean(i,:) = mean(trial_data_avg_force(i).force.filt,1);
end
EMGscf = max(EMGmean,[],1)';

EMGmean = [];
forcemean = [];
for i = 1:length(trial_data_avg_EMG)
    EMGmean(i,:) = mean(trial_data_avg_EMG(i).EMG.rect,1);
end
EMGscE = max(EMGmean,[],1)';

EMGDataBuffer = cell2mat(EMGDataOut_EMGCO(2:end));%trial_data_EMG(itrial).EMG.raw;
avgrectEMG = runningAvg(EMGDataBuffer(:,1:end-1),window,overlap,fchEMG,fclEMG);
avgrectEMG = avgrectEMG./repmat(EMGScal,[1 size(avgrectEMG,1)])';

EMGDataBufferF = cell2mat(EMGDataOut_ForceCO(2:end));%trial_data_EMG(itrial).EMG.raw;
avgrectEMGF = runningAvg(EMGDataBufferF(:,1:end-1),window,overlap,fchEMG,fclEMG);
avgrectEMGF = avgrectEMGF./repmat(EMGScal,[1 size(avgrectEMGF,1)])';

% wnh = (2/sampleRateEMG)*fchEMG;
% wnl = (2/sampleRateEMG)*fclEMG;
% [b,a] = butter(2,wnh,'high');
% [d,c] = butter(2,wnl,'low');
% 
% filtEMGBuffer = filtfilt(b,a,EMGDataBuffer);
% filtEMGBuffer = filter(d,c,filtEMGBuffer,[],1);
% rectEMG = abs(filtEMGBuffer);
% 
% avgg = movingAvg(rectEMG,600);

% Rectified and smoothed EMG
figure;
set(gcf,'Name','Run avg EMG - EMGCO');
for j = 1:length(EMGparams.channelControl)
    %subplot(1,length(EMGparams.channelControl),j);
    plot(avgrectEMG(10:end,EMGparams.channelControl(j)));
    hold on;
    xlabel('Time [s]'); ylabel('EMG [-]');
end
legend(EMGparams.channelName(EMGparams.channelControl))

figure;
set(gcf,'Name','Run avg EMG - ForceCO');
for j = 1:length(EMGparams.channelControl)
    %subplot(1,length(EMGparams.channelControl),j);
    plot(avgrectEMGF(10:end,EMGparams.channelControl(j)));
    hold on;
    xlabel('Time [s]'); ylabel('EMG [-]');
end
legend(EMGparams.channelName(EMGparams.channelControl))

% figure;
% set(gcf,'Name','Rectified EMG');
% for j = 1:length(EMGparams.channelControl)
%     %subplot(1,length(EMGparams.channelControl),j);
%     plot(trial_data_EMG(itrial).EMG.avg(:,EMGparams.channelControl(j))./EMGparams.EMGScale(EMGparams.channelControl(j)));
%     hold on;
%     xlabel('Time [s]'); ylabel('EMG [-]');
% end

% figure;
% set(gcf,'Name','Rectified EMG');
% for j = 1:length(EMGparams.channelControl)
%     subplot(1,length(EMGparams.channelControl),j);
%     plot(rectEMG(:,EMGparams.channelControl(j)));
%     xlabel('Time [s]'); ylabel('EMG [-]');
% end
%% 
fchEMG = [10,20,30];
fclEMG = [30,40,60];
window = 800;
overlap = 100;

EMGDataBuffer = cell2mat(EMGDataOut_EMGCO(2:end));%trial_data_EMG(itrial).EMG.raw;

for i = 1:length(fchEMG)
avgrectEMG = runningAvg(EMGDataBuffer(:,1:end-1),window,overlap,fchEMG(i),fclEMG(i));

figure('Name',['fch = ',num2str(fchEMG(i))]);
for j = 1:length(EMGparams.channelControl)
    %subplot(1,length(EMGparams.channelControl),j);
    plot(avgrectEMG(10:end,EMGparams.channelControl(j)));
    hold on;
    xlabel('Time [s]'); ylabel('EMG [-]');
end
end
%%
figure
plot(emg_save(1,:)*EMGparams.EMGScale(1))
hold on
plot(emg_save(2,:)*EMGparams.EMGScale(2),'r')
title('Raw')


%% 
figure
plot(emg_save(1,:))
hold on
plot(emg_save(2,:),'r')
title('Scaled')

%% 
EMGraw = cell2mat(EMGDataOut_EMGCO(2:end))';
figure
plot(EMGraw([3 6],1:2000));

