itrial = 8;
sampleRateEMG = 1024;
fchEMG = 10;
fclEMG = 30;
window = 600;
overlap = 200;

EMGDataBuffer = trial_data(itrial).EMG.raw;

wnh = (2/sampleRateEMG)*fchEMG;
wnl = (2/sampleRateEMG)*fclEMG;
[b,a] = butter(2,wnh,'high');
[d,c] = butter(2,wnl,'low');

filtEMGBuffer = filtfilt(b,a,EMGDataBuffer);
filtEMGBuffer = filter(d,c,filtEMGBuffer,[],1);
rectEMG = abs(filtEMGBuffer);
avgrectEMG = runningAvg(rectEMG,window,overlap);

% Rectified and smoothed EMG
figure;
set(gcf,'Name','Rectified EMG');
for j = 1:nmusc
    subplot(1,nmusc,j);
    plot(avgrectEMG(:,j));
    xlabel('Time [s]'); ylabel('EMG [-]');
end

%% 
figure
plot(emg_save(1,:)*EMGparams.EMGScale(1))
hold on
plot(emg_save(2,:)*EMGparams.EMGScale(2),'r')

%% 
figure
plot(emg_save(1,:))
hold on
plot(emg_save(2,:),'r')

