function[trial_data] = procEMG(trial_data,Aparams)
% Initialize variables
fchEMG = 5;
sampleRateEMG = 1024;
avgWindow = 200;
struct2vars(who,Aparams)

% Filtering parameters
wn = (2/sampleRateEMG)*fchEMG;
[b,a] = butter(2,wn,'high');

for i = 1:length(trial_data)
    EMGfilt = filtfilt(b,a,trial_data(i).EMG.raw);
    EMGrect = abs(EMGfilt);
    EMGavg = movingAvg(EMGrect,avgWindow);
    
    trial_data(i).EMG.filt = EMGfilt;
    trial_data(i).EMG.rect = EMGrect;
    trial_data(i).EMG.avg = EMGavg;
    
    trial_data(i).EMG.raw_fft = fft(trial_data(i).EMG.raw);
    trial_data(i).EMG.filt_fft = fft(EMGfilt);
    trial_data(i).EMG.rect_fft = fft(EMGrect);
end
end