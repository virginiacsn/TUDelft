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
    EMGfilt = filtfilt(b,a,trial_data(i).EMG);
    EMGrect = abs(EMGfilt);
    EMGavg = movingavg(EMGrect,avgWindow);
    
    trial_data(i).EMGfilt = EMGfilt;
    trial_data(i).EMGrect = EMGrect;
    trial_data(i).EMGavg = EMGavg;
end
end