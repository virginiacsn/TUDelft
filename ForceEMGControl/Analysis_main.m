%% Data Analysis
date =      '20180720';
task =      'EMGCO';
code =      '002';
EMG =       1;
filenameforce =  [date,'_',task,'_Force_',code,'.mat'];
filenameEMG = [date,'_',task,'_EMG_',code,'.mat'];
filepath =  ['D:\Student_experiments\Virginia\Data\' date '\'];

load([filepath,filenameforce]);
if EMG
    load([filepath,filenameEMG]);
end 

forceEMGData = {forceDataOut,EMGDataOut};
params.downsample = 2;
trial_data = trialCO(forceEMGData,params);
%trial_data = removeFailTrials(trial_data);

%% Plot
itrial = 2;
wn = (2/2048)*5;
[b,a] = butter(2,wn,'low');
forceDataFilt = filter(b,a,trial_data(itrial).force);

figure;
plot(forceDataFilt(:,1));
hold on
plot(forceDataFilt(:,2));

figure;
plot(trial_data(itrial).trigger);
hold on
plot(trial_data(itrial).EMGtrigger);