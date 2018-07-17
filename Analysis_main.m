%% Data Analysis
date =      '20180716';
task =      'CO';
code =      '006';
EMG =       0;
filenameforce =  [date,'_',task,'_',code,'.mat'];
filenameEMG = [date,'_',task,'_EMG_',code,'.mat'];
filepath =  [pwd '\Data\' date '\'];

load([filepath,filenameforce]);
if EMG
    load([filepath,filenameEMG]);
end

forceEMGData = {forceDataOut,EMGDataOut};
params.downsample = 2;
trial_data = trialCO(forceEMGData,params);
%trial_data = removeFailTrials(trial_data);

%% Plot
itrial = 1;
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