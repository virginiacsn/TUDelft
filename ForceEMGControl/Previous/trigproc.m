EMGDataOut = cell2mat(EMGDataOut_EMGCO(2:end));%trial_data_EMG(itrial).EMG.raw;
EMGdata = EMGDataOut(:,1:end-1);
EMG_trigger = EMGDataOut(:,end);
EMG_trigger(EMG_trigger ~= 6) = 0;
EMG_trigger(EMG_trigger > 0) = 1;

EMGdiff = diff(EMG_trigger);
EMGpulse = sum(EMGdiff(EMGdiff>0));
EMGpulseind = find(EMGdiff>0);
EMGindext = EMGpulseind(diff(EMGpulseind)>1024)+1;
EMGsamppulse = sum(EMG_trigger>0)/EMGpulse;

fDataOut = (forceDataOut_EMGCO(2:end,:));%trial_data_EMG(itrial).EMG.raw;
f_trigger = cell2mat(fDataOut(:,end));
f_trigger(f_trigger < 3) = 0;
f_trigger(f_trigger > 3) = 1;

fdiff = diff(f_trigger);
fpulse = sum(fdiff(fdiff>0));
fpulseind = find(fdiff>0);
fsamppulse = sum(f_trigger>0)/fpulse;

figure
plot(f_trigger)

%%

figure;
plot(trial_data_force(10).trigger.EMG)
hold on;
plot(trial_data_force(10).trigger.force)

%% Calib testing
code = '001';

filenameforce =  [date,'_s',subject,'_',task,'_Force_',code,'.mat'];
filenameEMG =    [date,'_s',subject,'_',task,'_EMG_',code,'.mat'];

load([fileparams.filepath,fileparams.filenameforce]);
load([fileparams.filepath,fileparams.filenameEMG]);

forceEMGData = {forceDataOut_ForceCO,EMGDataOut_ForceCO};
PreAparams.downsample = 2;
PreAparams.target_angles = taskparams.targetAnglesForce;
PreAparams.avgWindow = 200;
PreAparams.fclF = forceparams.fclF;
PreAparams.fchEMG = 20;

trial_data = trialCO(forceEMGData,PreAparams);

trial_data = removeFailTrials(trial_data);

trial_data = procEMG(trial_data,PreAparams);
trial_data = procForce(trial_data,PreAparams);

epoch = {'ihold','iend'};
fields = {'EMG.raw','force.filtmag','force.rawmag','force.filt','force.raw'};
trial_data_avg = trialAngleAvg(trial_data, epoch, fields);
trial_data_avg = procEMG(trial_data_avg,PreAparams);

EMGmean = zeros(length(trial_data_avg),length(EMGparams.channelSubset)-1);
forcemean = zeros(length(trial_data_avg),1);
for i = 1:length(trial_data_avg)
    EMGmean(i,:) = mean(trial_data_avg(i).EMG.rect,1);
    forcemean(i) = trial_data_avg(i).force.filtmag_mean;
end

EMGScale = max(EMGmean,[],1)';
targetForce = round(mean(forcemean))*0.8;
