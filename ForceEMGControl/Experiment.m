%% Experiment protocol
close all;
%clear all;

addpath(genpath('Tools'));

% Parameters
fileparams = struct(...
    'saveforce',    1,...
    'saveEMG',      1,...
    'date',         '20180925',...
    'subject',      '01');

if ~exist(['D:\Student_experiments\Virginia\Data\',fileparams.date],'dir') && (fileparams.saveEMG || fileparams.saveforce)
    mkdir(['D:\Student_experiments\Virginia\Data\',fileparams.date])
end
if ~exist(['D:\Student_experiments\Virginia\Data\',fileparams.date,'\s',fileparams.subject],'dir') && (fileparams.saveEMG || fileparams.saveforce)
    mkdir(['D:\Student_experiments\Virginia\Data\',fileparams.date,'\s',fileparams.subject])
end
fileparams.filepath =       ['D:\Student_experiments\Virginia\Data\',fileparams.date,'\s',fileparams.subject,'\'];

if ~exist([fileparams.filepath,'Parameters/'],'dir')
    mkdir([fileparams.filepath,'Parameters/'])
end

taskparams = struct(...
    'numTargetsForce',  7,...
    'numTargetsEMG',    3,...
    'targetForce',      10,... % [N]
    'targetForceCal',   30,... % [N]
    'targetEMG',        0.2,...
    'targetEMGCal',     40,...
    'targetTolForce',   0.1,...
    'targetTolEMG',     0.2,...
    'cursorTol',        2,...
    'movemtime',        5,... % sec
    'holdtime',         6,... % sec
    'timeout',          1,... % sec
    'relaxtime',        2,... % sec
    'setFig',           0);

if taskparams.numTargetsForce == 4
    taskparams.targetAnglesForce = [pi/4:pi/2:7*pi/4];
elseif taskparams.numTargetsForce == 7
    taskparams.targetAnglesForce = [pi/4:pi/4:7*pi/4];
end

forceparams = struct(...
    'fclF',             2,... % [Hz]
    'scanRate',         2048,... % [scans/sec]
    'availSamples',     200,... % [samples]
    'availSamplesEMG',  200,... % [samples]
    'bufferWin',        500,... % [samples]
    'iterUpdatePlot',   2,...
    'rotAngPlot',       -45); % [deg]

EMGparams = struct(...
    'plotEMG',          0,...
    'channelSubset',    [1 2 3 4 17],...
    'channelName',      {{'BB','TLH','DA','DP','Trigger'}},...
    'channelAngle',     [5*pi/4,pi/4,3*pi/4,7*pi/4],...%[0 0 5*pi/4 pi/4 pi/4 pi/4 7*pi/4],...
    'sampleRateEMG',    1024,... % [samples/sec]
    'fchEMG',           30,... % [Hz]
    'fclEMG',           60,... % [Hz]
    'smoothWin',        800,... % [samples]
    'iterUpdatePlotEMG',1);

%% Calibration with force-control task
fileparams.code = 'calib';
fileparams.task = 'ForceCO';

fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

ForceControl_Cal(fileparams,taskparams,forceparams,EMGparams);

%% Calibration with EMG-control task
fileparams.code = 'calib';
fileparams.task = 'EMGCO';

fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

EMGControl_Cal(fileparams,taskparams,forceparams,EMGparams);

%% Calibration with MVC task
EMGparams.EMGScaleMVC = MVCtest(EMGparams);

%% Pre-analysis
load([fileparams.filepath,fileparams.filenameforce]);
load([fileparams.filepath,fileparams.filenameEMG]);

if strcmp(fileparams.task,'ForceCO')
    forceEMGData = {forceDataOut_ForceCO,EMGDataOut_ForceCO};
    
    PreAparams.target_angles = taskparams.targetAnglesForce;
elseif strcmp(fileparams.task,'EMGCO')
    forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};
    
    PreAparams.target_angles = taskparams.targetAnglesEMG;
end

PreAparams.downsample = 2;
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

if strcmp(fileparams.task,'ForceCO')
    EMGparams.EMGScaleMVF = max(EMGmean,[],1)';
elseif strcmp(fileparams.task,'EMGCO')
    EMGparams.EMGScaleCC = max(EMGmean,[],1)';
end

taskparams.targetForce = round(mean(forcemean));

%% Force-control task
fileparams.code = '001';
fileparams.task = 'ForceCO';

fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

if strcmp(fileparams.task,'ForceCO')
    ForceControl_CO(fileparams,taskparams,forceparams,EMGparams);
end

%% EMG-control task
fileparams.code = '002';
fileparams.task = 'EMGCO';

fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

EMGparams.fchEMG = 30;
EMGparams.fclEMG = 60;
EMGparams.smoothWin = 800;
EMGparams.EMGScale = EMGparams.EMGScaleMVF;
EMGparams.EMGScaleType = 'MVF';
EMGparams.channelControl = [2 3];

taskparams.targetEMG = 1;
taskparams.targetAnglesEMG = sort([EMGparams.channelAngle(EMGparams.channelControl) mean(EMGparams.channelAngle(EMGparams.channelControl))]);
    
if strcmp(fileparams.task,'EMGCO') 
    EMGControl_CO(fileparams,taskparams,forceparams,EMGparams);
end

% Save params for each EMGCO file - different code
paramsfilename =  [fileparams.date,'_s',fileparams.subject,'_params_',fileparams.code,'.mat'];
save([fileparams.filepath,'Parameters/',paramsfilename],'taskparams','forceparams','EMGparams');
