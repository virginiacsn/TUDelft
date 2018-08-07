%% Experiment protocol
close all;
%clear all;

addpath(genpath('Tools'));

% Parameters
fileparams = struct(...
    'saveforce',    1,...
    'saveEMG',      1,...
    'date',         '20180807',...
    'subject',      '01',...
    'task',         'ForceCO',...
    'code',         '001');

fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];
if ~exist(['D:\Student_experiments\Virginia\Data\',fileparams.date],'dir') && (fileparams.saveEMG || fileparams.saveforce)
    mkdir(['D:\Student_experiments\Virginia\Data\',fileparams.date])
end
if ~exist(['D:\Student_experiments\Virginia\Data\',fileparams.date,'\s',fileparams.subject],'dir') && (fileparams.saveEMG || fileparams.saveforce)
    mkdir(['D:\Student_experiments\Virginia\Data\',fileparams.date,'\s',fileparams.subject])
end
fileparams.filepath =       ['D:\Student_experiments\Virginia\Data\',fileparams.date,'\s',fileparams.subject,'\'];

EMGparams = struct(...
    'plotEMG',          0,...
    'EMGEnabled',       0,...
    'channelSubset',    [1 2 3 17],...
    'channelControl',   [1 2],...
    'channelName',      {{'BB','TL','TLH','Trigger'}},...
    'sampleRateEMG',    1024,... % [samples/sec]
    'fchEMG',           10,... % [Hz]
    'fclEMG',           30,... % [Hz]
    'smoothWin',        600,... % [samples]
    'pauseSamp',        0.04,... % [s]
    'iterUpdatePlotEMG',1);

forceparams = struct(...
    'fclF',             2,... % [Hz]
    'scanRate',         2048,... % [scans/sec]
    'availSamples',     200,... % [samples]
    'availSamplesEMG',  200,... % [samples]
    'bufferWin',        500,... % [samples]
    'iterUpdatePlot',   2,...
    'rotAngPlot',       -45); % [deg]

taskparams = struct(...
    'numTargets',       4,...
    'numTargetsEMG',    3,...
    'targetAnglesForce',[pi/4:pi/2:7*pi/4],...
    'targetAnglesEMG',  [pi/4:pi/4:3*pi/4],... % [rad]
    'numTrials',        3,...
    'numTrialsEMG',     30,...
    'targetForce',      10,... % [N]
    'targetEMG',        0.2,...
    'targetTol',        0.1,...
    'targetTolEMG',     0.2,...
    'cursorTol',        1.5,...
    'movemtime',        5,... % sec
    'holdtime',         3,... % sec
    'timeout',          1,... % sec
    'relaxtime',        2); % sec

%% Force-control task
if strcmp(fileparams.task,'ForceCO')
    ForceControl_CO(fileparams,taskparams,forceparams,EMGparams);
end

%% Pre-analysis
load([fileparams.filepath,fileparams.filenameforce]);
load([fileparams.filepath,fileparams.filenameEMG]);

forceEMGData = {forceDataOut_ForceCO,EMGDataOut_ForceCO};

PreAparams.downsample = 2;
PreAparams.target_angles = taskparams.targetAnglesForce;
PreAparams.avgWindow = 200;
PreAparams.fclF = 5;
PreAparams.fchEMG = 10;

trial_data = trialCO(forceEMGData,PreAparams);

trial_data = removeFailTrials(trial_data);

trial_data = procEMG(trial_data,PreAparams);
trial_data = procForce(trial_data,PreAparams);

epoch = {'ihold','iend'};
fields = {'EMG.rect'};
trial_data_avg = trialAngleAvg(trial_data, epoch, fields);

EMGmean = zeros(length(trial_data_avg),length(EMGparams.channelControl));
for i = 1:length(trial_data_avg)
    EMGmean(i,:) = trial_data_avg(i).EMG.rect_mean(EMGparams.channelControl);
end

EMGparams.EMGScale = max(EMGmean,[],1)';

%% EMG-control task
fileparams.task = 'EMGCO';

fileparams.filenameforce =  [fileparams.date,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

if strcmp(fileparams.task,'EMGCO') 
    EMGControl_CO(fileparams,taskparams,forceparams,EMGparams);
end

%% Save params
if ~exist([fileparams.filepath,'Parameters/'],'dir')
    mkdir([fileparams.filepath,'Parameters/'])
end

paramsfilename =  [fileparams.date,'_s',fileparams.subject,'_params_',fileparams.code,'.mat'];
save([fileparams.filepath,'Parameters/',paramsfilename],'taskparams','forceparams','EMGparams');
