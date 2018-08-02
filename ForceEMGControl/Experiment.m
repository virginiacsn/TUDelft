%% Experiment protocol
close all;
%clear all;

addpath('Tools');

% Parameters
fileparams = struct(...
    'saveforce',    0,...
    'saveEMG',      0,...
    'date',         '20180731',...
    'task',         'ForceCO',...
    'code',         '001');

fileparams.filenameforce =  [fileparams.date,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];
if ~exist(['D:\Student_experiments\Virginia\Data\' fileparams.date],'dir') && (fileparams.saveEMG || fileparams.saveforce)
    mkdir(['D:\Student_experiments\Virginia\Data\' fileparams.date])
end
fileparams.filepath =       ['D:\Student_experiments\Virginia\Data\' fileparams.date '\'];

EMGparams = struct(...
    'plotEMG',          0,...
    'EMGEnabled',       0,...
    'channelSubset',    [1 2 17],...
    'channelControl',   [1 2],...
    'channelName',      {{'BB','TL','Trigger'}},...
    'sampleRateEMG',    1024,... [samples/sec]
    'fchEMG',           10,... % [Hz]
    'smoothWin',        600,...
    'pauseSamp',        0.04,... % [s]
    'iterUpdatePlotEMG',1);

forceparams = struct(...
    'fclF',             2,... % [Hz]
    'scanRate',         2048,... % [scans/sec]
    'availSamples',     200,... % [samples]
    'availSamplesEMG',  200,... % [samples]
    'bufferWin',        500,... % [samples]
    'iterUpdatePlot',   2);

taskparams = struct(...
    'numTargets',       8,...
    'numTrials',        3,...
    'numTrialsEMG',     30,...
    'targetForce',      10,... % [N]
    'targetEMG',        0.2,...
    'tolTarget',        0.1,...
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
PreAparams.target_angles = [0:2*pi/8:2*pi-pi/8];
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

EMGparams.EMGScale = max(reshape([trial_data_avg.EMGrect],2,length(trial_data_avg)),[],2);

%% EMG-control task
fileparams.task = 'EMGCO';

fileparams.filenameforce =       [fileparams.date,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

if strcmp(fileparams.task,'EMGCO')
    taskparams.numTargets =     3;
    EMGControl_CO(fileparams,taskparams,forceparams,EMGparams);
end

%% Save params
if ~exist([fileparams.filepath,'Parameters/'],'dir')
    mkdir([fileparams.filepath,'Parameters/'])
end

paramsfilename =  [fileparams.date,'_params_',fileparams.code,'.mat'];
save([fileparams.filepath,'Parameters/',paramsfilename],'taskparams','forceparams','EMGparams');
