%% Experiment protocol
close all;
%clear all;

addpath('Tools');

% Parameters
fileparams = struct(...
    'saveforce',    1,...
    'saveEMG',      1,...
    'date',         '20180726',...
    'task',         'ForceCO',...
    'code',         '002');

fileparams.filenameforce =       [fileparams.date,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];
if ~exist(['D:\Student_experiments\Virginia\Data\' fileparams.date],'dir')
    mkdir(['D:\Student_experiments\Virginia\Data\' fileparams.date])
end
fileparams.filepath =       ['D:\Student_experiments\Virginia\Data\' fileparams.date '\'];

EMGparams = struct(...
    'plotEMG',          0,...
    'EMGEnabled',       0,...
    'channelSubset',    [1 2 17],...
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
    'movemtime',        5,... % sec
    'targetForce',      10,... % [N]
    'holdtime',         1,... % sec
    'timeout',          1,... % sec
    'relaxtime',        2 ); % sec

taskparams.rCirTarget =     taskparams.targetForce/10; % [N]
taskparams.rCirCursor =     taskparams.targetForce/20; % [N]

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
fields = {'EMGrect'};
trial_data_avg = trialangleavg(trial_data, epoch, fields);

EMGparams.EMGScale = max(reshape([trial_data_avg.EMGrect],2,length(trial_data_avg)),[],2);

%% EMG-control task
fileparams.task = 'EMGCO';

fileparams.filenameforce =       [fileparams.date,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

if strcmp(fileparams.task,'EMGCO')
    taskparams.numTargets =     3;
    taskparams.targetForce =    0.2;
    taskparams.rCirTarget =     taskparams.targetForce/10; % [N]
    taskparams.rCirCursor =     taskparams.targetForce/20; % [N]
    
    EMGControl_CO(fileparams,taskparams,forceparams,EMGparams);
end

%% Save params

if ~exist([fileparams.filepath,'Parameters/'],'dir')
    mkdir([fileparams.filepath,'Parameters/'])
end

paramsfilename =  [fileparams.date,'_params_',fileparams.code,'.mat'];
save([fileparams.filepath,'Parameters/',paramsfilename],'taskparams','forceparams','EMGparams');
