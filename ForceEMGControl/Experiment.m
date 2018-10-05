%% Experiment protocol
close all;
clear all;

addpath(genpath('Tools'));

% Parameters
fileparams = struct(...
    'saveforce',    1,...
    'saveEMG',      1,...
    'date',         '20181005',...
    'subject',      '04');

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
    'numTargetsForce',      7,...
    'numTargetsEMG',        3,...
    'targetForce',          8,... % [N]
    'targetForceCal',       30,... % [N]
    'targetEMG',            1,... % [% EMGScale]
    'targetEMGCal',         1,...
    'targetTolForce',       0.1,... % targetForce*targetTol
    'targetTolEMG',         0.2,... % targetEMG*targetTol
    'cursorTol',            2,...   % targetTol/cursorTol
    'targetAnglesForceCal', [pi/4:pi/2:7*pi/4],...
    'movemtime',            3,... % [sec]
    'holdtimeForce',        6,... % [sec]
    'holdtimeEMG',          5,... % [sec]
    'timeout',              1,... % [sec]
    'relaxtime',            2,... % [sec]
    'setFig',               0); % fig display for test/experiment
 
if taskparams.numTargetsForce == 4
    taskparams.targetAnglesForce = [pi/4:pi/2:7*pi/4]; % [rad]
elseif taskparams.numTargetsForce == 7
    taskparams.targetAnglesForce = [pi/4:pi/4:7*pi/4]; % [rad]
end

forceparams = struct(...
    'fclF',             2,... % [Hz]
    'scanRate',         2048,... % [scans/sec]
    'availSamples',     200,... % [samples]
    'availSamplesEMG',  200,... % [samples]
    'bufferWin',        500,... % [samples]
    'iterUpdatePlot',   2,...   % [iterations of listener call]
    'rotAnglePlot',     -pi/4); % [rad]

EMGparams = struct(...
    'plotEMG',          0,... % plot EMG 
    'channelSubset',    [1 2 3 4 5 6 7 8 17],... %,...[1 2 3 4 17] 
    'channelSubsetCal', [1 2 3 4 17],... %[1 2 3 4 17],...
    'channelName',      {{'BB','TLH','DA','DP','ECRB','FCR','Br','TLat','Trigger'}},...%{{'BB','TLH','DA','DP','Trigger'}},... 
    'channelNameCal',   {{'BB','TLH','DA','DP','Trigger'}},...%{{'BB','TLH','DA','DP','Trigger'}},...
    'channelAngle',     [5*pi/4,pi/4,3*pi/4,7*pi/4,0,0,0,0],...
    'channelAngleCal',  [5*pi/4,pi/4,3*pi/4,7*pi/4],...%[5*pi/4,pi/4,3*pi/4,7*pi/4],...
    'sampleRateEMG',    1024,... % [samples/sec]
    'fchEMG',           30,... % [Hz]
    'fclEMG',           60,... % [Hz]
    'fnEMG',            [],... % [Hz]
    'smoothWin',        800,... % [samples]
    'iterUpdatePlotEMG',1); % [iterations of listener call]
 
%% EMG offset test
EMGparams.EMGOffset = EMGOffsettest(EMGparams);
%EMGparams.EMGScaleMVC_start = MVCtest(EMGparams);

%% Force-control task
fileparams.code = '002';
fileparams.task = 'ForceCO';

fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

if strcmp(fileparams.task,'ForceCO')
    ForceControl_CO(fileparams,taskparams,forceparams,EMGparams);
end

%% Pre-analysis for EMGCO calibration
codeF = {'001','002'};

trial_data = [];

for i = 1:length(codeF)
    
    fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',codeF{i},'.mat'];
    fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',codeF{i},'.mat'];
    
    load([fileparams.filepath,fileparams.filenameforce]);
    load([fileparams.filepath,fileparams.filenameEMG]);
    
    if strcmp(fileparams.task,'ForceCO')
        forceEMGData = {forceDataOut_ForceCO,EMGDataOut_ForceCO};
        
        PreAparams.targetAngles = taskparams.targetAnglesForce;
    end
    
    PreAparams.fs = min(forceparams.scanRate,EMGparams.sampleRateEMG);
    PreAparams.downsamp = 2;
    PreAparams.fclF = forceparams.fclF;
    PreAparams.fchEMG = EMGparams.fchEMG;
    PreAparams.fclEMG = EMGparams.fclEMG;
    PreAparams.fnEMG = EMGparams.fnEMG;
    PreAparams.avgWindow = 200;
    
    trial_data = trialCO(forceEMGData,PreAparams);
    trial_data = removeFailTrials(trial_data(5:end));
    
    PreAparams.targetAngles = sort(unique([trial_data.angle]));
    Aparams.targetAnglesForce = PreAparams.targetAngles;
    
    trial_data = procEMG(trial_data,PreAparams);
    trial_data = [trial_data, procForce(trial_data,PreAparams)];
end

epoch = {'ihold',1,'iend',0};
fields = {'EMG.rect','force.filtmag'};
trial_data_avg = trialAngleAvg(trial_data, epoch, fields);

EMGmean = zeros(length(trial_data_avg),length(EMGparams.channelSubset)-1);
forcemean = zeros(length(trial_data_avg),1);
for i = 1:length(trial_data_avg)
    EMGmean(i,:) = mean(trial_data_avg(i).EMG.rect,1);
    forcemean(i) = trial_data_avg(i).force.filtmag_mean;
end

EMGparams.EMGScaleForce = max(EMGmean,[],1)'; 
%EMGparams.EMGScaleForce = mean(EMGmean,1)';

fprintf('\nEMG scaling values: \n')
for k = 1:length(EMGparams.channelSubsetCal)-1
    fprintf('%s: %1.3f\n',EMGparams.channelName{EMGparams.channelSubset == EMGparams.channelSubsetCal(k)},EMGparams.EMGScaleForce(EMGparams.channelSubset == EMGparams.channelSubsetCal(k)))
end

fprintf('\nRecorded mean force value: %1.3f\n',round(mean(forcemean)))

% Saving trial data for force-control for comparison
PreAparams.fclF = 5;
PreAparams.fchEMG = 20;
PreAparams.fclEMG = 500;
PreAparams.fcnEMG = [];

trial_data = procEMG(trial_data,PreAparams);
trial_data = procForce(trial_data,PreAparams);

epoch = {'ihold',1,'iend',0};
fields = {'EMG.rect','EMG.avg','force.filt'};
trial_data_avg_force = trialAngleAvg(trial_data, epoch, fields);

% taskparams.targetForce = round(mean(forcemean))*0.5;

% Deleting variables
if strcmp(fileparams.task,'ForceCO')
    clear forceEMGData forceDataOut_ForceCO EMGDataOut_ForceCO
end

%% Calibration with EMG-control task
fileparams.code = 'calib';
fileparams.task = 'EMGCO';

EMGparams.EMGScale = EMGparams.EMGScaleForce;

fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

EMGControl_Cal(fileparams,taskparams,forceparams,EMGparams);

%% Pre-analysis for EMG comparison Force-control and EMG-control
load([fileparams.filepath,fileparams.filenameforce]);
load([fileparams.filepath,fileparams.filenameEMG]);

if strcmp(fileparams.task,'EMGCO')
    forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};
 
    PreAparams.targetAngles = sort(EMGparams.channelAngleCal);
end

trial_data = trialCO(forceEMGData,PreAparams);
trial_data = removeFailTrials(trial_data(5:end));

PreAparams.targetAngles = sort(unique([trial_data.angle]));
Aparams.targetAnglesEMG = PreAparams.targetAngles;

trial_data = procEMG(trial_data,PreAparams);
trial_data = procForce(trial_data,PreAparams);

epoch = {'ihold',1,'iend',0};
fields = {'EMG.rect','EMG.avg','force.filt'};
trial_data_avg_EMG = trialAngleAvg(trial_data, epoch, fields);

% Deleting variables
if strcmp(fileparams.task,'EMGCO')
    clear forceEMGData forceDataOut_EMGCO EMGDataOut_EMGCO
end

%% Check force and EMG for EMG-control calibration task
% Limits for force and EMG plots
Flim = 2; EMGlim = 100;

plotForceEMGtimeComp;

%% CONTINUE? [Y/N]
%% Force-control task 
fileparams.code = '003';
fileparams.task = 'ForceCO';

fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

if strcmp(fileparams.task,'ForceCO')
    ForceControl_CO(fileparams,taskparams,forceparams,EMGparams);
end

%% EMG-control task
fileparams.code = '001';
fileparams.task = 'EMGCO';

fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

EMGparams.fchEMG = 30;
EMGparams.fclEMG = 60;
EMGparams.fnEMG = [];
EMGparams.smoothWin = 800;
EMGparams.EMGScale = EMGparams.EMGScaleForce; %EMGparams.EMGScaleMVC_start(:,1);
EMGparams.EMGScaleType = 'Force';
EMGparams.channelControl = [1 3];

taskparams.numTargetsEMG = 3;
taskparams.targetEMG = 1;
taskparams.targetAnglesEMG = sort([EMGparams.channelAngle(EMGparams.channelControl) mean(EMGparams.channelAngle(EMGparams.channelControl))]);
    
if strcmp(fileparams.task,'EMGCO') 
    EMGControl_CO(fileparams,taskparams,forceparams,EMGparams);
end

% Save params for each EMGCO file - different code
paramsfilename =  [fileparams.date,'_s',fileparams.subject,'_params_',fileparams.code,'.mat'];
save([fileparams.filepath,'Parameters/',paramsfilename],'taskparams','forceparams','EMGparams');

%% End MVC test
EMGparams.EMGScaleMVC_end = MVCtest(EMGparams);

% Save params for each EMGCO file (save end MVC test in last params code)
paramsfilename =  [fileparams.date,'_s',fileparams.subject,'_params_',fileparams.code,'.mat'];
save([fileparams.filepath,'Parameters/',paramsfilename],'taskparams','forceparams','EMGparams');
