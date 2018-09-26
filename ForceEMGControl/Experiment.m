%% Experiment protocol
close all;
%clear all;

addpath(genpath('Tools'));

% Parameters
fileparams = struct(...
    'saveforce',    1,...
    'saveEMG',      1,...
    'date',         '20180926',...
    'subject',      '02');

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
    'targetForce',          10,... % [N]
    'targetForceCal',       30,... % [N]
    'targetEMG',            0.2,... % [% EMGScale]
    'targetEMGCal',         40,...
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
    'channelSubset',    [1 2 3 4 5 6 7 8 17],...
    'channelSubsetCal', [1 2 3 4 17],...
    'channelName',      {{'BB','TLH','DA','DP','ECRB','FCR','Br','TLat','Trigger'}},...
    'channelNameCal',   {{'BB','TLH','DA','DP','Trigger'}},...
    'channelAngle',     [5*pi/4,pi/4,3*pi/4,7*pi/4,0,0,0,0],...
    'channelAngleCal',  [5*pi/4,pi/4,3*pi/4,7*pi/4],...
    'sampleRateEMG',    1024,... % [samples/sec]
    'fchEMG',           30,... % [Hz]
    'fclEMG',           60,... % [Hz]
    'smoothWin',        800,... % [samples]
    'iterUpdatePlotEMG',1); % [iterations of listener call]

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

%% Calibration with MVC task  - start MVC test
EMGparams.EMGScaleMVC_start = MVCtest(EMGparams);

%% Pre-analysis for calibration
load([fileparams.filepath,fileparams.filenameforce]);
load([fileparams.filepath,fileparams.filenameEMG]);

if strcmp(fileparams.task,'ForceCO')
    forceEMGData = {forceDataOut_ForceCO,EMGDataOut_ForceCO};
    
    PreAparams.targetAngles = taskparams.targetAnglesForce;
elseif strcmp(fileparams.task,'EMGCO')
    forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};
    
    PreAparams.targetAngles = sort(EMGparams.channelAngle);
end

PreAparams.downsample = 2;
PreAparams.avgWindow = 200;
PreAparams.fclF = forceparams.fclF;
PreAparams.fchEMG = 20;

trial_data = trialCO(forceEMGData,PreAparams);

trial_data = removeFailTrials(trial_data);

PreAparams.targetAngles = sort(unique([trial_data.angle]));


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

taskparams.targetForce = round(mean(forcemean))*0.5;

if strcmp(fileparams.task,'ForceCO')
    clear forceEMGData forceDataOut_ForceCO EMGDataOut_ForceCO
elseif strcmp(fileparams.task,'EMGCO')
    clear forceEMGData forceDataOut_EMGCO EMGDataOut_EMGCO
end

%% Check force and EMG 
plotForceEMGtime;

%% Force-control task
fileparams.code = '002';
fileparams.task = 'ForceCO';

fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

if strcmp(fileparams.task,'ForceCO')
    ForceControl_CO(fileparams,taskparams,forceparams,EMGparams);
end

%% EMG-control task
fileparams.code = '006';
fileparams.task = 'EMGCO';

fileparams.filenameforce =  [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_Force_',fileparams.code,'.mat'];
fileparams.filenameEMG =    [fileparams.date,'_s',fileparams.subject,'_',fileparams.task,'_EMG_',fileparams.code,'.mat'];

EMGparams.fchEMG = 30;
EMGparams.fclEMG = 60;
EMGparams.smoothWin = 800;
EMGparams.EMGScale = EMGparams.EMGScaleMVF; %EMGparams.EMGScaleMVC_start(:,1);
EMGparams.EMGScaleType = 'MVF';
EMGparams.channelControl = [1 3];

taskparams.numTargetsEMG = 3;
taskparams.targetEMG = 0.6;
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
