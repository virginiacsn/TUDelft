%% Main for DAQ from F/T sensor
addpath('Tools');

% Parameters
fileparams = struct(...
    'saveforce',    1,...
    'saveEMG',      1,...
    'date',         '20180724',...
    'task',         'ForceCO',...
    'code',         '002');

if fileparams.saveEMG
    fileparams.saveformat = '.mat';
else
    fileparams.saveformat = '.txt';
end

fileparams.filenameforce =       [fileparams.date,'_',fileparams.task,'_Force_',fileparams.code,fileparams.saveformat];
fileparams.filenameEMG =    [fileparams.date,'_',fileparams.task,'_EMG_',fileparams.code,fileparams.saveformat];
if ~exist(['D:\Student_experiments\Virginia\Data\' fileparams.date],'dir')
    mkdir(['D:\Student_experiments\Virginia\Data\' fileparams.date])
end
fileparams.filepath =       ['D:\Student_experiments\Virginia\Data\' fileparams.date '\'];

EMGparams = struct(...
    'EMG',              1,...
    'plotEMG',          0,...
    'EMGEnabled',       0,...
    'channelSubset',    [1 2 17],...
    'channelName',      {{'BB','TL','Trigger'}},...
    'sampleRateEMG',    1024,...
    'fclH',             5,...
    'smoothWin',        600,...
    'pauseSamp',        0.04,...
    'iterUpdatePlotEMG',1);

forceparams = struct(...
    'fclF',             2,... % [Hz]
    'scanRate',         2048,... % [scans/sec]
    'availSamples',     200,... % [samples]
    'bufferWin',        500,... % [samples]
    'iterUpdatePlot',   2);

taskparams = struct(...
    'numTargets',       8,...
    'numTrials',        3,...
    'movemtime',        5,... % sec
    'targetForce',      5,... % [N]
    'holdtime',         5,... % sec
    'timeout',          1,... % sec
    'relaxtime',        2 ); % sec

taskparams.rCirTarget =     taskparams.targetForce/10; % [N]
taskparams.rCirCursor =     taskparams.targetForce/20; % [N]

% Data acquisition
if strcmp(fileparams.task,'BAR')
    DAQ_ForceBAR(fileparams,taskparams,forceparams);
elseif strcmp(fileparams.task,'ForceCO')
    if EMGparams.EMG
        ForceControl_CO(fileparams,taskparams,forceparams,EMGparams);
    else
        ForceControl_CO(fileparams,taskparams,forceparams);
    end
elseif strcmp(fileparams.task,'EMGCO')
    taskparams.numTargets = 3;
    EMGControl_CO(fileparams,taskparams,forceparams,EMGparams);
end
