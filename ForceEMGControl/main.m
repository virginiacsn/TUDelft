%% Main for DAQ from F/T sensor
addpath('Tools');

% Parameters
fileparams =    struct(...
    'saveforce', 0,...
    'saveEMG', 1,...
    'date', '20180719',...
    'task', 'EMGCO',...
    'code', '001');

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


EMGparams =     struct(...
    'EMG', 1,...
    'plotEMG', 0,...
    'EMGEnabled', 0,...
    'channelSubset', [1 2],...
    'channelName', {{'BB','TBL'}},...
    'sampleRateEMG', 1024,...
    'fclH', 10,...
    'smoothWin', 800,...
    'pauseSamp', 0.01);

taskparams =    struct(...
    'fclF', 2,... % [Hz]
    'scanRate', 2048,... % [scans/sec]
    'availSamples', 200,... % [samples]
    'bufferWin', 500,... % [samples]
    'iterUpdatePlot', 2,...
    'movemtime', 3,... % sec
    'holdtime', 1,... % sec
    'timeout', 1,... % sec
    'relaxtime', 2,... % sec
    'targetForce', 0.1,... % [N]
    'numTargets', 3,...
    'numTrials', 3);

taskparams.rCirTarget =     taskparams.targetForce/10; % [N]
taskparams.rCirCursor =     taskparams.targetForce/20; % [N]

% Data acquisition
if strcmp(fileparams.task,'BAR')  
    DAQ_ForceBAR(fileparams,taskparams);
elseif strcmp(fileparams.task,'ForceCO')
    if EMGparams.EMG
        ForceControl_CO(fileparams,taskparams,EMGparams);
    else
        ForceControl_CO(fileparams,taskparams);
    end
elseif strcmp(fileparams.task,'EMGCO')
    EMGControl_CO(fileparams,taskparams,EMGparams);
end

