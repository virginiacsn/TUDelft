%% Main for DAQ from F/T sensor
addpath('Tools');

% Parameters
fileparams =    struct(...
    'saveforce', 1,...
    'saveEMG', 1,...
    'date', '20180717',...
    'task', 'EMGCO',...
    'code', '001');

if fileparams.saveEMG
    fileparams.saveformat = '.mat';
else
     fileparams.saveformat = '.txt';
end

fileparams.filenameforce =       [fileparams.date,'_',fileparams.task,'_',fileparams.code,fileparams.saveformat];
fileparams.filenameEMG =    [fileparams.date,'_',fileparams.task,'_EMG_',fileparams.code,fileparams.saveformat];
if ~exist([pwd '\Data\' fileparams.date],'dir')
    mkdir([pwd '\Data\' fileparams.date])
end
fileparams.filepath =       [pwd '\Data\' fileparams.date '\'];

EMGparams =     struct(...
    'EMG', 1,...
    'plotEMG', 0,...
    'EMGEnabled', 0,...
    'channel_subset', [1 2],...
    'channel_name', {{'BB','TBL'}},...
    'sample_rate', 1024);

taskparams =    struct(...
    'fc', 3,... % [Hz]
    'scanRate', 2048,... % [scans/sec]
    'availSamples', 200,... % [samples]
    'movemtime', 3,... % sec
    'holdtime', 1,... % sec
    'timeout', 1,... % sec
    'relaxtime', 2,... % sec
    'targetForce', 5000,... % [N]
    'numTargets', 3);

taskparams.rCirTarget =     taskparams.targetForce/10; % [N]
taskparams.rCirCursor =     taskparams.targetForce/20; % [N]

% Data acquisition
if strcmp(fileparams.task,'BAR')
    
    DAQ_ForceBAR(fileparams,taskparams);
    
elseif strcmp(fileparams.task,'CO')
 
    
    if EMGparams.EMG
        DAQ_ForceCO_EMG(fileparams,taskparams,EMGparams);
    else
        DAQ_ForceCO(fileparams,taskparams);
    end
elseif strcmp(fileparams.task,'EMGCO')
    DAQ_EMGCO(fileparams,taskparams,EMGparams);
end

