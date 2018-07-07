% Main for DAQ from F/T sensor

% Parameters
fileparams = struct();
fileparams.savefile = 1;
fileparams.filename = '20180706_FS_test.txt';
fileparams.filepath = ['D:\Student_experiments\Virginia\FTSensor\Data\',fileparams.filename];

taskparams = struct();
taskparams.targetForce = 5; % [N]
taskparams.targetAngles = [0:pi/4:2*pi]; % [rad]
taskparams.rCirCursor = 0.25; % [N]
taskparams.rCirTarget = 0.5; % [N]

taskparams.movemtime = 5; % sec
taskparams.holdtime = 0.5; % sec
taskparams.timeout = 1; % sec
taskparams.relaxtime = 1; % sec

% Data acquisition
DataAcquisition(fileparams,taskparams);

