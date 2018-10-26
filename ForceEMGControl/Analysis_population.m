%% Data Analysis for population
% clear all;
addpath(genpath('Tools'));

switch computer
    case 'PCWIN64'
        filepath =  ['D:\Student_experiments\Virginia\Data\TD\'];
    case 'MACI64'
        filepath =  ['/Users/virginia/Documents/MATLAB/Thesis/Data/TD/'];
end

Aparams_pp = [];
trial_pp_force = [];
trial_pp_EMG = [];

dirs = dir(filepath);
files = dirs(~[dirs.isdir]);
for i = 1:length(files)
    load([filepath,files(i).name]);
    Aparams_pp = [Aparams_pp, trial_pp.Aparams];
    trial_pp_force = [trial_pp_force, trial_pp.forceCO];
    trial_pp_EMG = [trial_pp_EMG, trial_pp.EMGCO];
end

