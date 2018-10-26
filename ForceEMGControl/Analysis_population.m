%% Data Analysis for population
addpath(genpath('Tools'));

switch computer
    case 'PCWIN64'
        filepath =  ['D:\Student_experiments\Virginia\Data\'];
    case 'MACI64'
        filepath =  ['/Users/virginia/Documents/MATLAB/Thesis/Data/'];
end

trial_ppl = struct();
h = 0;
datedir = dir(filepath);
datedir = datedir([datedir.isdir]);
for i = 3:length(datedir)
    subjectdir = dir([filepath,datedir(i).name]);
    subjectdir = subjectdir([subjectdir.isdir]);
    for j = 3:length(subjectdir)
        if exist([subjectdir(j).folder,'/',subjectdir(j).name,'/TD/'],'dir')
            load([subjectdir(j).folder,'/',subjectdir(j).name,'/TD/',datedir(i).name,'_',subjectdir(j).name,'trial_pp.mat']);
            h = h+1;
            trial_ppl(h) = trial_pp;
            %
            %         load([subjectdir(i).folder,'/',subjectdir(i).name,'/trial_data/trial_data_coh_EMG.mat']);
            %         trial_data_coh_EMG_pop = [trial_data_coh_EMG_pop, trial_data_coh_EMG];            
        end
    end
end

% tdc_force_pop_avg = trialAngleAvg(trial_data_coh_force_pop,[],{'rect.coh','rect.z'});
% tdc_EMG_pop_avg = trialAngleAvg(trial_data_coh_EMG_pop,[],{'rect.coh','rect.z'});

