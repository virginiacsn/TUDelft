%% Data analysis for individual subject
% close all
clear all
addpath(genpath('Tools'));

date =      '20181010';
subject =   '06'; 

savepp = 1;
 
switch computer
    case 'PCWIN'
        datapath = 'D:\Student_experiments\Virginia\Data\';
        filepath =  [datapath,'Exp\',date,'\s',subject,'\'];
        paramfolder = 'Parameters\';
    case 'PCWIN64'
        datapath = 'D:\Student_experiments\Virginia\Data\';
        filepath =  [datapath,'Exp\',date,'\s',subject,'\'];
        paramfolder = 'Parameters\';
        if ~exist([datapath,'TD'],'dir') && exist(datapath,'dir')
            mkdir([datapath,'TD'])
        end
        filepathpp = [datapath,'TD\'];
    case 'MACI64'
        datapath = '/Users/virginia/Documents/MATLAB/Thesis/Data/';
        filepath =  [datapath,'Exp/',date,'/s',subject,'/'];
        paramfolder = 'Parameters/';
        if ~exist([datapath,'TD'],'dir') && exist(datapath,'dir')
            mkdir([datapath,'TD'])
        end
        filepathpp = [datapath,'TD/'];
end

%% DATA LOADING AND PREPROC
%% All blocks for force-control task, blocks corresponding to each muscle control pair for EMG-control
codeF = {'001','002','003','004'};

% Force-control
task =  'ForceCO';
filenameparams = [date,'_s',subject,'_params_','001','.mat'];
load([filepath,paramfolder,filenameparams]);

% Trial number to start analysis. Leave 10 initial trials for learning task
startTrial = 10;

% Parameters for analysis
Aparams.fs = min(forceparams.scanRate,EMGparams.sampleRateEMG);
Aparams.downsamp = forceparams.scanRate/EMGparams.sampleRateEMG;
Aparams.channelNameEMG = EMGparams.channelName;
Aparams.targetAngles = taskparams.targetAnglesForce;
Aparams.fclF = 5;
Aparams.fchEMG = 10;
Aparams.fclEMG = 500;
Aparams.fnEMG = 50;
Aparams.avgWindow = 200;

% Epoch interval {epoch start, time start, epoch end, time end} and fields to trial average
Aparams.epoch = {'ihold',1,'iend',0};
fields_avg = {'EMG.raw','EMG.filt','EMG.rect','EMG.avg','force.filt','force.filtmag'};
% Fields to trial append
fields_app = {'EMG.raw','EMG.filt','EMG.rect','force.filt'};

% Trial data struct for force-control task, will append trial data for each
% block (code)
trial_data_force = [];

for i = 1:length(codeF)
    
    code = codeF{i};
    
    filenameforce =  [date,'_s',subject,'_',task,'_Force_',code,'.mat'];
    filenameEMG = [date,'_s',subject,'_',task,'_EMG_',code,'.mat'];
    
    if exist([filepath,filenameforce],'file')&&exist([filepath,filenameEMG],'file')
        load([filepath,filenameforce]);
        load([filepath,filenameEMG]);
        
        forceEMGData = {forceDataOut_ForceCO,EMGDataOut_ForceCO};
        
        Aparams.block = str2double(code);
        if exist('EMGOffset','var')
            Aparams.EMGOffset = EMGOffset;
        end
        
        trial_data = trialCO(forceEMGData,Aparams);
        
        trial_data = removeFailTrials(trial_data(startTrial:end));
        
        trial_data = procEMG(trial_data,Aparams);
        trial_data_force = [trial_data_force, procForce(trial_data,Aparams)];
        
        % Check EMG offset values for each file
        if exist('EMGOffset','var')
            fprintf(['\nEMG offset values (code ',code,'): \n'])
            for k = 1:length(EMGparams.channelSubsetCal)-1
                fprintf('%s: %1.3f\n',EMGparams.channelName{EMGparams.channelSubset == EMGparams.channelSubsetCal(k)},Aparams.EMGOffset(EMGparams.channelSubset == EMGparams.channelSubsetCal(k)))
            end
            fprintf('\n');
        end
    end
end

% Final force-control target angles based on successful trials
Aparams.targetAnglesForce = sort(unique(extractfield(trial_data_force,'angle')));

% Trial average by angle
trial_avg_force = trialAngleAvg(trial_data_force, Aparams.epoch, fields_avg);
% FFT of averaged EMG data
trial_avg_force = fftEMG(trial_avg_force);

% Select number of trials to append by finding minimum number of trials per
% target
mintrialsForce = min([trial_avg_force.ntrials]);

% Trial append by angle. Limit to 20 trials (100 sec of data)
trial_app_force = trialAngleApp(trial_data_force, Aparams.epoch, fields_app,[],20);

% Compute mean rectified EMG and filtered force magnitude for each angle.
% Check calibration values
EMGmean = zeros(length(trial_avg_force),length(EMGparams.channelSubset)-1);
forcemean = zeros(length(trial_avg_force),1);
for i = 1:length(trial_avg_force)
    EMGmean(i,:) = mean(trial_avg_force(i).EMG.rect,1);
    forcemean(i) = trial_avg_force(i).force.filtmag_mean;
end
EMGmeanForce = EMGmean;

% EMG scaling through max EMG for force-control
EMGScaleForce = max(EMGmean,[],1)';
% EMG target tolerance values
EMGTolForce = (min(EMGmean,[],1))'./EMGScaleForce;

fprintf('EMG mean values: \n')
for k = 1:length(EMGparams.channelSubsetCal)-1
    fprintf('%s: %1.3f\n',EMGparams.channelName{EMGparams.channelSubset == EMGparams.channelSubsetCal(k)},EMGScaleForce(EMGparams.channelSubset == EMGparams.channelSubsetCal(k)))
end

% Mean filtered force magnitude value
fprintf('\nRecorded mean force value: %1.3f\n',round(mean(forcemean)))

%% EMG-control
task = 'EMGCO';
codeE = {'001','002','003','004','005','006'};

% Trial number to start analysis. Leave 10 initial trials for learning task
startTrial = 5;

% Trial data struct for EMG-control task, will append trial data for each
% block (code)
trial_data_EMG = [];

for i = 1:length(codeE)
    
    code = codeE{i};
    
    filenameforce =  [date,'_s',subject,'_',task,'_Force_',code,'.mat'];
    filenameEMG = [date,'_s',subject,'_',task,'_EMG_',code,'.mat'];
    filenameparams = [date,'_s',subject,'_params_',code,'.mat'];
    
    load([filepath,filenameforce]);
    load([filepath,filenameEMG]);
    load([filepath,paramfolder,filenameparams]);
    
    forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};

    Aparams.block = str2double(code);
    Aparams.targetAngles = taskparams.targetAnglesEMG;
    Aparams.chanControlPair{i} = EMGparams.channelControl;
    if exist('EMGOffset','var')
        if size(EMGOffset,1) < size(EMGOffset,2)
            Aparams.EMGOffset = EMGOffset;
        else
            Aparams.EMGOffset = EMGOffset';
        end
    end
    
    trial_data = trialCO(forceEMGData,Aparams);
    
    trial_data = removeFailTrials(trial_data(startTrial:end));

    trial_data = procEMG(trial_data,Aparams);
    trial_data_EMG = [trial_data_EMG, procForce(trial_data,Aparams)];
    
    if exist('EMGOffset','var')
        % Check EMG offset values for each file
        if length(Aparams.EMGOffset) == length(EMGparams.channelSubset)-1
            fprintf(['\nEMG offset values (code ',code,'): \n']);
            for k = 1:length(EMGparams.channelSubsetCal)-1
                fprintf('%s: %1.3f\n',EMGparams.channelName{EMGparams.channelSubset == EMGparams.channelSubsetCal(k)},Aparams.EMGOffset(EMGparams.channelSubset == EMGparams.channelSubsetCal(k)))
            end
            fprintf('\n');
        end
    end
end

% Final EMG-control target angles based on successful trials
Aparams.targetAnglesEMG = sort(unique(extractfield(trial_data_EMG,'angle')));

% Trial average by angle
trial_avg_EMG = trialAngleAvg(trial_data_EMG, Aparams.epoch, fields_avg);
% FFT of averaged EMG data
trial_avg_EMG = fftEMG(trial_avg_EMG);

% Select number of trials to append by finding minimum number of trials per
% target
mintrialsEMG = min([trial_avg_EMG.ntrials]);

% Trial append by angle. Limit to 25 trials (100 sec of data)
trial_app_EMG = trialAngleApp(trial_data_EMG, Aparams.epoch, fields_app,[],25);

% Compute mean rectified EMG and filtered force magnitude for each angle.
% Check calibration values
EMGmean = zeros(length(trial_avg_EMG),length(EMGparams.channelSubset)-1);
forcemean = zeros(length(trial_avg_EMG),1);
for i = 1:length(trial_avg_EMG)
    EMGmean(i,:) = mean(trial_avg_EMG(i).EMG.rect,1);
    forcemean(i) = trial_avg_EMG(i).force.filtmag_mean;
end
EMGmeanEMG = EMGmean;

% Angles and corresponding muscle names to compare between tasks
[muscAngle,ichan] = sort(EMGparams.channelAngle(EMGparams.channelAngle>0));
muscAngles = sort([muscAngle,mean([muscAngle(1:end-1);muscAngle(2:end)])]);

chanControl = EMGparams.channelSubset(EMGparams.channelAngle>0);
channelName = EMGparams.channelName(chanControl(ichan));
Aparams.chanControl = chanControl(ichan);
Aparams.chanControlName = channelName;

channelNames = {};
k = 1;
for i = 1:length(channelName)
    channelNames{k} = channelName{i};
    if i ~= length(channelName)
        channelNames{k+1} = [channelName{i},',',channelName{i+1}];
    end
    k = k+2;
end

% Target angles to compare
angComp = muscAngles(ismember(muscAngles,Aparams.targetAnglesEMG));
Aparams.angCompUni = angComp;
Aparams.angComp = {};
for i = 1:length(angComp)
    Aparams.angComp{i} = angComp(i)*[1 1];
end

% Muscles names to compare
Aparams.muscComp = channelNames(ismember(muscAngles,Aparams.targetAnglesEMG));

% Compare by control muscle pair
siP = [0 0];
j = 0;
for i = 1:length(Aparams.chanControlPair)
    iC1 = find(angComp==EMGparams.channelAngle(Aparams.chanControlPair{i}(1)));
    iC2 = find(angComp==EMGparams.channelAngle(Aparams.chanControlPair{i}(2)));
    siC = sort([iC1,iC2]);
    if sum(siP == siC) == 2
        j = j+1;
        Aparams.angCompPair{j} = angComp(siC(1):siC(2));
        Aparams.muscCompPair{j} = Aparams.muscComp(siC(1):siC(2));
    end
    siP = siC;
end

%% LIMITS for force and EMG plots
Flim = round(max(forcemean))/taskparams.targetForce; EMGlim = 60;

% EMG_lim(i,j) will be max EMG level between force-control and EMG-control
% where rows are target angles to compare and cols muscles (all)
EMG_lim = [];
for i = 1:length(Aparams.angComp)
    for j = 1:length(EMGparams.channelName)-1
        EMG_lim(i,j) = max(max(trial_avg_force(Aparams.targetAnglesForce == Aparams.angComp{i}(1)).EMG.rect(:,j)),...
            max(trial_avg_EMG(Aparams.targetAnglesEMG == Aparams.angComp{i}(2)).EMG.rect(:,j)));
    end
end

% FFT lim
fc = 100;

% Limits for fft of rectified EMG
fft_fields = {'raw','filt','rect'};
for h = 1:length(fft_fields)
    for i = 1:length(Aparams.angComp)
        for j = 1:length(EMGparams.channelName)-1
            fft_lim.(fft_fields{h})(i,j) = max([max(abs(trial_avg_force...
                (Aparams.targetAnglesForce == Aparams.angComp{i}(1)).EMG.(['fft',fft_fields{h}])(3:end-2,j))),...
                max(abs(trial_avg_EMG(Aparams.targetAnglesEMG == Aparams.angComp{i}(2)).EMG.(['fft',fft_fields{h}])(3:end-2,j)))]);
        end
    end
end

%% COHERENCE ANALYSIS
Aparams.cohparams.data = 'app';
Aparams.cohparams.tseg = 1;
Aparams.cohparams.nseg = 10;
Aparams.cohparams.my_nseg = 10;
Aparams.cohparams.window = @(N) hanning(N);
Aparams.cohparams.CLoverlap = 1;
fields_coh = {'filt','rect'};
fc = 100;

if strcmp(Aparams.cohparams.data,'avg')
    trial_coh_force = cohStruct(trial_avg_force,EMGparams.channelName,fields_coh,Aparams.cohparams);
    trial_coh_EMG = cohStruct(trial_avg_EMG,EMGparams.channelName,fields_coh,Aparams.cohparams);
else
    trial_coh_force = cohStruct(trial_app_force,EMGparams.channelName,fields_coh,Aparams.cohparams);
    trial_coh_EMG = cohStruct(trial_app_EMG,EMGparams.channelName,fields_coh,Aparams.cohparams);
end

%% Save vars for population average
if savepp
    fprintf('Finished coherence analysis.\nSaving data.\n');
    trial_pp.Aparams = Aparams;
    trial_pp.Aparams.EMGScale = EMGparams.EMGScale';
    
    fstruct = [trial_avg_force.force];
    EMGstruct = [trial_avg_force.EMG];
    trial_pp.forceCO.angle = Aparams.targetAnglesForce';
    trial_pp.forceCO.force.mag_mean = [fstruct.filtmag_mean]';
    trial_pp.forceCO.force.mag_std = [fstruct.filtmag_std]';
    trial_pp.forceCO.force.mag_pstd = [fstruct.filtmag_pstd]';
    trial_pp.forceCO.force.mag_sem = [fstruct.filtmag_sem]';
    trial_pp.forceCO.force.mag_psem = [fstruct.filtmag_psem]';
    trial_pp.forceCO.force.filt_mean = [fstruct.filt_mean]';
    trial_pp.forceCO.force.filt_std = [fstruct.filt_std]';
    trial_pp.forceCO.force.filt_pstd = [fstruct.filt_pstd]';
    trial_pp.forceCO.force.filt_sem = [fstruct.filt_sem]';
    trial_pp.forceCO.force.filt_psem = [fstruct.filt_psem]';
    trial_pp.forceCO.EMG.rect = cat(1,EMGstruct.rect_mean);
    trial_pp.forceCO.EMG.filt = cat(1,EMGstruct.filt_mean);
    trial_pp.forceCO.trial_coh = trial_coh_force;
    
    fstruct = [trial_avg_EMG.force];
    EMGstruct = [trial_avg_EMG.EMG];
    trial_pp.EMGCO.angle = Aparams.targetAnglesEMG';
    trial_pp.EMGCO.force.mag_mean = [fstruct.filtmag_mean]';
    trial_pp.EMGCO.force.mag_std = [fstruct.filtmag_std]';
    trial_pp.EMGCO.force.mag_pstd = [fstruct.filtmag_pstd]';
    trial_pp.EMGCO.force.mag_sem = [fstruct.filtmag_sem]';
    trial_pp.EMGCO.force.mag_psem = [fstruct.filtmag_psem]';
    trial_pp.EMGCO.force.filt_mean = [fstruct.filt_mean]';
    trial_pp.EMGCO.force.filt_std = [fstruct.filt_std]';
    trial_pp.EMGCO.force.filt_pstd = [fstruct.filt_pstd]';
    trial_pp.EMGCO.force.filt_sem = [fstruct.filt_sem]';
    trial_pp.EMGCO.force.filt_psem = [fstruct.filt_psem]';
    trial_pp.EMGCO.EMG.rect = cat(1,EMGstruct.rect_mean);
    trial_pp.EMGCO.EMG.filt = cat(1,EMGstruct.filt_mean);
    trial_pp.EMGCO.trial_coh = trial_coh_EMG;
    
    save([filepathpp,date,'_s',subject,'_TDP','.mat'],'trial_pp');
end

% %% Calibration
% calibtype = 'EMGCO';
% 
% filenameforce =  [date,'_s',subject,'_',calibtype,'_Force_calib.mat'];
% filenameEMG = [date,'_s',subject,'_',calibtype,'_EMG_calib.mat'];
% filenameparams = [date,'_s',subject,'_params_','001','.mat'];
% 
% load([filepath,filenameforce]);
% load([filepath,filenameEMG]);
% load([filepath,paramfolder,filenameparams]);
% 
% forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};
% 
% % Parameters for analysis
% Aparams.fs = min(forceparams.scanRate,EMGparams.sampleRateEMG);
% Aparams.downsamp = forceparams.scanRate/EMGparams.sampleRateEMG;
% Aparams.channelNameEMG = EMGparams.channelName;
% Aparams.targetAngles = sort(EMGparams.channelAngleCal);
% Aparams.fclF = forceparams.fclF;
% Aparams.fchEMG = 20;
% Aparams.fclEMG = 500;
% Aparams.avgWindow = 200;
% Aparams.EMGOffset = EMGOffset;
% 
% % Create trial data struct and remove failed or incomplete trials
% trial_data = trialCO(forceEMGData,Aparams);
% trial_data = removeFailTrials(trial_data(5:end));
% 
% % Process EMG and force data and add in struct
% trial_data = procEMG(trial_data,Aparams);
% trial_data_EMG_calib = procForce(trial_data,Aparams);
% 
% % Actual target angles (should be the same as taskparams.targetAnglesForce)
% Aparams.targetAnglesForce = sort(unique(extractfield(trial_data_EMG_calib,'angle')));
% 
% % Epoch interval and signals to trial average
% Aparams.epoch = {'ihold',1,'iend',0};
% fields_avg = {'EMG.rect','EMG.avg','force.filt','force.filtmag'};
% 
% % Trial average by angle
% trial_data_avg_calib = trialAngleAvg(trial_data_EMG_calib, Aparams.epoch, fields_avg);
% 
% % Compute mean rectified EMG and filtered force magnitude for each angle.
% % Check calibration values
% EMGmean = zeros(length(trial_data_avg_calib),length(EMGparams.channelSubset)-1);
% forcemean = zeros(length(trial_data_avg_calib),1);
% for i = 1:length(trial_data_avg_calib)
%     EMGmean(i,:) = mean(trial_data_avg_calib(i).EMG.rect,1);
%     forcemean(i) = trial_data_avg_calib(i).force.filtmag_mean;
% end
% EMGScaleCalib = max(EMGmean,[],1);
% EMGmeanCalib = EMGmean;
% 
% % EMG offset values from EMG calibration
% fprintf('\nEMG offset values: \n')
% for k = 1:length(EMGparams.channelSubsetCal)-1
%     fprintf('%s: %1.3f\n',EMGparams.channelName{EMGparams.channelSubset == EMGparams.channelSubsetCal(k)},Aparams.EMGOffset(EMGparams.channelSubset == EMGparams.channelSubsetCal(k)))
% end
% 
% % EMG scaling values from EMG calibration
% fprintf('\nEMG mean values: \n')
% for k = 1:length(EMGparams.channelSubsetCal)-1
%     fprintf('%s: %1.3f\n',EMGparams.channelName{EMGparams.channelSubset == EMGparams.channelSubsetCal(k)},EMGScaleCalib(EMGparams.channelSubset == EMGparams.channelSubsetCal(k)))
% end
% 
% % Mean filtered force magnitude value
% fprintf('\nRecorded mean force value: %1.3f\n\n',round(mean(forcemean)))
