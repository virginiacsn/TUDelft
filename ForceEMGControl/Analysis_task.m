%% Data analysis for individual subject
% close all
clear all
addpath(genpath('Tools'));

date =      '20181022';
subject =   '09';

switch computer
    case 'PCWIN'
        filepath =  ['D:\Student_experiments\Virginia\Data\',date,'\s',subject,'\'];
        paramfolder = 'Parameters\';
    case 'PCWIN64'
        filepath =  ['D:\Student_experiments\Virginia\Data\',date,'\s',subject,'\'];
        paramfolder = 'Parameters\';
    case 'MACI64'
        filepath =  ['/Users/virginia/Documents/MATLAB/Thesis/Data/',date,'/s',subject,'/'];
        paramfolder = 'Parameters/';
end

%% DATA LOADING AND PREPROC
%% Calibration
calibtype = 'EMGCO';

filenameforce =  [date,'_s',subject,'_',calibtype,'_Force_calib.mat'];
filenameEMG = [date,'_s',subject,'_',calibtype,'_EMG_calib.mat'];
filenameparams = [date,'_s',subject,'_params_','001','.mat'];

load([filepath,filenameforce]);
load([filepath,filenameEMG]);
load([filepath,paramfolder,filenameparams]);

forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};

% Parameters for analysis
Aparams.fs = min(forceparams.scanRate,EMGparams.sampleRateEMG);
Aparams.downsamp = forceparams.scanRate/EMGparams.sampleRateEMG;
Aparams.channelNameEMG = EMGparams.channelName;
Aparams.targetAngles = sort(EMGparams.channelAngleCal);
Aparams.fclF = forceparams.fclF;
Aparams.fchEMG = 20;
Aparams.fclEMG = 500;
Aparams.avgWindow = 200;
Aparams.EMGOffset = EMGOffset;

% Create trial data struct and remove failed or incomplete trials
trial_data = trialCO(forceEMGData,Aparams);
trial_data = removeFailTrials(trial_data(5:end));

% Process EMG and force data and add in struct
trial_data = procEMG(trial_data,Aparams);
trial_data_EMG_calib = procForce(trial_data,Aparams);

% Actual target angles (should be the same as taskparams.targetAnglesForce)
Aparams.targetAnglesForce = sort(unique(extractfield(trial_data_EMG_calib,'angle')));

% Epoch interval and signals to trial average
Aparams.epoch = {'ihold',1,'iend',0};
fields_avg = {'EMG.rect','EMG.avg','force.filt','force.filtmag'};

% Trial average by angle
trial_data_avg_calib = trialAngleAvg(trial_data_EMG_calib, Aparams.epoch, fields_avg);

% Compute mean rectified EMG and filtered force magnitude for each angle.
% Check calibration values
EMGmean = zeros(length(trial_data_avg_calib),length(EMGparams.channelSubset)-1);
forcemean = zeros(length(trial_data_avg_calib),1);
for i = 1:length(trial_data_avg_calib)
    EMGmean(i,:) = mean(trial_data_avg_calib(i).EMG.rect,1);
    forcemean(i) = trial_data_avg_calib(i).force.filtmag_mean;
end
EMGScaleCalib = max(EMGmean,[],1);
EMGmeanCalib = EMGmean;

% EMG offset values from EMG calibration
fprintf('\nEMG offset values: \n')
for k = 1:length(EMGparams.channelSubsetCal)-1
    fprintf('%s: %1.3f\n',EMGparams.channelName{EMGparams.channelSubset == EMGparams.channelSubsetCal(k)},Aparams.EMGOffset(EMGparams.channelSubset == EMGparams.channelSubsetCal(k)))
end

% EMG scaling values from EMG calibration
fprintf('\nEMG mean values: \n')
for k = 1:length(EMGparams.channelSubsetCal)-1
    fprintf('%s: %1.3f\n',EMGparams.channelName{EMGparams.channelSubset == EMGparams.channelSubsetCal(k)},EMGScaleCalib(EMGparams.channelSubset == EMGparams.channelSubsetCal(k)))
end

% Mean filtered force magnitude value
fprintf('\nRecorded mean force value: %1.3f\n\n',round(mean(forcemean)))

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
Aparams.fchEMG = 20;
Aparams.fclEMG = 500;
Aparams.fnEMG = 50;
Aparams.avgWindow = 200;

% Epoch interval {epoch start, time start, epoch end, time end} and fields to trial average
Aparams.epoch = {'ihold',1,'iend',0};
fields_avg = {'EMG.raw','EMG.filt','EMG.rect','EMG.avg','force.filt','force.filtmag'};
% Fields to trial append
fields_app = {'EMG.raw','EMG.filt','EMG.rect'};

% Trial data struct for force-control task, will append trial data for each
% block (code)
trial_data_force = [];

for i = 1:length(codeF)
    
    code = codeF{i};
    
    filenameforce =  [date,'_s',subject,'_',task,'_Force_',code,'.mat'];
    filenameEMG = [date,'_s',subject,'_',task,'_EMG_',code,'.mat'];
    
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

% Final force-control target angles based on successful trials
Aparams.targetAnglesForce = sort(unique(extractfield(trial_data_force,'angle')));

% Trial average by angle
trial_data_avg_force = trialAngleAvg(trial_data_force, Aparams.epoch, fields_avg);
% FFT of averaged EMG data
trial_data_avg_force = fftEMG(trial_data_avg_force);

% Select number of trials to append by finding minimum number of trials per
% target
mintrialsForce = min([trial_data_avg_force.ntrials]);

% Trial append by angle. Limit to 20 trials (100 sec of data)
trial_data_app_force = trialAngleApp(trial_data_force, Aparams.epoch, fields_app,[],20);

% Compute mean rectified EMG and filtered force magnitude for each angle.
% Check calibration values
EMGmean = zeros(length(trial_data_avg_force),length(EMGparams.channelSubset)-1);
forcemean = zeros(length(trial_data_avg_force),1);
for i = 1:length(trial_data_avg_force)
    EMGmean(i,:) = mean(trial_data_avg_force(i).EMG.rect,1);
    forcemean(i) = trial_data_avg_force(i).force.filtmag_mean;
end
EMGmeanForce = EMGmean;

% EMG scaling through max EMG for force-control
EMGScaleForce = max(EMGmean,[],1)';
% EMG target tolerance values
EMGTolForce = min(EMGmean,[],1)'./EMGScaleForce;

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
trial_data_avg_EMG = trialAngleAvg(trial_data_EMG, Aparams.epoch, fields_avg);
% FFT of averaged EMG data
trial_data_avg_EMG = fftEMG(trial_data_avg_EMG);

% Select number of trials to append by finding minimum number of trials per
% target
mintrialsEMG = min([trial_data_avg_EMG.ntrials]);

% Trial append by angle. Limit to 25 trials (100 sec of data)
trial_data_app_EMG = trialAngleApp(trial_data_EMG, Aparams.epoch, fields_app,[],25);

% Compute mean rectified EMG and filtered force magnitude for each angle.
% Check calibration values
EMGmean = zeros(length(trial_data_avg_EMG),length(EMGparams.channelSubset)-1);
forcemean = zeros(length(trial_data_avg_EMG),1);
for i = 1:length(trial_data_avg_EMG)
    EMGmean(i,:) = mean(trial_data_avg_EMG(i).EMG.rect,1);
    forcemean(i) = trial_data_avg_EMG(i).force.filtmag_mean;
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
Aparams.angComp = {};
for i = 1:length(angComp)
    Aparams.angComp{i} = angComp(i)*[1 1];
end

% Muscles names to compare
Aparams.muscComp = channelNames(ismember(muscAngles,Aparams.targetAnglesEMG));

%% LIMITS for force and EMG plots
Flim = round(max(forcemean))/taskparams.targetForce; EMGlim = 60;

% EMG_lim(i,j) will be max EMG level between force-control and EMG-control
% where rows are target angles to compare and cols muscles (all)
EMG_lim = [];
for i = 1:length(Aparams.angComp)
    for j = 1:length(EMGparams.channelName)-1
        EMG_lim(i,j) = max(max(trial_data_avg_force(Aparams.targetAnglesForce == Aparams.angComp{i}(1)).EMG.rect(:,j)),...
            max(trial_data_avg_EMG(Aparams.targetAnglesEMG == Aparams.angComp{i}(2)).EMG.rect(:,j)));
    end
end

% FFT lim
fc = 100;

% Limits for fft of rectified EMG
fft_fields = {'raw','filt','rect'};
for h = 1:length(fft_fields)
    for i = 1:length(Aparams.angComp)
        for j = 1:length(EMGparams.channelName)-1
            fft_lim.(fft_fields{h})(i,j) = max([max(abs(trial_data_avg_force...
                (Aparams.targetAnglesForce == Aparams.angComp{i}(1)).EMG.(['fft',fft_fields{h}])(3:end-2,j))),...
                max(abs(trial_data_avg_EMG(Aparams.targetAnglesEMG == Aparams.angComp{i}(2)).EMG.(['fft',fft_fields{h}])(3:end-2,j)))]);
        end
    end
end

%% FORCE FIGURES (low-pass filtered)
%% Force in time. Fig per task. Subplot per target angle.
figure('Name','ForceCO; Force in time');
for i = 1:length(Aparams.targetAnglesForce)
    if rem(length(Aparams.targetAnglesForce),2) == 0
        subplot(2,length(Aparams.targetAnglesForce)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesForce),i);
    end
    for j = 1:2
        plot(trial_data_avg_force(i).ts,trial_data_avg_force(i).force.filt(:,j));
        hold on;
    end
    ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Time [s]'); ylabel('Force [N]');
    title(['Target: ',num2str(rad2deg(trial_data_avg_force(i).angle)),' deg']);
    legend('Fx','Fy')
end

figure('Name','EMGCO; Force in time');
for i = 1:length(Aparams.targetAnglesEMG)
    if rem(length(Aparams.targetAnglesEMG),2) == 0
        subplot(2,length(Aparams.targetAnglesEMG)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesEMG),i);
    end
    for j = 1:2
        plot(trial_data_avg_EMG(i).ts,trial_data_avg_EMG(i).force.filt(:,j));
        hold on;
    end
    ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Time [s]'); ylabel('Force [N]');
    title(['Target: ',num2str(rad2deg(trial_data_avg_EMG(i).angle)),' deg']);
    legend('Fx','Fy')
end

%% Force trajectory. Fig per task. Subplot per target angle.
figure('Name','ForceCO; Force trajectory');
for i = 1:length(Aparams.targetAnglesForce)
    if rem(length(Aparams.targetAnglesForce),2) == 0
        subplot(2,length(Aparams.targetAnglesForce)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesForce),i);
    end
    plot(trial_data_avg_force(i).force.filt(:,1),trial_data_avg_force(i).force.filt(:,2));
    hold on;
    h1 = plot(trial_data_avg_force(i).force.filt(1,1),trial_data_avg_force(i).force.filt(1,2),'go');
    h2 = plot(trial_data_avg_force(i).force.filt(end,1),trial_data_avg_force(i).force.filt(end,2),'ro');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]);ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    title(['Target: ',num2str(rad2deg(trial_data_avg_force(i).angle)),' deg']);
    legend([h1,h2],'Start','End');
end

figure('Name','EMGCO; Force trajectory');
for i = 1:length(Aparams.targetAnglesEMG)
    if rem(length(Aparams.targetAnglesEMG),2) == 0
        subplot(2,length(Aparams.targetAnglesEMG)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesEMG),i);
    end
    plot(trial_data_avg_EMG(i).force.filt(:,1),trial_data_avg_EMG(i).force.filt(:,2));
    hold on;
    h1 = plot(trial_data_avg_EMG(i).force.filt(1,1),trial_data_avg_EMG(i).force.filt(1,2),'go');
    h2 = plot(trial_data_avg_EMG(i).force.filt(end,1),trial_data_avg_EMG(i).force.filt(end,2),'ro');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]);ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    title(['Target: ',num2str(rad2deg(trial_data_avg_EMG(i).angle)),' deg']);
    legend([h1,h2],'Start','End');
end

%% Force trajectory. Subplot per task. Target angles in same subplot.
figure('Name','Force trajectory');
subplot(1,2,1)
for i = 1:length(Aparams.targetAnglesForce)
    plot(trial_data_avg_force(i).force.filt(:,1),trial_data_avg_force(i).force.filt(:,2));
    hold on;
    h1 = plot(trial_data_avg_force(i).force.filt(1,1),trial_data_avg_force(i).force.filt(1,2),'go');
    h2 = plot(trial_data_avg_force(i).force.filt(end,1),trial_data_avg_force(i).force.filt(end,2),'ro');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]);ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    text(trial_data_avg_force(i).force.filt(end,1),trial_data_avg_force(i).force.filt(end,2),num2str(rad2deg(Aparams.targetAnglesForce(i))),'VerticalAlignment','top','HorizontalAlignment','left');
    legend([h1 h2],'Start','End');
end
title('ForceCO');

subplot(1,2,2)
for i = 1:length(Aparams.targetAnglesEMG)
    plot(trial_data_avg_EMG(i).force.filt(:,1),trial_data_avg_EMG(i).force.filt(:,2));
    hold on;
    h1 = plot(trial_data_avg_EMG(i).force.filt(1,1),trial_data_avg_EMG(i).force.filt(1,2),'go');
    h2 = plot(trial_data_avg_EMG(i).force.filt(end,1),trial_data_avg_EMG(i).force.filt(end,2),'ro');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-Flim Flim]);ylim(taskparams.targetForce*[-Flim Flim]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    text(trial_data_avg_EMG(i).force.filt(end,1),trial_data_avg_EMG(i).force.filt(end,2),num2str(rad2deg(Aparams.targetAnglesEMG(i))),'VerticalAlignment','top','HorizontalAlignment','left');
    legend([h1 h2],'Start','End');
end
title('EMGCO');

%% EMG FIGURES
%% TIME analysis
%% Fig per muscle (all) and task. Subplot per target. 
% Compare targets for each task and muscle.
for j = 1:length(EMGparams.channelName)-1
    figure('Name',['ForceCO; EMG ',EMGparams.channelName{j}]);
    for i = 1:length(Aparams.targetAnglesForce)
        if rem(length(Aparams.targetAnglesForce) ,2) == 0
            subplot(2,length(Aparams.targetAnglesForce) /2,i);
        else
            subplot(1,length(Aparams.targetAnglesForce) ,i);
        end
        plot(trial_data_avg_force(i).ts,trial_data_avg_force(i).EMG.rect(:,j));
        hold on;
        plot(trial_data_avg_force(i).ts,trial_data_avg_force(i).EMG.avg(:,j));
        ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_force(i).EMG.rect(:))+50]);
        xlim([0 trial_data_avg_force(i).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['Target: ',num2str(rad2deg(trial_data_avg_force(i).angle)),' deg']);
    end
end

for j = 1:length(EMGparams.channelName)-1
    figure('Name',['EMGCO; EMG ',EMGparams.channelName{j}]);
    for i = 1:length(Aparams.targetAnglesEMG)
        if rem(length(Aparams.targetAnglesEMG),2) == 0
            subplot(2,length(Aparams.targetAnglesEMG)/2,i);
        else
            subplot(1,length(Aparams.targetAnglesEMG),i);
        end
        plot(trial_data_avg_EMG(i).ts,trial_data_avg_EMG(i).EMG.rect(:,j));
        hold on;
        plot(trial_data_avg_EMG(i).ts,trial_data_avg_EMG(i).EMG.avg(:,j));
        ylim([0 EMGlim]); %ylim([0 EMG_lim(i,j)+10]);
        xlim([0 trial_data_avg_EMG(i).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['Target: ',num2str(rad2deg(trial_data_avg_EMG(i).angle)),' deg']);
    end
end

%% Fig per task and target. Subplot per muscle (all). 
% Compare all muscles for each task and target.
for j = 1:length(Aparams.targetAnglesForce)
    figure('Name',['ForceCO; Target: ',num2str(rad2deg(Aparams.targetAnglesForce(j))), ' deg']);
    for i = 1:length(EMGparams.channelName)-1
        if rem(length(EMGparams.channelName)-1,2) == 0
            subplot(2,(length(EMGparams.channelName)-1)/2,i);
        else
            subplot(1,length(EMGparams.channelName)-1,i);
        end
        plot(trial_data_avg_force(j).ts,trial_data_avg_force(j).EMG.rect(:,i));
        hold on;
        plot(trial_data_avg_force(j).ts,trial_data_avg_force(j).EMG.avg(:,i));
        ylim([0 EMG_lim(j,i)+5]);%ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_force(j).EMG.rect(:))+50]);
        xlim([0 trial_data_avg_force(j).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['Musc: ',EMGparams.channelName{i}]);
    end
end

for j = 1:length(Aparams.targetAnglesEMG)
    figure('Name',['EMGCO; Target: ',num2str(rad2deg(Aparams.targetAnglesEMG(j))), ' deg']);
    for i = 1:length(EMGparams.channelName)-1
        if rem(length(EMGparams.channelName)-1,2) == 0
            subplot(2,(length(EMGparams.channelName)-1)/2,i);
        else
            subplot(1,length(EMGparams.channelName)-1,i);
        end
        plot(trial_data_avg_EMG(j).ts,trial_data_avg_EMG(j).EMG.rect(:,i));
        hold on;
        plot(trial_data_avg_EMG(j).ts,trial_data_avg_EMG(j).EMG.avg(:,i));
        ylim([0 EMG_lim(j,i)+5]);%ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_EMG(j).EMG.rect(:))+50]);
        xlim([0 trial_data_avg_EMG(j).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['Musc: ',EMGparams.channelName{i}]);
    end
end

%% Fig per task. Subplot per muscle (all) (col) and target (row). 
% Compare all muscles and targets for each task.
h = 0;
figure('Name','ForceCO');
for j = 1:length(Aparams.targetAnglesForce)
    for i = 1:length(EMGparams.channelName)-1
        h = h+1;
        subplot(length(Aparams.targetAnglesForce),length(EMGparams.channelName)-1,h);
        plot(trial_data_avg_force(j).ts,trial_data_avg_force(j).EMG.rect(:,i));
        hold on;
        plot(trial_data_avg_force(j).ts,trial_data_avg_force(j).EMG.avg(:,i));
        ylim([0 max(EMG_lim(:,i))]);%ylim([0 EMGlim])%ylim([0 max(trial_data_avg_force(j).EMG.rect(:))+10]);%
        xlim([0 trial_data_avg_force(j).ts(end)]);
        if j == length(Aparams.targetAnglesForce)
            xlabel('Time [s]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.targetAnglesForce(j))) ,' deg']);
        end
        if j == 1
            title(['Musc: ',EMGparams.channelName{i}]);
        end
    end
end

h = 0;
figure('Name','EMGCO');
for j = 1:length(Aparams.targetAnglesEMG)
    for i = 1:length(EMGparams.channelName)-1
        h = h+1;
        subplot(length(Aparams.targetAnglesEMG),length(EMGparams.channelName)-1,h);
        
        plot(trial_data_avg_EMG(j).ts,trial_data_avg_EMG(j).EMG.rect(:,i));
        hold on;
        plot(trial_data_avg_EMG(j).ts,trial_data_avg_EMG(j).EMG.avg(:,i));
        ylim([0 max(EMG_lim(:,i))]);%ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_EMG(j).EMG.rect(:))+50]);
        xlim([0 trial_data_avg_EMG(j).ts(end)]);
        if j == length(Aparams.targetAnglesEMG)
            xlabel('Time [s]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.targetAnglesEMG(j))) ,' deg']);
        end
        if j == 1
            title(['Musc: ',EMGparams.channelName{i}]);
        end
    end
end

%% Fig per task. Subplot per muscle (control) (col) and target (row). 
% Compare control muscles and targets for each task.
h = 0;
figure('Name','ForceCO');
for j = 1:length(Aparams.targetAnglesForce)
    for i = 1:length(Aparams.chanControl)
        h = h+1;
        subplot(length(Aparams.targetAnglesForce),length(Aparams.chanControl),h);
        plot(trial_data_avg_force(j).ts,trial_data_avg_force(j).EMG.rect(:,Aparams.chanControl(i)));
        hold on;
        plot(trial_data_avg_force(j).ts,trial_data_avg_force(j).EMG.avg(:,Aparams.chanControl(i)));
        ylim([0 EMGlim]); %ylim([0 max(trial_data_avg_force(j).EMG.rect(:))+10]);%
        xlim([0 trial_data_avg_force(j).ts(end)]);
        if j == length(Aparams.targetAnglesForce)
            xlabel('Time [s]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.targetAnglesForce(j))) ,' deg']);
        end
        if j == 1
            title(['Musc: ',Aparams.chanControlName{i}]);
        end
    end
end

h = 0;
figure('Name','EMGCO');
for j = 1:length(Aparams.targetAnglesEMG)
    for i = 1:length(Aparams.chanControl)
        h = h+1;
        subplot(length(Aparams.targetAnglesEMG),length(Aparams.chanControl),h);
        
        plot(trial_data_avg_EMG(j).ts,trial_data_avg_EMG(j).EMG.rect(:,Aparams.chanControl(i)));
        hold on;
        plot(trial_data_avg_EMG(j).ts,trial_data_avg_EMG(j).EMG.avg(:,Aparams.chanControl(i)));
        ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_EMG(j).EMG.rect(:))+50]);
        xlim([0 trial_data_avg_EMG(j).ts(end)]);
        if j == length(Aparams.targetAnglesEMG)
            xlabel('Time [s]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.targetAnglesEMG(j))) ,' deg']);
        end
        if j == 1
            title(['Musc: ',Aparams.chanControlName{i}]);
        end
    end
end

%% Fig per muscle (all). Subplot per task (row) and target from angComp (col). 
% Compare task and target for each muscle (all).
for i = 1:length(EMGparams.channelName)-1
    figure('Name',EMGparams.channelName{i});
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);

    h = 1;
    for j = 1:length(Aparams.angComp)
        subplot(length(Aparams.angComp),2,h)
        iangf = find([trial_data_avg_force.angle] == Aparams.angComp{j}(1));
        plot(trial_data_avg_force(iangf).ts,trial_data_avg_force(iangf).EMG.rect(:,i))
        hold on
        plot(trial_data_avg_force(iangf).ts,trial_data_avg_force(iangf).EMG.avg(:,i),'r')
        xlim([0 trial_data_avg_force(iangf).ts(end)]);
        if EMG_lim(j,i)>0
            ylim([0 EMG_lim(j,i)]);
        else
            ylim([0 EMGlim])
        end
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['ForceCO; Target: ',num2str(rad2deg(Aparams.angComp{j}(1))),' deg']);
        
        subplot(length(Aparams.angComp),2,h+1)
        iangE = find([trial_data_avg_EMG.angle] == Aparams.angComp{j}(2));
        plot(trial_data_avg_EMG(iangE).ts,trial_data_avg_EMG(iangE).EMG.rect(:,i))
        hold on
        plot(trial_data_avg_EMG(iangE).ts,trial_data_avg_EMG(iangE).EMG.avg(:,i),'r')
        xlim([0 trial_data_avg_EMG(iangE).ts(end)]);
        if EMG_lim(j,i)>0
            ylim([0 EMG_lim(j,i)]);
        else
            ylim([0 EMGlim])
        end
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['EMGCO; Target: ',num2str(rad2deg(Aparams.angComp{j}(2))),' deg (',Aparams.muscComp{j},')']);
        
        h = h+2;
    end
end

%% Fig per target from angComp. Subplot per task (row) per muscle (all) (col). 
% Compare muscles (all) and tasks for each target from angComp.
for j = 1:length(Aparams.angComp)
    figure('Name',['Targets: Force-',num2str(rad2deg(Aparams.angComp{j}(1))),'; EMG-',num2str(rad2deg(Aparams.angComp{j}(2))),' (',Aparams.muscComp{j},')'])
    
    for i = 1:length(EMGparams.channelName)-1
        subplot(2,length(EMGparams.channelName)-1,i)
        iangf = find([trial_data_avg_force.angle] == Aparams.angComp{j}(1));
        plot(trial_data_avg_force(iangf).ts,trial_data_avg_force(iangf).EMG.rect(:,i))
        hold on
        plot(trial_data_avg_force(iangf).ts,trial_data_avg_force(iangf).EMG.avg(:,i),'r')
        ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_force(iangf).EMG.rect(:))+10]);
        xlim([0 trial_data_avg_force(iangf).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title([EMGparams.channelName{i}]);
        
        subplot(2,length(EMGparams.channelName)-1,i+length(EMGparams.channelName)-1)
        iangE = find([trial_data_avg_EMG.angle] == Aparams.angComp{j}(2));
        plot(trial_data_avg_EMG(iangE).ts,trial_data_avg_EMG(iangE).EMG.rect(:,i))
        hold on
        plot(trial_data_avg_EMG(iangE).ts,trial_data_avg_EMG(iangE).EMG.avg(:,i),'r')
        ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_EMG(iangE).EMG.rect(:))+10]);
        xlim([0 trial_data_avg_EMG(iangE).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title([EMGparams.channelName{i}]);
    end
end

%% EMG MEAN/VAR/SNR
%% Fig of EMG mean of all muscles. Subplot per target from angComp. 
% Compare means of all muscles for each target.
figure('Name','EMG mean')
for j = 1:length(Aparams.angComp)
    for i = 1:length(EMGparams.channelName)-1
        subplot(length(Aparams.angComp),1,j)
        
        iangf = find([trial_data_avg_force.angle] == Aparams.angComp{j}(1));
        stem(i,mean(trial_data_avg_force(iangf).EMG.rect(:,i)),'b')
        hold on
        iangE = find([trial_data_avg_EMG.angle] == Aparams.angComp{j}(2));
        stem(i,mean(trial_data_avg_EMG(iangE).EMG.rect(:,i)),'r')
        
        xlim([0 length(EMGparams.channelName)]);
        ylabel('Mean EMG [-]');
        title(['Targets: Force-',num2str(rad2deg(Aparams.angComp{j}(1))),'; EMG-',num2str(rad2deg(Aparams.angComp{j}(2))),' (',Aparams.muscComp{j},')']);
        legend('ForceCO','EMGCO')
    end
    xticklabels([{''},EMGparams.channelName(1:end-1)])
end

%% Fig of EMG mean of control muscles. Subplot per target from angComp. 
% Compare means of control muscles for each target.
figure('Name','EMG mean')
for j = 1:length(Aparams.angComp)
    for i = 1:length(Aparams.chanControl)
        subplot(length(Aparams.angComp),1,j)
        
        iangf = find([trial_data_avg_force.angle] == Aparams.angComp{j}(1));
        stem(i,mean(trial_data_avg_force(iangf).EMG.rect(:,Aparams.chanControl(i))),'b')
        hold on
        iangE = find([trial_data_avg_EMG.angle] == Aparams.angComp{j}(2));
        stem(i,mean(trial_data_avg_EMG(iangE).EMG.rect(:,Aparams.chanControl(i))),'r')
        
        xlim([0 length(Aparams.chanControl)+1]);
        ylabel('Mean EMG [-]');
        title(['Targets: Force-',num2str(rad2deg(Aparams.angComp{j}(1))),'; EMG-',num2str(rad2deg(Aparams.angComp{j}(2))),' (',Aparams.muscComp{j},')']);
        legend('ForceCO','EMGCO')
    end
    xticks([0:length(Aparams.chanControl)]);
    xticklabels([{''},Aparams.chanControlName]);  
end

%% Fig of EMG variance of all muscles. Subplot per target from angComp.
figure('Name','EMG var')
for j = 1:length(Aparams.angComp)
    for i = 1:length(EMGparams.channelName)-1
        subplot(length(Aparams.angComp),1,j)
        
        iangf = find([trial_data_avg_force.angle] == Aparams.angComp{j}(1));
        stem(i,var(trial_data_avg_force(iangf).EMG.rect(:,i)),'b')
        hold on
        iangE = find([trial_data_avg_EMG.angle] == Aparams.angComp{j}(2));
        stem(i,var(trial_data_avg_EMG(iangE).EMG.rect(:,i)),'r')
        
        xlim([0 length(EMGparams.channelName)]);
        ylabel('var [-]');
        title(['Targets: Force-',num2str(rad2deg(Aparams.angComp{j}(1))),'; EMG-',num2str(rad2deg(Aparams.angComp{j}(2))),' (',Aparams.muscComp{j},')']);
        legend('ForceCO','EMGCO')
    end
    xticklabels([{''},EMGparams.channelName(1:end-1)])
end

%% Fig of EMG SNR of all muscles. Subplot per target from angComp.
figure('Name','EMG SNR')
for j = 1:length(Aparams.angComp)
    for i = 1:length(EMGparams.channelName)-1
        subplot(length(Aparams.angComp),1,j)
        
        iangf = find([trial_data_avg_force.angle] == Aparams.angComp{j}(1));
        stem(i,mean(trial_data_avg_force(iangf).EMG.rect(:,i))./std(trial_data_avg_force(iangf).EMG.rect(:,i)),'b')
        hold on
        iangE = find([trial_data_avg_EMG.angle] == Aparams.angComp{j}(2));
        stem(i,mean(trial_data_avg_EMG(iangE).EMG.rect(:,i))./std(trial_data_avg_EMG(iangE).EMG.rect(:,i)),'r')
        
        xlim([0 length(EMGparams.channelName)]);
        ylabel('SNR [-]');
        title(['Targets: Force-',num2str(rad2deg(Aparams.angComp{j}(1))),'; EMG-',num2str(rad2deg(Aparams.angComp{j}(2))),' (',Aparams.muscComp{j},')']);
        legend('ForceCO','EMGCO')
    end
    xticklabels([{''},EMGparams.channelName(1:end-1)])
end

%% FREQUENCY analysis
%% FFT
%% Fig per muscle (all). Subplot per task (row) and target from angComp (col).
% Compare targets from angComp between tasks for each muscle (all).
for i = 1:length(EMGparams.channelName)-1
    figure('Name',EMGparams.channelName{i});
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    
    h = 1;
    for j = 1:length(Aparams.angComp)
        subplot(length(Aparams.angComp),2,h)
        iangf = find([trial_data_avg_force.angle] == Aparams.angComp{j}(1));
        plot(trial_data_avg_force(iangf).fv,abs(trial_data_avg_force(iangf).EMG.fftraw(:,i)));
        hold on;
        plot(trial_data_avg_force(iangf).fv,abs(trial_data_avg_force(iangf).EMG.fftfilt(:,i)));
        plot(trial_data_avg_force(iangf).fv,abs(trial_data_avg_force(iangf).EMG.fftrect(:,i)));
        xlim([1 fc]); ylim([0 max([fft_lim.raw(j,i),fft_lim.filt(j,i),fft_lim.rect(j,i)])]);
        if j == length(Aparams.angComp)
            xlabel('Frequency [Hz]');
        end
        ylabel('a.u. [-]');
        legend('Raw','Filt','Rect')
        title(['ForceCO; Target: ',num2str(rad2deg(Aparams.angComp{j}(1))),' deg']);
        
        subplot(length(Aparams.angComp),2,h+1)
        iangE = find([trial_data_avg_EMG.angle] == Aparams.angComp{j}(2));
        plot(trial_data_avg_EMG(iangE).fv,abs(trial_data_avg_EMG(iangE).EMG.fftraw(:,i)));
        hold on;
        plot(trial_data_avg_EMG(iangE).fv,abs(trial_data_avg_EMG(iangE).EMG.fftfilt(:,i)));
        plot(trial_data_avg_EMG(iangE).fv,abs(trial_data_avg_EMG(iangE).EMG.fftrect(:,i)));
        xlim([1 fc]); ylim([0 max([fft_lim.raw(j,i),fft_lim.filt(j,i),fft_lim.rect(j,i)])]);
        if j == length(Aparams.angComp)
            xlabel('Frequency [Hz]');
        end
        ylabel('a.u. [-]');
        legend('Raw','Filt','Rect')
        title(['EMGCO; Target: ',num2str(rad2deg(Aparams.angComp{j}(2))),' deg (',Aparams.muscComp{j},')']);
        
        h = h+2;
    end
end

%% Fig per task. Subplot per muscle (control (col) and target from angComp (row).
% Compare target and muscle (control) for each task.
h = 0;
figure('Name','ForceCO');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for j = 1:length(Aparams.targetAnglesForce)
    for i = 1:length(Aparams.chanControl)
        h = h+1;
        subplot(length(Aparams.targetAnglesForce),length(Aparams.chanControl),h);
        plot(trial_data_avg_force(j).fv,abs(trial_data_avg_force(j).EMG.fftraw(:,Aparams.chanControl(i))));
        hold on;
        plot(trial_data_avg_force(j).fv,abs(trial_data_avg_force(j).EMG.fftfilt(:,Aparams.chanControl(i))));
        plot(trial_data_avg_force(j).fv,abs(trial_data_avg_force(j).EMG.fftrect(:,Aparams.chanControl(i))));
        xlim([1 fc]); ylim([0 max([max(fft_lim.raw(:,Aparams.chanControl(i))),max(fft_lim.filt(:,Aparams.chanControl(i))),max(fft_lim.rect(:,Aparams.chanControl(i)))])]);
        legend('Raw','Filt','Rect');
        if j == length(Aparams.targetAnglesForce)
            xlabel('Frequency [Hz]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.targetAnglesForce(j))) ,' deg']);
        end
        if j == 1
            title(['Musc: ',Aparams.chanControlName{i}]);
        end
    end
end

h = 0;
figure('Name','EMGCO');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for j = 1:length(Aparams.targetAnglesEMG)
    for i = 1:length(Aparams.chanControl)
        h = h+1;
        subplot(length(Aparams.targetAnglesEMG),length(Aparams.chanControl),h);
        plot(trial_data_avg_EMG(j).fv,abs(trial_data_avg_EMG(j).EMG.fftraw(:,Aparams.chanControl(i))));
        hold on;
        plot(trial_data_avg_EMG(j).fv,abs(trial_data_avg_EMG(j).EMG.fftfilt(:,Aparams.chanControl(i))));
        plot(trial_data_avg_EMG(j).fv,abs(trial_data_avg_EMG(j).EMG.fftrect(:,Aparams.chanControl(i))));
        xlim([1 fc]); ylim([0 max([max(fft_lim.raw(:,Aparams.chanControl(i))),max(fft_lim.filt(:,Aparams.chanControl(i))),max(fft_lim.rect(:,Aparams.chanControl(i)))])]);
        legend('Raw','Filt','Rect');
        if j == length(Aparams.targetAnglesEMG)
            xlabel('Frequency [Hz]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.targetAnglesEMG(j))) ,' deg']);
        end
        if j == 1
            title(['Musc: ',Aparams.chanControlName{i}]);
        end
    end
end
