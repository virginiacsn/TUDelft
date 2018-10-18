%% Data analysis for individual subject
% close all
clear all
addpath(genpath('Tools'));

date =      '20181017';
subject =   '07';

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

%% Limits for force and EMG plots
Flim = round(max(forcemean))/taskparams.targetForce; EMGlim = 60;

% EMG_lim(i,j) will be max EMG level between force-control and EMG-control
% where rows are target angles to compare and cols muscles (all)
EMG_lim = [];
for i = 1:length(Aparams.angComp)
    for j = 1:length(EMGparams.channelName)-1
        EMG_lim(i,j) = max(max(trial_data_avg_force(Aparams.targetAnglesForce == Aparams.angComp{i}(1)).EMG.rect(:,j)),max(trial_data_avg_EMG(Aparams.targetAnglesEMG == Aparams.angComp{i}(2)).EMG.rect(:,j)));
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
        plot(trial_data_avg_force(iangf).fv,abs(fft(trial_data_avg_force(iangf).EMG.raw(:,i)))/length(trial_data_avg_force(iangf).fv));
        hold on;
        plot(trial_data_avg_force(iangf).fv,abs(fft(trial_data_avg_force(iangf).EMG.filt(:,i)))/length(trial_data_avg_force(iangf).fv));
        plot(trial_data_avg_force(iangf).fv,abs(fft(trial_data_avg_force(iangf).EMG.rect(:,i)))/length(trial_data_avg_force(iangf).fv));
        xlim([1 fc]);  ylim([0 0.6]);
        xlabel('Frequency [Hz]'); ylabel('a.u. [-]');
        legend('Raw','Filt','Rect')
        title(['ForceCO; Target: ',num2str(rad2deg(Aparams.angComp{j}(1))),' deg']);
        
        subplot(length(Aparams.angComp),2,h+1)
        iangE = find([trial_data_avg_EMG.angle] == Aparams.angComp{j}(2));
        plot(trial_data_avg_EMG(iangE).fv,abs(fft(trial_data_avg_EMG(iangE).EMG.raw(:,i)))/length(trial_data_avg_EMG(iangE).fv));
        hold on;
        plot(trial_data_avg_EMG(iangE).fv,abs(fft(trial_data_avg_EMG(iangE).EMG.filt(:,i)))/length(trial_data_avg_EMG(iangE).fv));
        plot(trial_data_avg_EMG(iangE).fv,abs(fft(trial_data_avg_EMG(iangE).EMG.rect(:,i)))/length(trial_data_avg_EMG(iangE).fv));
        xlim([1 fc]); ylim([0 0.6]);
        xlabel('Frequency [Hz]'); ylabel('a.u. [-]');
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
        iangf = j;
        subplot(length(Aparams.targetAnglesForce),length(Aparams.chanControl),h);
        plot(trial_data_avg_force(iangf).fv,abs(fft(trial_data_avg_force(iangf).EMG.raw(:,Aparams.chanControl(i))))/length(trial_data_avg_force(iangf).fv));
        hold on;
        plot(trial_data_avg_force(iangf).fv,abs(fft(trial_data_avg_force(iangf).EMG.filt(:,Aparams.chanControl(i))))/length(trial_data_avg_force(iangf).fv));
        plot(trial_data_avg_force(iangf).fv,abs(fft(trial_data_avg_force(iangf).EMG.rect(:,Aparams.chanControl(i))))/length(trial_data_avg_force(iangf).fv));
        xlim([1 fc]); ylim([0 0.8]);
        legend('Raw','Filt','Rect');
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
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for j = 1:length(Aparams.targetAnglesEMG)
    for i = 1:length(Aparams.chanControl)
        h = h+1;
        iangE = j;
        subplot(length(Aparams.targetAnglesEMG),length(Aparams.chanControl),h);
        plot(trial_data_avg_EMG(iangE).fv,abs(fft(trial_data_avg_EMG(iangE).EMG.raw(:,Aparams.chanControl(i))))/length(trial_data_avg_EMG(iangE).fv));
        hold on;
        plot(trial_data_avg_EMG(iangE).fv,abs(fft(trial_data_avg_EMG(iangE).EMG.filt(:,Aparams.chanControl(i))))/length(trial_data_avg_EMG(iangE).fv));
        plot(trial_data_avg_EMG(iangE).fv,abs(fft(trial_data_avg_EMG(iangE).EMG.rect(:,Aparams.chanControl(i))))/length(trial_data_avg_EMG(iangE).fv));
        xlim([1 fc]); ylim([0 0.8]);
        legend('Raw','Filt','Rect');
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

%% COHERENCE
Aparams.cohparams.data = 'app';
Aparams.cohparams.tseg = 1;
Aparams.cohparams.nseg = 10;
Aparams.cohparams.my_nseg = 10;
Aparams.cohparams.window = @(N) hanning(N);
fields_coh = {'raw','filt','rect'};
fc = 200;

if strcmp(Aparams.cohparams.data,'avg')
    trial_data_coh_force = cohStruct(trial_data_avg_force,EMGparams.channelName,fields_coh,Aparams.cohparams);
    trial_data_coh_EMG = cohStruct(trial_data_avg_EMG,EMGparams.channelName,fields_coh,Aparams.cohparams);
else
    trial_data_coh_force = cohStruct(trial_data_app_force,EMGparams.channelName,fields_coh,Aparams.cohparams);
    trial_data_coh_EMG = cohStruct(trial_data_app_EMG,EMGparams.channelName,fields_coh,Aparams.cohparams);
end

% Number of muscles
nmusc = length(EMGparams.channelName)-1;
% Number of muscle pair combinations
nmusccomb = (length(EMGparams.channelName)-1)*(length(EMGparams.channelName)-2)/2; % n*(n-1)/2

%% Field and channels to plot
field = 'filt';

% Find muscle in all muscle pair combinations
plotBB = sum(reshape(contains([trial_data_coh_force(1).(field).muscles{:}],'BB'),[2,nmusccomb]));
plotTLH = sum(reshape(contains([trial_data_coh_force(1).(field).muscles{:}],'TLH'),[2,nmusccomb]));
plotDA = sum(reshape(contains([trial_data_coh_force(1).(field).muscles{:}],'DA'),[2,nmusccomb]));
plotDP = sum(reshape(contains([trial_data_coh_force(1).(field).muscles{:}],'DP'),[2,nmusccomb]));

% Find indexes of muscle pair combinations to plot. Combinations of control
% muscles.
plotmusc = find((plotBB&plotDA)|(plotDP&plotBB)|(plotDA&plotTLH)|(plotDA&plotDP)|(plotBB&plotTLH)|(plotTLH&plotDP));

% Coherence (my_coh) limits for each field in stuct. coh_lim(i,j) will be max coh 
% level between force-control and EMG-control where rows are target angles 
% to compare and cols muscle combinations (all)
for h = 1:length(fields_coh)
    for i = 1:length(Aparams.angComp)
        for j = 1:nmusccomb
            coh_lim.(fields_coh{h})(i,j) = max(max(trial_data_coh_force(Aparams.targetAnglesForce == Aparams.angComp{i}(1)).(field).my_coh(3:end-3,j)),max(trial_data_coh_EMG(Aparams.targetAnglesEMG == Aparams.angComp{i}(2)).(field).my_coh(3:end-3,j)));
        end
    end
end

%% Coherence 
%% Fig per muscle pair in plotmusc. Subplot per target from angComp. Line per task.
% Compare coherence of each muscle combination (plotmusc) between tasks and targets from angComp
for j = 1:length(plotmusc)
    figure('Name',['My Coherence: Musc: ',trial_data_coh_force(1).(field).muscles{plotmusc(j)}{1},...
        ',',trial_data_coh_force(1).(field).muscles{plotmusc(j)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(length(Aparams.angComp),1,i);
        end
        h1 = plot(trial_data_coh_force(iangf).(field).my_fcoh(:,plotmusc(j)),trial_data_coh_force(iangf).(field).my_coh(:,plotmusc(j)));
        hold on;
        line(xlim,trial_data_coh_force(iangf).(field).my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_data_coh_EMG(iangE).(field).my_fcoh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).my_coh(:,plotmusc(j)),'r');
        hold on;
        line(xlim,trial_data_coh_EMG(iangE).(field).my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','-.');
        
        xlim([2 fc]); ylim([0 coh_lim.(field)(i,plotmusc(j))]);
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),' deg (',Aparams.muscComp{i},')']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end

%% Fig per muscle pair (all). Subplot per target from angComp. Line per task.
% Compare coherence of each muscle combination (all) between tasks and targets from angComp
for j = 1:nmusccomb
    figure('Name','My Coherence');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(length(Aparams.angComp),1,i);
        end
        h1 = plot(trial_data_coh_force(iangf).(field).my_fcoh(:,j),trial_data_coh_force(iangf).(field).my_coh(:,j));
        hold on;
        line(xlim,trial_data_coh_force(iangf).(field).my_CL(j)*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_data_coh_EMG(iangE).(field).my_fcoh(:,j),trial_data_coh_EMG(iangE).(field).my_coh(:,j),'r');
        hold on;
        line(xlim,trial_data_coh_EMG(iangE).(field).my_CL(j)*[1 1],'Color','k','LineStyle','-.');
                
        xlim([trial_data_coh_force(iangf).(field).my_fcoh(2,j) fc]);
        %ylim([0 1]);
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_data_coh_force(iangf).(field).muscles{j}{1},...
            ',',trial_data_coh_force(iangf).(field).muscles{j}{2},...
            '; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end

%% Fig per muscle pair (plotmusc). Subplot per task (col) per target from angComp (row). Line per processing method.
% Compare processing method (filt/rect) for each muscle combination between
% task and target from angComp
for j = 1:length(plotmusc)
    figure('Name',['My Coherence - Comparison Proc; Musc: ',trial_data_coh_force(1).filt.muscles{plotmusc(j)}{1},...
        ',',trial_data_coh_force(1).filt.muscles{plotmusc(j)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    h = 1;
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        subplot(length(Aparams.angComp),2,h);
        h1 = plot(trial_data_coh_force(iangf).filt.my_fcoh(:,plotmusc(j)),trial_data_coh_force(iangf).filt.my_coh(:,plotmusc(j)),'b');
        hold on;
        line(xlim,trial_data_coh_force(iangf).filt.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','--');
        h2 = plot(trial_data_coh_force(iangf).rect.my_fcoh(:,plotmusc(j)),trial_data_coh_force(iangf).rect.my_coh(:,plotmusc(j)),'c');
        line(xlim,trial_data_coh_force(iangf).rect.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','-.');
        xlim([2 fc]); ylim([0 max(coh_lim.filt(i,plotmusc(j)),coh_lim.rect(i,plotmusc(j)))]);        
        if i == length(Aparams.angComp)
            xlabel('Frequency [Hz]'); 
        end
        ylabel('Coh [-]');
        title(['ForceCO; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'Filt','Rect')
        
        subplot(length(Aparams.angComp),2,h+1);
        h1 = plot(trial_data_coh_EMG(iangE).filt.my_fcoh(:,plotmusc(j)),trial_data_coh_EMG(iangE).filt.my_coh(:,plotmusc(j)),'r');
        hold on;
        line(xlim,trial_data_coh_EMG(iangE).filt.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','--');
        h2 = plot(trial_data_coh_EMG(iangE).rect.my_fcoh(:,plotmusc(j)),trial_data_coh_EMG(iangE).rect.my_coh(:,plotmusc(j)),'m');
        line(xlim,trial_data_coh_EMG(iangE).rect.my_CL(plotmusc(j))*[1 1],'Color','k','LineStyle','-.');
        xlim([2 fc]); ylim([0 max(coh_lim.filt(i,plotmusc(j)),coh_lim.rect(i,plotmusc(j)))]);        
        if i == length(Aparams.angComp)
            xlabel('Frequency [Hz]'); 
        end
        ylabel('Coh [-]');
        title(['EMGCO; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'Filt','Rect')
        
        h = h+2;
    end
end

%% Significant coherence
%% Fig of mean significant coherence per muscle pair (plotmusc). Subplot per target from angComp. Bar per frequency band and task.
% Compare mean significant coherence (SEM for CI) for each muscle
% combination (plotmusc) between task and frequency bands 
bandCol = {'c','g','r'};
freqBand = {'Alpha','Beta','Gamma'};

for j = 1:length(plotmusc)
    figure('Name',['Significant coherence: Musc: ',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{1},...
        ',',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));

        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(1,length(Aparams.angComp),i);
        end
        bar([1 2 3],[trial_data_coh_force(iangf).(field).my_mean_sig_coh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).my_mean_sig_coh(:,plotmusc(j))],'barwidth',0.9);
        hold on;
        errorbar([1-0.15 1+0.15; 2-0.15 2+0.15; 3-0.15 3+0.15],[trial_data_coh_force(iangf).(field).my_mean_sig_coh(:,plotmusc(j)),...
            trial_data_coh_EMG(iangE).(field).my_mean_sig_coh(:,plotmusc(j))],...
            [trial_data_coh_force(iangf).(field).my_SEM_sig_coh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).my_SEM_sig_coh(:,plotmusc(j))],'k.');%,'facecolor',bandCol{i},'barwidth',0.9);
        ylim([0 1]);
        xlabel('Frequency [Hz]'); ylabel('Sig Coh [-]');
        title(['Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),' deg (',Aparams.muscComp{i},')']);
        set(gca,'XTick',1:3,'XTickLabel',freqBand);
        legend('ForceCO','EMGCO');
    end
end

%% Fig of percent of significant coherence per muscle pair (plotmusc). Subplot per target from angComp. Bar per frequency band and task.
% Compare percent of significant coherence of each frequency band for each 
% muscle pair (plotmusc) between task and frequency band
freqBand = {'Alpha','Beta','Gamma'};

for j = 1:length(plotmusc)
    figure('Name',['Significant coherence: Musc: ',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{1},...
        ',',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{2}]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        %mean_coh_alp(i,:) = [trial_data_coh_force(iangf).(field).sig_coh(1,plotmusc(j))
        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(1,length(Aparams.angComp),i);
        end
        bar([1 2 3],[trial_data_coh_force(iangf).(field).my_perc_sig_coh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).my_perc_sig_coh(:,plotmusc(j))],'barwidth',0.9);
        hold on;
        errorbar([1-0.15 1+0.15; 2-0.15 2+0.15; 3-0.15 3+0.15],[trial_data_coh_force(iangf).(field).my_perc_sig_coh(:,plotmusc(j)),...
            trial_data_coh_EMG(iangE).(field).my_perc_sig_coh(:,plotmusc(j))],...
            [trial_data_coh_force(iangf).(field).my_perc_CI_sig_coh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).my_perc_CI_sig_coh(:,plotmusc(j))],'k.');%,'facecolor',bandCol{i},'barwidth',0.9);
        %ylim([0 1]);
        xlabel('Frequency [Hz]'); ylabel('Perc Sig Coh [-]');
        title(['Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),' deg (',Aparams.muscComp{i},')']);
        set(gca,'XTick',1:3,'XTickLabel',freqBand);
        legend('ForceCO','EMGCO');
    end
end

%% Correlation
%% Fig per muscle pair (plotmusc) and task. Subplot per target from angComp.
% Compare correlation between targets from angComp for each task and muscle
% combination (plotmusc)
for i = 1:length(plotmusc)
    figure('Name','ForceCO');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for j = 1:length(Aparams.targetAnglesForce)
        if rem(length(Aparams.targetAnglesForce),2) == 0
            subplot(2,length(Aparams.targetAnglesForce)/2,j);
        else
            subplot(1,length(Aparams.targetAnglesForce),j);
        end
        plot(trial_data_coh_force(j).(field).lags(:,plotmusc(i))*trial_data_avg_force(j).ts(2),trial_data_coh_force(j).(field).corr(:,plotmusc(i)));
        hold on;
        xlim([-0.1 0.1])
        xlabel('Time lag [s]'); ylabel('Corr [-]');
        title(['Musc: ',trial_data_coh_force(j).(field).muscles{plotmusc(i)}{1},...
            ',',trial_data_coh_force(j).(field).muscles{plotmusc(i)}{2},...
            ' ; Target: ',num2str(rad2deg(Aparams.targetAnglesForce(j))),' deg']);
    end
end

for i = 1:length(plotmusc)
    figure('Name','EMGCO');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for j = 1:length(Aparams.targetAnglesEMG)
        if rem(length(Aparams.targetAnglesEMG),2) == 0
            subplot(2,length(Aparams.targetAnglesEMG)/2,j);
        else
            subplot(1,length(Aparams.targetAnglesEMG),j);
        end
        plot(trial_data_coh_EMG(j).(field).lags(:,plotmusc(i))*trial_data_avg_EMG(j).ts(2),trial_data_coh_EMG(j).(field).corr(:,plotmusc(i)));
        hold on;
        xlim([-0.1 0.1])
        xlabel('Time lag [s]'); ylabel('Corr [-]');
        title(['Musc: ',trial_data_coh_EMG(j).(field).muscles{plotmusc(i)}{1},...
            ',',trial_data_coh_EMG(j).(field).muscles{plotmusc(i)}{2},...
            ' ; Target: ',num2str(rad2deg(Aparams.targetAnglesEMG(j))),' deg']);
    end
end

%% MATLAB coherence
%% Fig per muscle pair (plotmusc). Subplot per target from angComp. Line per task.
% Compare coherence of each muscle combination (plotmusc) between tasks and 
% targets from angComp
for j = 1:length(plotmusc)
    figure('Name','MATLAB Coherence');
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(length(Aparams.angComp),1,i);
        end
        h1 = plot(trial_data_coh_force(iangf).(field).fcoh(:,plotmusc(j)),trial_data_coh_force(iangf).(field).coh(:,plotmusc(j)));
        hold on;
        line(xlim,trial_data_coh_force(iangf).(field).CL(plotmusc(j))*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_data_coh_EMG(iangE).(field).fcoh(:,plotmusc(j)),trial_data_coh_EMG(iangE).(field).coh(:,plotmusc(j)),'r');
        hold on;
        line(xlim,trial_data_coh_EMG(iangE).(field).CL(plotmusc(j))*[1 1],'Color','k','LineStyle','-.');
        
        xlim([trial_data_coh_force(iangf).(field).fcoh(2,plotmusc(j)) fc]);
        ylim([0 0.8]);
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{1},...
            ',',trial_data_coh_force(iangf).(field).muscles{plotmusc(j)}{2},...
            '; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end

%% Fig per muscle pair (all). Subplot per target from angComp. Line per task.
% Compare coherence of each muscle combination (all) between tasks and 
% targets from angComp
for j = 1:nmusccomb
    figure('Name','Coherence');
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));

        if rem(length(Aparams.angComp),2) == 0
            subplot(length(Aparams.angComp)/2,2,i);
        else
            subplot(1,length(Aparams.angComp),i);
        end
        h1 = plot(trial_data_coh_force(iangf).(field).fcoh(:,j),trial_data_coh_force(iangf).(field).coh(:,j));
        hold on;
        line(xlim,trial_data_coh_force(iangf).(field).CL(j)*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_data_coh_EMG(iangE).(field).fcoh(:,j),trial_data_coh_EMG(iangE).(field).coh(:,j),'r');
        hold on;
        line(xlim,trial_data_coh_EMG(iangE).(field).CL(j)*[1 1],'Color','k','LineStyle','-.');
        
        xlim([trial_data_coh_force(iangf).(field).fcoh(2,j) fc])
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_data_coh_force(iangf).(field).muscles{j}{1},...
            ',',trial_data_coh_force(iangf).(field).muscles{j}{2},...
            '; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',Aparams.muscComp{i},') deg']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end

%% Scrap
% %% Save trial data structs
% if ~exist([filepath,'TrialData/'],'dir')
%     mkdir([filepath,'TrialData/']);
% end
% if exist('trial_data_EMG_calib','var')
%     save([filepath,'TrialData/trial_data_calib'],'trial_data_EMG_calib','trial_data_avg_calib');
% end
% save([filepath,'TrialData/trial_data'],'trial_data_force','trial_data_EMG','Aparams');
% save([filepath,'TrialData/trial_data_avg'],'trial_data_avg_force','trial_data_avg_EMG');
% save([filepath,'TrialData/trial_data_app'],'trial_data_app_force','trial_data_app_EMG');

% for i = 1:6
%     for kk = 1:4
%         subplot(2,3,i)
%         h(kk)=bar(kk,mpt(kk,i),'facecolor',cols{kk},'barwidth',0.9); hold on;
%         errorbar(kk,mpt(kk,i),spt(kk,i),'k.')
%         title(strcat(header{i+4},' (\mu,\sigma)'));ylabel('Score');xlim([0 5]);
%         set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
%     end
% end