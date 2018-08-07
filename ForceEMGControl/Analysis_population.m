%% Data Analysis for population
addpath(genpath('Tools'));

switch computer
    case 'PCWIN'
        filepath =  ['D:\Student_experiments\Virginia\Data\'];
    case 'MACI64'
        filepath =  ['/Users/virginia/Documents/MATLAB/Thesis/Data/'];
end

datedir = dir(filepath);
for i = 2:length(datedir)
    subjectdir = dir([filepath,datedir(i).name]);
    for j = 1:length(subjectdir)
        
    end
end
date =      '20180806';
subject =   '01';
task =      'ForceCO';
code =      '003';
EMG =       1;
filenameforce =  [date,'_s',subject,'_',task,'_Force_',code,'.mat'];
filenameEMG = [date,'_s',subject,'_',task,'_EMG_',code,'.mat'];
filenameparams = [date,'_s',subject,'_','_params_',code,'.mat'];



load([filepath,filenameforce]);
if EMG
    load([filepath,filenameEMG]);
end 
load([filepath,'Parameters/',filenameparams]);

Aparams.downsample = forceparams.scanRate/EMGparams.sampleRateEMG;
Aparams.channelNameEMG = EMGparams.channelName;
if strcmp(task,'ForceCO')
    forceEMGData = {forceDataOut_ForceCO,EMGDataOut_ForceCO};
    Aparams.target_angles = taskparams.targetAnglesForce;
elseif strcmp(task,'EMGCO')
    forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};
    Aparams.target_angles = taskparams.targetAnglesEMG;
end
Aparams.avgWindow = 100;
Aparams.fclF = 5;
Aparams.fchEMG = 10;
Aparams.fs = min(forceparams.scanRate,EMGparams.sampleRateEMG);

trial_data = trialCO(forceEMGData,Aparams);

trial_data = removeFailTrials(trial_data);

trial_data = procEMG(trial_data,Aparams);
trial_data = procForce(trial_data,Aparams);

%% Trial-average
epoch = {'ihold','iend'};
fields = {'EMG.filt','EMG.rect','EMG.avg','force.filt'};
trial_data_avg = trialAngleAvg(trial_data, epoch, fields);
trial_data_app = trialAngleApp(trial_data, epoch, fields);

%% EMG 
data_analysis = trial_data_avg;

%% Frequency analysis
cohparams.tseg = 1;
cohparams.nseg = 10;
cohparams.my_nseg = 10;
cohparams.window = @(N) hanning(N);
field = 'rect';
trial_data_coh = cohStruct(data_analysis,EMGparams.channelName,{field},cohparams);

target_angles = sort(unique(extractfield(data_analysis,'angle')));
nangles = length(target_angles);
nmusc = length(EMGparams.channelName)-1;
nmusccomb = (length(EMGparams.channelName)-1)*(length(EMGparams.channelName)-2)/2; % n*(n-1)/2

fc = 80;

% FFT
figure;
set(gcf,'Name','FFT');
for i = 1:nangles
    if rem(nangles,2) == 0
        subplot(2,nangles/2,i);
    else
        subplot(1,nangles,i);
    end
    for j = 1:nmusc
        plot(data_analysis(i).fv,abs(data_analysis(i).EMG.([field,'_fft'])(:,j))/length(data_analysis(i).EMG.([field,'_fft'])(:,j)));
        hold on;
    end
    xlim([data_analysis(i).fv(2) fc])% data_analysis(i).fv(round(end/2))])
    ylim([0 max(max(abs(data_analysis(i).EMG.([field,'_fft'])(2:end,:))/size(data_analysis(i).EMG.([field,'_fft']),1)))+1]);
    xlabel('Frequency [Hz]'); ylabel('FFT [-]');
    title(['Target: ',num2str(rad2deg(trial_data_coh(i).angle)),' deg']);
    legend(EMGparams.channelName{1:end-1});
end

% Coherence
for j = 1:nmusccomb
    figure;
    set(gcf,'Name','Coherence');
    for i = 1:nangles
        if rem(nangles,2) == 0
            subplot(2,nangles/2,i);
        else
            subplot(1,nangles,i);
        end
        plot(trial_data_coh(i).(field).fcoh(:,j),trial_data_coh(i).(field).coh(:,j));
        hold on;
        line(xlim,trial_data_coh(i).(field).CL(j)*[1 1]);
        xlim([trial_data_coh(i).(field).fcoh(2,j) fc])
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_data_coh(i).(field).muscles{j}{1},',',trial_data_coh(i).(field).muscles{j}{2},'; Target: ',num2str(rad2deg(trial_data_coh(i).angle)),' deg']);
    end
end

% My coherence
for j = 1:nmusccomb
    figure;
    
    set(gcf,'Name','My coherence');
    for i = 1:nangles
        if rem(nangles,2) == 0
            subplot(2,nangles/2,i);
        else
            subplot(1,nangles,i);
        end
        plot(trial_data_coh(i).(field).my_fcoh(:,j),trial_data_coh(i).(field).my_coh(:,j));
        hold on;
        %plot(trial_data_coh(i).(field).my_fcoh(:,j),movingAvg(trial_data_coh(i).(field).my_coh(:,j),3));
        line(xlim,trial_data_coh(i).(field).my_CL(j)*[1 1]);
        xlim([trial_data_coh(i).(field).my_fcoh(2,j) fc])
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_data_coh(i).(field).muscles{j}{1},',',trial_data_coh(i).(field).muscles{j}{2},'; Target: ',num2str(rad2deg(trial_data_coh(i).angle)),' deg']);
    end
end

%% Time analysis
% Rectified and smoothed EMG
for j = 1:nmusc
    figure;
    set(gcf,'Name',['EMG ',EMGparams.channelName{j}]);
    for i = 1:nangles
        if rem(nangles,2) == 0
            subplot(2,nangles/2,i);
        else
            subplot(1,nangles,i);
        end
        plot(data_analysis(i).ts,data_analysis(i).EMG.rect(:,j));
        hold on;
        plot(data_analysis(i).ts,data_analysis(i).EMG.avg(:,j));
        ylim([0 max(data_analysis(i).EMG.rect(:))+50]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['Target: ',num2str(rad2deg(data_analysis(i).angle)),' deg']);
    end
end

%% Correlation
for j = 1:nmusccomb
    figure;
    set(gcf,'Name','Unbiased correlation');
    for i = 1:nangles
        if rem(nangles,2) == 0
            subplot(2,nangles/2,i);
        else
            subplot(1,nangles,i);
        end
        plot(trial_data_coh(i).(field).lags(:,j)*trial_data_avg(i).ts(2),trial_data_coh(i).(field).corr(:,j));
        hold on;
        xlim([-0.1 0.1])
        xlabel('Time lag [s]'); ylabel('Corr [-]');
        title(['Musc: ',trial_data_coh(i).(field).muscles{j}{1},',',trial_data_coh(i).(field).muscles{j}{2},'; Target: ',num2str(rad2deg(trial_data_coh(i).angle)),' deg']);
    end
end

%% Force analysis - only with trial_data_avg
% LPF force in time
figure;
set(gcf,'Name','Force in time');
for i = 1:nangles
    if rem(nangles,2) == 0
        subplot(2,nangles/2,i);
    else
        subplot(1,nangles,i);
    end
    for j = 1:2
        plot(trial_data_avg(i).ts,trial_data_avg(i).force.filt(:,j));
        hold on;
    end
    ylim(taskparams.targetForce*[-2 2]);
    xlabel('Time [s]'); ylabel('Force [N]');
    title(['Target: ',num2str(rad2deg(trial_data_avg(i).angle)),' deg']);
    legend('Fx','Fy')
end

% LPF force trajectory
figure;
set(gcf,'Name','Force trajectory');
for i = 1:nangles
    if rem(nangles,2) == 0
        subplot(2,nangles/2,i);
    else
        subplot(1,nangles,i);
    end
    plot(trial_data_avg(i).force.filt(:,1),trial_data_avg(i).force.filt(:,2));
    hold on;
    plot(trial_data_avg(i).force.filt(1,1),trial_data_avg(i).force.filt(1,2),'go');
    plot(trial_data_avg(i).force.filt(end,1),trial_data_avg(i).force.filt(end,2),'ro');
    xlim(taskparams.targetForce*[-2 2]);ylim(taskparams.targetForce*[-2 2]);
    grid on;
    xlabel('Fx [N]'); ylabel('Fy [N]');
    title(['Target: ',num2str(rad2deg(trial_data_avg(i).angle)),' deg']);
end