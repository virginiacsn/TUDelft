%% Data Analysis for individual subject
addpath(genpath('Tools'));

date =      '20180810';
subject =   '01';
task =      'ForceCO';

switch computer
    case 'PCWIN'
        filepath =  ['D:\Student_experiments\Virginia\Data\',date,'\s',subject,'\'];
        paramfolder = 'Parameters\';
    case 'MACI64'
        filepath =  ['/Users/virginia/Documents/MATLAB/Thesis/Data/',date,'/s',subject,'/'];
        paramfolder = 'Parameters/';
end

%% Individual block
code =      '001';

filenameforce =  [date,'_s',subject,'_',task,'_Force_',code,'.mat'];
filenameEMG = [date,'_s',subject,'_',task,'_EMG_',code,'.mat'];
filenameparams = [date,'_s',subject,'_params_',code,'.mat'];

load([filepath,filenameforce]);
load([filepath,filenameEMG]);
load([filepath,paramfolder,filenameparams]);

Aparams.downsample = forceparams.scanRate/EMGparams.sampleRateEMG;
Aparams.channelNameEMG = EMGparams.channelName;
if strcmp(task,'ForceCO')
    forceEMGData = {forceDataOut_ForceCO,EMGDataOut_ForceCO};
    Aparams.target_angles =  [pi/4:pi/2:7*pi/4];%taskparams.targetAnglesForce;
elseif strcmp(task,'EMGCO')
    forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};
    Aparams.target_angles = [pi/4:pi/4:3*pi/4];%taskparams.targetAnglesEMG;
end
Aparams.avgWindow = 100;
Aparams.fclF = 5;
Aparams.fchEMG = 10;
Aparams.fclEMG = 40;
Aparams.fs = min(forceparams.scanRate,EMGparams.sampleRateEMG);

trial_data = trialCO(forceEMGData,Aparams);

trial_data = removeFailTrials(trial_data);

trial_data = procEMG(trial_data,Aparams);
trial_data = procForce(trial_data,Aparams);

%% All blocks
code =      {'002','003'};

trial_data_block = [];

for i = 1:length(code)
    filenameforce =  [date,'_s',subject,'_',task,'_Force_',code{i},'.mat'];
    filenameEMG = [date,'_s',subject,'_',task,'_EMG_',code{i},'.mat'];
    filenameparams = [date,'_s',subject,'_params_',code{i},'.mat'];
    
    load([filepath,filenameforce]);
    load([filepath,filenameEMG]);
    load([filepath,'Parameters/',filenameparams]);
    
    Aparams.downsample = forceparams.scanRate/EMGparams.sampleRateEMG;
    Aparams.channelNameEMG = EMGparams.channelName;
    if strcmp(task,'ForceCO')
        forceEMGData = {forceDataOut_ForceCO,EMGDataOut_ForceCO};
        Aparams.target_angles = [pi/4:pi/2:7*pi/4];%taskparams.targetAnglesForce;
    elseif strcmp(task,'EMGCO')
        forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};
        Aparams.target_angles = [pi/4:pi/4:3*pi/4];%taskparams.targetAnglesEMG;
    end
    Aparams.avgWindow = 100;
    Aparams.fclF = 5;
    Aparams.fchEMG = 10;
    Aparams.fs = min(forceparams.scanRate,EMGparams.sampleRateEMG);
    Aparams.block = str2double(code{i});
    
    trial_data = trialCO(forceEMGData,Aparams);
    
    trial_data = removeFailTrials(trial_data);
    
    trial_data = procEMG(trial_data,Aparams);
    trial_data = procForce(trial_data,Aparams);
    
    trial_data_block =  [trial_data_block, trial_data];
end

%% Trial-average
epoch = {'ihold','iend'};
fields = {'EMG.raw','EMG.filt','EMG.rect','EMG.avg','force.filt'};
trial_data_avg = trialAngleAvg(trial_data, epoch, fields);
trial_data_app = trialAngleApp(trial_data, epoch, fields);

data_analysis = trial_data_avg;

target_angles = sort(unique(extractfield(data_analysis,'angle')));
nangles = length(target_angles);

%% EMG 
nmusc = length(EMGparams.channelName)-1;
musccont = EMGparams.channelControl;
nmusccomb = (length(EMGparams.channelName)-1)*(length(EMGparams.channelName)-2)/2; % n*(n-1)/2
fc = 80;

%% Frequency analysis
cohparams.tseg = 1;
cohparams.nseg = 10;
cohparams.my_nseg = 10;
cohparams.window = @(N) hanning(N);
field = 'rect';
trial_data_coh = cohStruct(data_analysis,EMGparams.channelName,{field},cohparams);

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
% for j = 1:nmusccomb
%     figure;
%     
%     set(gcf,'Name','My coherence');
%     for i = 1:nangles
%         if rem(nangles,2) == 0
%             subplot(2,nangles/2,i);
%         else
%             subplot(1,nangles,i);
%         end
%         plot(trial_data_coh(i).(field).my_fcoh(:,j),trial_data_coh(i).(field).my_coh(:,j));
%         hold on;
%         %plot(trial_data_coh(i).(field).my_fcoh(:,j),movingAvg(trial_data_coh(i).(field).my_coh(:,j),3));
%         line(xlim,trial_data_coh(i).(field).my_CL(j)*[1 1]);
%         xlim([trial_data_coh(i).(field).my_fcoh(2,j) fc])
%         xlabel('Frequency [Hz]'); ylabel('Coh [-]');
%         title(['Musc: ',trial_data_coh(i).(field).muscles{j}{1},',',trial_data_coh(i).(field).muscles{j}{2},'; Target: ',num2str(rad2deg(trial_data_coh(i).angle)),' deg']);
%     end
% end

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

%% Single-trial analysis
itrial = 8;

%% EMG
% Rectified and smoothed EMG
figure;
set(gcf,'Name','Rectified EMG');
for j = 1:nmusc
    subplot(1,nmusc,j);
    plot(trial_data(itrial).ts,trial_data(itrial).EMG.rect(:,j));
    hold on;
    plot(trial_data(itrial).ts,trial_data(itrial).EMG.avg(:,j));
    xlabel('Time [s]'); ylabel('EMG [-]');
    ylim([0 max(trial_data(itrial).EMG.rect(:))+50]);
    title(['Musc: ',EMGparams.channelName{j},'; Target: ',num2str(rad2deg(trial_data(itrial).angle)),' deg']);
end

% EMG trajectory
figure;
set(gcf,'Name','Trajectory Rectified EMG');
plot(trial_data(itrial).EMG.avg(:,musccont(1)),trial_data(itrial).EMG.avg(:,musccont(2)));
hold on;
h1 = plot(trial_data(itrial).EMG.avg(1,musccont(1)),trial_data(itrial).EMG.avg(1,musccont(2)),'go');
h2 = plot(trial_data(itrial).EMG.avg(end,musccont(1)),trial_data(itrial).EMG.avg(end,musccont(2)),'ro');
grid on;
axis square;
xlabel(EMGparams.channelName{musccont(1)}); ylabel(EMGparams.channelName{musccont(2)});
title(['Musc: ',EMGparams.channelName{musccont(1)},',',EMGparams.channelName{musccont(2)},'; Target: ',num2str(rad2deg(trial_data(itrial).angle)),' deg']);
legend([h1,h2],'Start','End');

% FFT of rectified EMG
figure;
set(gcf,'Name','FFT Rectified EMG');
for j = 1:nmusc
    plot(trial_data(itrial).fv,abs(trial_data(itrial).EMG.rect_fft(:,j))/length(trial_data(itrial).EMG.rect_fft(:,j)));
    hold on;
end
xlabel('Frequency [Hz]'); ylabel('FFT EMG [-]');
xlim([trial_data(itrial).fv(2) trial_data(itrial).fv(round(end/2))]);
title(['Target: ',num2str(rad2deg(trial_data(itrial).angle)),' deg']);
legend(EMGparams.channelName{1:end-1});

%% Force
% LPF force in time
figure;
set(gcf,'Name','Force in time');
for j = 1:2
    plot(trial_data(itrial).ts,trial_data(itrial).force.filt(:,j));
    hold on;
end
%ylim(taskparams.targetForce*[-1.5 1.5]);
xlabel('Time [s]'); ylabel('Force [N]');
title(['Trial: ',num2str(itrial),'; Target: ',num2str(rad2deg(trial_data(itrial).angle)),' deg']);
legend('Fx','Fy')

% LPF force trajectory
figure;
set(gcf,'Name','Force trajectory');
plot(trial_data(itrial).force.filt(:,1),trial_data(itrial).force.filt(:,2));
hold on;
h1 = plot(trial_data(itrial).force.filt(1,1),trial_data(itrial).force.filt(1,2),'go');
h2 = plot(trial_data(itrial).force.filt(end,1),trial_data(itrial).force.filt(end,2),'ro');
xlim(taskparams.targetForce*[-2 2]);ylim(taskparams.targetForce*[-2 2]);
grid on;
axis square;
xlabel('Fx [N]'); ylabel('Fy [N]');
title(['Trial: ',num2str(itrial),'; Target: ',num2str(rad2deg(trial_data(itrial).angle)),' deg']);
legend([h1,h2],'Start','End');

%% Mapping
% N = size(trial_data(4).EMGrect,1);
% tv = (0:N-1)*dt;
% 
% FD = trial_data(4).forcefilt;
% ED = trial_data(4).EMGavg;
% 
% bx = regress(FD(:,1),[ones(size(ED,1),1) ED]); 
% by = regress(FD(:,2),[ones(size(ED,1),1) ED]); 
% 
% fitx = bx(1)+bx(2)*ED(:,1)+bx(3)*ED(:,2);
% fity = by(1)+by(2)*ED(:,1)+by(3)*ED(:,2);
% 
% figure 
% subplot(121)
% plot(tv,FD(:,1))
% hold on
% plot(tv,fitx)
% 
% subplot(122)
% plot(tv,FD(:,2))
% hold on
% plot(tv,fity)
