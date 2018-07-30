%% Data Analysis
addpath('Tools');

date =      '20180726';
task =      'EMGCO';
code =      '002';
EMG =       1;
filenameforce =  [date,'_',task,'_Force_',code,'.mat'];
filenameEMG = [date,'_',task,'_EMG_',code,'.mat'];
filenameparams = [date,'_params_',code,'.mat'];

switch computer
    case 'PCWIN'
        filepath =  ['D:\Student_experiments\Virginia\Data\' date '\'];
    case 'MACI64'
        filepath =  ['/Users/virginia/Documents/MATLAB/Thesis/Data/' date '/'];
end

load([filepath,filenameforce]);
if EMG
    load([filepath,filenameEMG]);
end 
load([filepath,'Parameters/',filenameparams]);

Aparams.downsample = forceparams.scanRate/EMGparams.sampleRateEMG;
Aparams.channelNameEMG = EMGparams.channelName;
if strcmp(task,'ForceCO')
    forceEMGData = {forceDataOut_ForceCO,EMGDataOut_ForceCO};
    Aparams.target_angles = [0:2*pi/8:2*pi-pi/8];
elseif strcmp(task,'EMGCO')
    forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};
    Aparams.target_angles = [pi/4:pi/4:3*pi/4];
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
fields = {'EMG.rect','EMG.avg','force.filt'};
trial_data_avg = trialAngleAvg(trial_data, epoch, fields);
trial_data_app = trialAngleApp(trial_data, epoch, fields);

%% EMG 
data_analysis = trial_data_app;

%% Frequency analysis
cohparams.nseg = 10;
cohparams.window = @(N) hanning(N);
cohparams.overlap = length(cohparams.window)/2;
trial_data_coh = cohStruct(data_analysis,EMGparams.channelName,{'rect'},cohparams);

target_angles = sort(unique(extractfield(data_analysis,'angle')));
nangles = length(target_angles);
nmusc = length(EMGparams.channelName)-1;
nmusccomb = (length(EMGparams.channelName)-1)*(length(EMGparams.channelName)-2)/2; % n*(n-1)/2

fc = 500;

figure;
set(gcf,'Name','FFT');
for i = 1:nangles
    if rem(nangles,2) == 0
        subplot(2,nangles/2,i);
    else
        subplot(1,nangles,i);
    end
    for j = 1:nmusc
        plot(data_analysis(i).fv,abs(data_analysis(i).EMG.rect_fft(:,j))/length(data_analysis(i).EMG.rect_fft(:,j)));
        hold on;
    end
    xlim([data_analysis(i).fv(2) fc])% data_analysis(i).fv(round(end/2))])
    ylim([0 max(max(abs(data_analysis(i).EMG.rect_fft(2:end,:))/size(data_analysis(i).EMG.rect_fft,1)))+1]);
    xlabel('Frequency [Hz]'); ylabel('FFT [-]');
    title(['Target: ',num2str(rad2deg(trial_data_coh(i).angle)),' deg']);
    legend(EMGparams.channelName{1:end-1});
end

for j = 1:nmusccomb
    figure;
    set(gcf,'Name','Coherence');
    for i = 1:nangles
        if rem(nangles,2) == 0
            subplot(2,nangles/2,i);
        else
            subplot(1,nangles,i);
        end
        plot(trial_data_coh(i).rect.fcoh(:,j),trial_data_coh(i).rect.coh(:,j));
        hold on;
        line(xlim,trial_data_coh(i).rect.CL(j)*[1 1]);
        xlim([trial_data_coh(i).rect.fcoh(1,j) fc])
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_data_coh(i).rect.muscles{j}{1},',',trial_data_coh(i).rect.muscles{j}{2},'; Target: ',num2str(rad2deg(trial_data_coh(i).angle)),' deg']);
    end
end

%% Time analysis
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

%% Force analysis - only with trial_data_avg
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
    ylim(taskparams.targetForce*[-1.5 1.5]);
    xlabel('Time [s]'); ylabel('Force [N]');
    title(['Target: ',num2str(rad2deg(trial_data_avg(i).angle)),' deg']);
    legend('Fx','Fy')
end

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
    xlim(taskparams.targetForce*[-1.5 1.5]);ylim(taskparams.targetForce*[-1.5 1.5]);
    grid on;
    xlabel('Fx [N]'); ylabel('Fy [N]');
    title(['Target: ',num2str(rad2deg(trial_data_avg(i).angle)),' deg']);
end

%% Single-trial analysis
itrial = 4;

%% EMG
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

figure;
set(gcf,'Name','Force trajectory');
plot(trial_data(itrial).force.filt(:,1),trial_data(itrial).force.filt(:,2));
hold on;
h1 = plot(trial_data(itrial).force.filt(1,1),trial_data(itrial).force.filt(1,2),'go');
h2 = plot(trial_data(itrial).force.filt(end,1),trial_data(itrial).force.filt(end,2),'ro');
%xlim(taskparams.targetForce*[-1.5 1.5]);ylim(taskparams.targetForce*[-1.5 1.5]);
grid on;
axis square;
xlabel('Fx [N]'); ylabel('Fy [N]');
title(['Trial: ',num2str(itrial),'; Target: ',num2str(rad2deg(trial_data(itrial).angle)),' deg']);
legend([h1,h2],'Start','End');

%% Mapping
N = size(trial_data(4).EMGrect,1);
tv = (0:N-1)*dt;

FD = trial_data(4).forcefilt;
ED = trial_data(4).EMGavg;

bx = regress(FD(:,1),[ones(size(ED,1),1) ED]); 
by = regress(FD(:,2),[ones(size(ED,1),1) ED]); 

fitx = bx(1)+bx(2)*ED(:,1)+bx(3)*ED(:,2);
fity = by(1)+by(2)*ED(:,1)+by(3)*ED(:,2);

figure 
subplot(121)
plot(tv,FD(:,1))
hold on
plot(tv,fitx)

subplot(122)
plot(tv,FD(:,2))
hold on
plot(tv,fity)

