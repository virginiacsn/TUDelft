%% Data Analysis
addpath('Tools');

date =      '20180723';
task =      'ForceCO';
code =      '001';
EMG =       1;
filenameforce =  [date,'_',task,'_Force_',code,'.mat'];
filenameEMG = [date,'_',task,'_EMG_',code,'.mat'];

switch computer
    case 'PCWIN'
        filepath =  ['D:\Student_experiments\Virginia\Data\' date '\'];
    case 'MACI64'
        filepath =  ['/Users/virginia/Documents/MATLAB/TUDELFT/Data/'];
end

load([filepath,filenameforce]);
if EMG
    load([filepath,filenameEMG]);
end 

Aparams.downsample = 2;
if strcmp(task,'ForceCO')
    forceEMGData = {forceDataOut,EMGDataOut};
    Aparams.target_angles = [0:2*pi/8:2*pi-pi/8];
elseif strcmp(task,'EMGCO')
    forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};
    Aparams.target_angles = [pi/4:pi/4:3*pi/4];
end
Aparams.avgWindow = 200;
Aparams.fclF = 5;
Aparams.fchEMG = 10;

trial_data = trialCO(forceEMGData,Aparams);

trial_data = removeFailTrials(trial_data);

trial_data = procEMG(trial_data,Aparams);
trial_data = procForce(trial_data,Aparams);

epoch = {'ihold','iend'};
fields = {'EMGrect','forcefilt','forcemag'};
trial_data_avg = trialangleavg(trial_data, epoch, fields);

%% Single-trial analysis
itrial = 1;

fs = 1024;
dt = 1/fs;
N = size(trial_data(itrial).EMGrect,1);
tv = (0:N-1)*dt;
fv = (0:N-1)/(N*dt);
fc = 50;

figure;
% subplot(1,3,1)
% plot(time,trial_data(itrial).forcemag)
% title('|F|')
% xlabel('Time [s]'); ylabel('Force [N]')

subplot(1,2,1)
plot(tv,trial_data(itrial).forcefilt(:,1));
hold on
plot(tv,trial_data(itrial).forcefilt(:,2));
legend('Fx','Fy')
title('Fx, Fy')
xlabel('Time [s]'); ylabel('Force [N]')

subplot(1,2,2)
plot(trial_data(itrial).forcefilt(:,1),trial_data(itrial).forcefilt(:,2));
hold on
plot(trial_data(itrial).forcefilt(1,1),trial_data(itrial).forcefilt(1,2),'go');
plot(trial_data(itrial).forcefilt(end,1),trial_data(itrial).forcefilt(end,2),'ro');
title('Trajectory')
xlabel('Fx [N]'); ylabel('Fy [N]')
axis equal

%% EMG
EMGfft = fft(trial_data(itrial).EMGrect)/N;
[coh,fcoh] = mscohere(trial_data(itrial).EMGrect(:,1),trial_data(itrial).EMGrect(:,2),hamming(128),[],1024);

figure
subplot(121)
plot(tv,trial_data(itrial).EMGrect(:,1));
hold on
plot(tv,trial_data(itrial).EMGavg(:,1));
legend('Rectified','Smoothed')
xlabel('Time [s]'); ylabel('EMG [mV]')
title('BB')

subplot(122)
plot(tv,trial_data(itrial).EMGrect(:,2));
hold on
plot(tv,trial_data(itrial).EMGavg(:,2));
legend('Rectified','Smoothed')
xlabel('Time [s]'); ylabel('EMG [mV]')
title('TL')

figure
plot(fv,abs(EMGfft))
legend('BB','TL')
title('FFT')
xlabel('Frequency [Hz]'); ylabel('Power [-]')
xlim([1/(N*dt) fc])

%spectrogram(trial_data(itrial).EMGrect(:,2),256);
alp = 0.05;
Z = 1-alp^(1/(length(coh)-1));
figure
plot(fcoh,coh)
hold on
line(xlim,[Z Z])
%xlim([1/(N*dt) fc])
title('Coherence')
xlabel('Frequency [Hz]'); ylabel('Coherence [-]')

%% Trial-average analysis
angle = 0;
iangle = find(extractfield(trial_data_avg,'angle') ==  angle);

fs = 1024;
dt = 1/fs;
N = size(trial_data_avg(iangle).EMGrect_time,1);
tv = (0:N-1)*dt;
fv = (0:N-1)/(N*dt);
fc = 50;

figure;
subplot(1,2,1)
plot(tv,trial_data_avg(iangle).forcefilt_time(:,1));
hold on
plot(tv,trial_data_avg(iangle).forcefilt_time(:,2));
legend('Fx','Fy')
title('Fx, Fy')
xlabel('Time [s]'); ylabel('Force [N]')

subplot(1,2,2)
plot(trial_data_avg(iangle).forcefilt_time(:,1),trial_data_avg(iangle).forcefilt_time(:,2));
hold on
plot(trial_data_avg(iangle).forcefilt_time(1,1),trial_data_avg(iangle).forcefilt_time(1,2),'go');
plot(trial_data_avg(iangle).forcefilt_time(end,1),trial_data_avg(iangle).forcefilt_time(end,2),'ro');
title('Trajectory')
xlabel('Fx [N]'); ylabel('Fy [N]')
axis equal

%% EMG
avgWindow = 100;
EMGavg_time = movingavg(trial_data_avg(iangle).EMGrect_time,avgWindow);
EMGfft_time = fft(trial_data_avg(iangle).EMGrect_time)/N;
[coh,fcoh] = mscohere(trial_data_avg(iangle).EMGrect_time(:,1),trial_data_avg(iangle).EMGrect_time(:,2),hamming(128),[],1024);

figure
subplot(121)
plot(tv,trial_data_avg(iangle).EMGrect_time(:,1));
hold on
plot(tv,EMGavg_time(:,1));
legend('Rectified','Smoothed')
xlabel('Time [s]'); ylabel('EMG [mV]')
title('BB')

subplot(122)
plot(tv,trial_data_avg(iangle).EMGrect_time(:,2));
hold on
plot(tv,EMGavg_time(:,2));
legend('Rectified','Smoothed')
xlabel('Time [s]'); ylabel('EMG [mV]')
title('TL')

figure
plot(fv,abs(EMGfft_time))
legend('BB','TL')
title('FFT')
xlabel('Frequency [Hz]'); ylabel('Power [-]')
xlim([1/(N*dt) fc])

%spectrogram(trial_data(itrial).EMGrect(:,2),256);
alp = 0.05;
Z = 1-alp^(1/(trial_data_avg(iangle).ntrials-1));
figure
plot(fcoh,coh)
hold on
line(xlim,[Z Z])
%xlim([1/(N*dt) fc])
title('Coherence')
xlabel('Frequency [Hz]'); ylabel('Coherence [-]')


