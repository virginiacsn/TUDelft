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

epoch = {'imove','iend'};
fields = {'EMGrect','forcefilt','forcemag'};
trial_data_avg = trialangleavg(trial_data, epoch, fields);

%% Single-trial analysis
itrial = 4;

fs = 1024;
dt = 1/fs;
N = size(trial_data(itrial).EMGrect,1);
tv = (0:N-1)*dt;
fv = (0:N-1)/(N*dt);
fc = 50;

forceFiltx = trial_data(itrial).forcefilt(:,1);
forceFilty = trial_data(itrial).forcefilt(:,2);

figure
% subplot(1,3,1)
% plot(time,trial_data(itrial).forcemag)
% title('|F|')
% xlabel('Time [s]'); ylabel('Force [N]')

subplot(1,2,1)
plot(tv,forceFiltx);
hold on
plot(tv,forceFilty);
legend('Fx','Fy')
title('Fx, Fy')
xlabel('Time [s]'); ylabel('Force [N]')

subplot(1,2,2)
plot(forceFiltx,forceFilty);
hold on
plot(forceFiltx(1),forceFilty(1),'go');
plot(forceFiltx(end),forceFilty(end),'ro');
title('Trajectory')
xlabel('Fx [N]'); ylabel('Fy [N]')
axis equal

%% EMG
EMGfft = fft(trial_data(itrial).EMGrect)/N;
[EMGcoh,fcoh] = mscohere(trial_data(itrial).EMGrect(:,1),trial_data(itrial).EMGrect(:,2),hamming(128),[],1024);
L = length(trial_data(itrial).EMGrect(:,1))/128;

figure
subplot(121)
plot(tv,trial_data(itrial).EMGrect(:,1));
hold on
plot(tv,trial_data(itrial).EMGavg(:,1));
legend('Rectified','Smoothed')
xlabel('Time [s]'); ylabel('EMG [-]')
title('BB')

subplot(122)
plot(tv,trial_data(itrial).EMGrect(:,2));
hold on
plot(tv,trial_data(itrial).EMGavg(:,2));
legend('Rectified','Smoothed')
xlabel('Time [s]'); ylabel('EMG [-]')
title('TL')

figure
plot(fv,abs(EMGfft))
legend('BB','TL')
title('FFT')
xlabel('Frequency [Hz]'); ylabel('Power [-]')
xlim([1/(N*dt) fc])

% Coherence
alp = 0.05;
Z = 1-alp^(1/(L-1));

figure
plot(fcoh,EMGcoh)
hold on
line(xlim,[Z Z])
%xlim([1/(N*dt) fc])
title('Coherence')
xlabel('Frequency [Hz]'); ylabel('Coherence [-]')

%% Trial-average analysis
target_angles = [0:2*pi/8:2*pi-pi/8];
nangles = length(target_angles);
j = 0;

figure

for i = 1:nangles
    angle = target_angles(i);
    iangle = find(extractfield(trial_data_avg,'angle') ==  angle);
    
    fs = 1024;
    dt = 1/fs;
    N = size(trial_data_avg(iangle).EMGrect_time,1);
    tv = (0:N-1)*dt;
    fv = (0:N-1)/(N*dt);
    fc = 50;
    
    subplot(nangles,2,i*2-1)
    plot(tv,trial_data_avg(iangle).forcefilt_time(:,1));
    hold on
    plot(tv,trial_data_avg(iangle).forcefilt_time(:,2));
    legend('Fx','Fy')
    title('Fx, Fy')
    xlabel('Time [s]'); ylabel('Force [N]')
    
    subplot(nangles,2,i*2)
    plot(trial_data_avg(iangle).forcefilt_time(:,1),trial_data_avg(iangle).forcefilt_time(:,2));
    hold on
    plot(trial_data_avg(iangle).forcefilt_time(1,1),trial_data_avg(iangle).forcefilt_time(1,2),'go');
    plot(trial_data_avg(iangle).forcefilt_time(end,1),trial_data_avg(iangle).forcefilt_time(end,2),'ro');
    title('Trajectory')
    xlabel('Fx [N]'); ylabel('Fy [N]')
    axis equal
end

%% EMG
avgWindow = 100;
figure

for i = 1:nangles
    angle = target_angles(i);
    
    iangle = find(extractfield(trial_data_avg,'angle') ==  angle);
    
    subplot(nangles,2,i*2-1)
    plot(tv,trial_data_avg(iangle).EMGrect_time(:,1));
    hold on
    plot(tv,EMGavg(:,1));
    %legend('Rectified','Smoothed')
    xlabel('Time [s]'); ylabel('EMG [-]')
    title('BB') 
    
    subplot(nangles,2,i*2)
    plot(tv,trial_data_avg(iangle).EMGrect_time(:,2));
    hold on
    plot(tv,EMGavg(:,2));
    %legend('Rectified','Smoothed')
    xlabel('Time [s]'); ylabel('EMG [-]')
    title('TL')
    
end

figure
plot(fv,abs(EMGfft))
legend('BB','TL')
title('FFT')
xlabel('Frequency [Hz]'); ylabel('Power [-]')
xlim([1/(N*dt) fc])

% Coherence
alp = 0.05;
Z = 1-alp^(1/(L-1));

EMGavg = movingavg(trial_data_avg(iangle).EMGrect_time,avgWindow);
EMGfft = fft(trial_data_avg(iangle).EMGrect_time)/N;
[EMGcoh,fcoh] = mscohere(trial_data_avg(iangle).EMGrect_time(:,1),trial_data_avg(iangle).EMGrect_time(:,2),hamming(128),[],1024);
L = length(trial_data_avg(iangle).EMGrect_time(:,1))/128;

figure
plot(fcoh/pi*fs,EMGcoh)
hold on
line(xlim,[Z Z])
%xlim([1/(N*dt) fc])
title('Coherence')
xlabel('Frequency [Hz]'); ylabel('Coherence [-]')

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


