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

forceEMGData = {forceDataOut,EMGDataOut};

Aparams.downsample = 2;
if strcmp(task,'ForceCO')
    Aparams.target_angles = [0:2*pi/8:2*pi-pi/8];
elseif strcmp(task,'EMGCO')
    Aparams.target_angles = [pi/4:pi/4:3*pi/4];
end
Aparams.avgWindow = 200;
Aparams.fclF = 5;
Aparams.fchEMG = 10;

trial_data = trialCO(forceEMGData,Aparams);

trial_data = removeFailTrials(trial_data);

trial_data = procEMG(trial_data,Aparams);
trial_data = procForce(trial_data,Aparams);

%% Plot
itrial = 2;
time = [trial_data(itrial).dt:trial_data(itrial).dt:trial_data(itrial).dt*trial_data(itrial).iend];
timeavg = [trial_data(itrial).dt:trial_data(itrial).dt:trial_data(itrial).dt*size(trial_data(itrial).EMGavg,1)];

figure;
% subplot(1,3,1)
% plot(time,trial_data(itrial).forcemag)
% title('|F|')
% xlabel('Time [s]'); ylabel('Force [N]')

subplot(1,2,1)
plot(time,trial_data(itrial).forcefilt(:,1));
hold on
plot(time,trial_data(itrial).forcefilt(:,2));
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
figure
plot(time,trial_data(itrial).EMGrect(:,1));
hold on
plot(timeavg,trial_data(itrial).EMGavg(:,1));
legend('Rectified','Smoothed')
xlabel('Time [s]'); ylabel('EMG [mV]')

