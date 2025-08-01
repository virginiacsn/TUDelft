%% Data Analysis for individual subject to compare tasks
% close all
% clear all
addpath(genpath('Tools'));

date =      '20180810';
subject =   '01';

switch computer
    case 'PCWIN'
        filepath =  ['D:\Student_experiments\Virginia\Data\',date,'\s',subject,'\'];
        paramfolder = 'Parameters\';
    case 'MACI64'
        filepath =  ['/Users/virginia/Documents/MATLAB/Thesis/Data/',date,'/s',subject,'/'];
        paramfolder = 'Parameters/';
end

%% All blocks for ForceCO, individual block for EMGCO
codeF = {'001','002'};
codeE = '002';
tasks = {'ForceCO','EMGCO'};

% ForceCO
task =  'ForceCO';
trial_data_force = [];

for i = 1:length(codeF)
    
    code = codeF{i};
    
    filenameforce =  [date,'_s',subject,'_',task,'_Force_',code,'.mat'];
    filenameEMG = [date,'_s',subject,'_',task,'_EMG_',code,'.mat'];
    filenameparams = [date,'_s',subject,'_params_',code,'.mat'];
    
    load([filepath,filenameforce]);
    load([filepath,filenameEMG]);
    load([filepath,paramfolder,filenameparams]);
    
    forceEMGData = {forceDataOut_ForceCO,EMGDataOut_ForceCO};
    
    Aparams.downsamp = forceparams.scanRate/EMGparams.sampleRateEMG;
    Aparams.channelNameEMG = EMGparams.channelName;
    Aparams.target_angles =  taskparams.targetAnglesForce;
    Aparams.avgWindow = 200;
    Aparams.fclF = 5;
    Aparams.fchEMG = 20;
    Aparams.fclEMG = 60;
    Aparams.fs = min(forceparams.scanRate,EMGparams.sampleRateEMG);
    Aparams.block = str2double(code);

    trial_data = trialCO(forceEMGData,Aparams);
    trial_data = removeFailTrials(trial_data);
    
    trial_data = procEMG(trial_data,Aparams);
    trial_data_force = [trial_data_force, procForce(trial_data,Aparams)];
end

Aparams.targetAnglesForce = sort(unique(extractfield(trial_data_force,'angle')));

Aparams.epoch = {'ihold','iend'};
fields = {'EMG.raw','EMG.filt','force.filt'};

trial_data_avg_force = trialAngleAvg(trial_data_force, Aparams.epoch, fields);
trial_data_avg_force = procEMG(trial_data_avg_force,Aparams);

trial_data_app_force = trialAngleApp(trial_data_force, Aparams.epoch, fields,[]);
trial_data_app_force = procEMG(trial_data_app_force,Aparams);

% EMGCO
task = 'EMGCO';
code = codeE;

filenameforce =  [date,'_s',subject,'_',task,'_Force_',code,'.mat'];
filenameEMG = [date,'_s',subject,'_',task,'_EMG_',code,'.mat'];
filenameparams = [date,'_s',subject,'_params_',code,'.mat'];

load([filepath,filenameforce]);
load([filepath,filenameEMG]);
load([filepath,paramfolder,filenameparams]);

forceEMGData = {forceDataOut_EMGCO,EMGDataOut_EMGCO};

Aparams.downsamp = forceparams.scanRate/EMGparams.sampleRateEMG;
Aparams.channelNameEMG = EMGparams.channelName;
Aparams.target_angles = taskparams.targetAnglesEMG;
Aparams.avgWindow = 200;
Aparams.fclF = 5;
Aparams.fchEMG = 30;
Aparams.fclEMG = 60;
Aparams.fs = min(forceparams.scanRate,EMGparams.sampleRateEMG);
Aparams.block = str2double(code);

trial_data = trialCO(forceEMGData,Aparams);
trial_data = removeFailTrials(trial_data);

trial_data = procEMG(trial_data,Aparams);
trial_data_EMG = procForce(trial_data,Aparams);

Aparams.targetAnglesEMG = sort(unique(extractfield(trial_data_EMG,'angle')));

trial_data_avg_EMG = trialAngleAvg(trial_data_EMG, Aparams.epoch, fields);
trial_data_avg_EMG = procEMG(trial_data_avg_EMG,Aparams);

trial_data_app_EMG = trialAngleApp(trial_data_EMG, Aparams.epoch, fields, []);
trial_data_app_EMG = procEMG(trial_data_app_EMG,Aparams);

%% Force - only with trial_data_avg
%% LPF force in time for both tasks, check target angles
figure;
set(gcf,'Name','ForceCO; Force in time');
for i = 1:length(Aparams.targetAnglesForce)
    if rem(length(Aparams.targetAnglesForce),2) == 0
        subplot(2,length(Aparams.targetAnglesForce)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesForce),i);
    end
    for j = 1:2
        plot(trial_data_app_force(i).ts,trial_data_app_force(i).force.filt(:,j));
        hold on;
    end
    ylim(taskparams.targetForce*[-5 5]);
    xlabel('Time [s]'); ylabel('Force [N]');
    title(['Target: ',num2str(rad2deg(trial_data_app_force(i).angle)),' deg']);
    legend('Fx','Fy')
end

figure;
set(gcf,'Name','EMGCO; Force in time');
for i = 1:length(Aparams.targetAnglesEMG)
    if rem(length(Aparams.targetAnglesEMG),2) == 0
        subplot(2,length(Aparams.targetAnglesEMG)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesEMG),i);
    end
    for j = 1:2
        plot(trial_data_app_EMG(i).ts,trial_data_app_EMG(i).force.filt(:,j));
        hold on;
    end
    ylim(taskparams.targetForce*[-5 5]);
    xlabel('Time [s]'); ylabel('Force [N]');
    title(['Target: ',num2str(rad2deg(trial_data_app_EMG(i).angle)),' deg']);
    legend('Fx','Fy')
end

%% LPF force trajectory for both tasks, check target angles
figure;
set(gcf,'Name','ForceCO; Force trajectory');
for i = 1:length(Aparams.targetAnglesForce)
    if rem(length(Aparams.targetAnglesForce),2) == 0
        subplot(2,length(Aparams.targetAnglesForce)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesForce),i);
    end
    plot(trial_data_app_force(i).force.filt(:,1),trial_data_app_force(i).force.filt(:,2));
    hold on;
    h1 = plot(trial_data_app_force(i).force.filt(1,1),trial_data_app_force(i).force.filt(1,2),'go');
    h2 = plot(trial_data_app_force(i).force.filt(end,1),trial_data_app_force(i).force.filt(end,2),'ro');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-5 5]);ylim(taskparams.targetForce*[-5 5]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    title(['Target: ',num2str(rad2deg(trial_data_app_force(i).angle)),' deg']);
    legend([h1,h2],'Start','End');
end

figure;
set(gcf,'Name','EMGCO; Force trajectory');
for i = 1:length(Aparams.targetAnglesEMG)
    if rem(length(Aparams.targetAnglesEMG),2) == 0
        subplot(2,length(Aparams.targetAnglesEMG)/2,i);
    else
        subplot(1,length(Aparams.targetAnglesEMG),i);
    end
    plot(trial_data_app_EMG(i).force.filt(:,1),trial_data_app_EMG(i).force.filt(:,2));
    hold on;
    h1 = plot(trial_data_app_EMG(i).force.filt(1,1),trial_data_app_EMG(i).force.filt(1,2),'go');
    h2 = plot(trial_data_app_EMG(i).force.filt(end,1),trial_data_app_EMG(i).force.filt(end,2),'ro');
    axis square;
    grid on;
    xlim(taskparams.targetForce*[-5 5]);ylim(taskparams.targetForce*[-5 5]);
    xlabel('Fx [N]'); ylabel('Fy [N]');
    title(['Target: ',num2str(rad2deg(trial_data_app_EMG(i).angle)),' deg']);
    legend([h1,h2],'Start','End');
end

%% EMG 
%% Time analysis
%% Fig per muscle, subplot per target
for j = 1:length(EMGparams.channelName)-1
    figure('Name',['ForceCO; EMG ',EMGparams.channelName{j}]);
    for i = 1:length(taskparams.targetAnglesForce) 
        if rem(length(taskparams.targetAnglesForce) ,2) == 0
            subplot(2,length(taskparams.targetAnglesForce) /2,i);
        else
            subplot(1,length(taskparams.targetAnglesForce) ,i);
        end
        plot(trial_data_app_force(i).ts,trial_data_app_force(i).EMG.rect(:,j));
        hold on;
        plot(trial_data_app_force(i).ts,trial_data_app_force(i).EMG.avg(:,j));
        ylim([0 40]);%ylim([0 max(trial_data_app_force(i).EMG.rect(:))+50]);
        xlim([0 trial_data_app_force(i).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['Target: ',num2str(rad2deg(trial_data_app_force(i).angle)),' deg']);
    end
end

for j = 1:length(EMGparams.channelName)-1
    figure('Name',['EMGCO; EMG ',EMGparams.channelName{j}]);
    for i = 1:length(taskparams.targetAnglesEMG) 
        if rem(length(taskparams.targetAnglesEMG),2) == 0
            subplot(2,length(taskparams.targetAnglesEMG)/2,i);
        else
            subplot(1,length(taskparams.targetAnglesEMG),i);
        end
        plot(trial_data_app_EMG(i).ts,trial_data_app_EMG(i).EMG.rect(:,j));
        hold on;
        plot(trial_data_app_EMG(i).ts,trial_data_app_EMG(i).EMG.avg(:,j));
        ylim([0 40]);%ylim([0 max(trial_data_avg_EMG(i).EMG.rect(:))+50]);
        xlim([0 trial_data_app_EMG(i).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['Target: ',num2str(rad2deg(trial_data_app_EMG(i).angle)),' deg']);
    end
end

%% Fig per target, subplot per muscle
for j = 1:length(Aparams.targetAnglesForce)
    figure('Name',['ForceCO; Target: ',num2str(rad2deg(Aparams.targetAnglesForce(j))), ' deg']);
    for i = 1:length(EMGparams.channelName)-1
        if rem(length(EMGparams.channelName)-1,2) == 0
            subplot(2,(length(EMGparams.channelName)-1)/2,i);
        else
            subplot(1,length(EMGparams.channelName)-1,i);
        end
        plot(trial_data_app_force(j).ts,trial_data_app_force(j).EMG.rect(:,i));
        hold on;
        plot(trial_data_app_force(j).ts,trial_data_app_force(j).EMG.avg(:,i));
        ylim([0 40]);%ylim([0 max(trial_data_app_force(j).EMG.rect(:))+50]);
        xlim([0 trial_data_app_force(j).ts(end)]);
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
        plot(trial_data_app_EMG(j).ts,trial_data_app_EMG(j).EMG.rect(:,i));
        hold on;
        plot(trial_data_app_EMG(j).ts,trial_data_app_EMG(j).EMG.avg(:,i));
        ylim([0 40]);%ylim([0 max(trial_data_app_EMG(j).EMG.rect(:))+50]);
        xlim([0 trial_data_app_EMG(j).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['Musc: ',EMGparams.channelName{i}]);
    end
end

%% Fig per task, subplot per muscle (col) and target (row)
h = 0;
figure('Name','ForceCO');
for j = 1:length(Aparams.targetAnglesForce)
    for i = 1:length(EMGparams.channelName)-1
        h = h+1;
        subplot(length(Aparams.targetAnglesForce),length(EMGparams.channelName)-1,h);
        plot(trial_data_app_force(j).ts,trial_data_app_force(j).EMG.rect(:,i));
        hold on;
        plot(trial_data_app_force(j).ts,trial_data_app_force(j).EMG.avg(:,i));
        ylim([0 60]);%ylim([0 max(trial_data_app_force(j).EMG.rect(:))+50]);
        xlim([0 trial_data_app_force(j).ts(end)]);
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
        
        plot(trial_data_app_EMG(j).ts,trial_data_app_EMG(j).EMG.rect(:,i));
        hold on;
        plot(trial_data_app_EMG(j).ts,trial_data_app_EMG(j).EMG.avg(:,i));
        ylim([0 60]);%ylim([0 max(trial_data_app_EMG(j).EMG.rect(:))+50]);
        xlim([0 trial_data_app_EMG(j).ts(end)]);
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

%% EMG comparison between tasks. Fig per muscle, subplot per task
Aparams.angComp = {deg2rad([225 45]),deg2rad([45 135])};

for i = 1:length(EMGparams.channelName)-1
    figure('Name',EMGparams.channelName{i})
    h = 1;
    
    for j = 1:length(Aparams.angComp)
        subplot(length(Aparams.angComp),2,h)
        iangf = find([trial_data_app_force.angle] == Aparams.angComp{j}(1));
        plot(trial_data_app_force(iangf).ts,trial_data_app_force(iangf).EMG.rect(:,i))
        hold on
        plot(trial_data_app_force(iangf).ts,trial_data_app_force(iangf).EMG.avg(:,i),'r')
        xlim([0 trial_data_app_force(iangf).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['ForceCO; Target: ',num2str(rad2deg(Aparams.angComp{j}(1))),' deg']);
        
        subplot(length(Aparams.angComp),2,h+1)
        iangE = find([trial_data_app_EMG.angle] == Aparams.angComp{j}(2));
        plot(trial_data_app_EMG(iangE).ts,trial_data_app_EMG(iangE).EMG.rect(:,i))
        hold on
        plot(trial_data_app_EMG(iangE).ts,trial_data_app_EMG(iangE).EMG.avg(:,i),'r')
        xlim([0 trial_data_app_EMG(iangE).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['EMGCO; Target: ',num2str(rad2deg(Aparams.angComp{j}(2))),' deg (',EMGparams.channelName{EMGparams.channelControl(j)},')']);
        
        h = h+2;
    end
end

%% EMG comparison between tasks. Fig per angle, subplot per muscle for both tasks
for j = 1:length(Aparams.angComp)
    figure('Name',['Targets: Force-',num2str(rad2deg(Aparams.angComp{j}(1))),'; EMG-',num2str(rad2deg(Aparams.angComp{j}(2))),' (',EMGparams.channelName{EMGparams.channelControl(j)},')'])
    
    for i = 1:length(EMGparams.channelName)-1
        subplot(2,length(EMGparams.channelName)-1,i)
        iangf = find([trial_data_app_force.angle] == Aparams.angComp{j}(1));
        plot(trial_data_app_force(iangf).ts,trial_data_app_force(iangf).EMG.rect(:,i))
        hold on
        plot(trial_data_app_force(iangf).ts,trial_data_app_force(iangf).EMG.avg(:,i),'r')
        ylim([0 60]);%ylim([0 max(trial_data_app_force(iangf).EMG.rect(:))+10]);
        xlim([0 trial_data_app_force(iangf).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title([EMGparams.channelName{i}]);
        
        subplot(2,length(EMGparams.channelName)-1,i+length(EMGparams.channelName)-1)
        iangE = find([trial_data_app_EMG.angle] == Aparams.angComp{j}(2));
        plot(trial_data_app_EMG(iangE).ts,trial_data_app_EMG(iangE).EMG.rect(:,i))
        hold on
        plot(trial_data_app_EMG(iangE).ts,trial_data_app_EMG(iangE).EMG.avg(:,i),'r')
        ylim([0 60]);%ylim([0 max(trial_data_app_EMG(iangE).EMG.rect(:))+10]);
        xlim([0 trial_data_app_EMG(iangE).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title([EMGparams.channelName{i}]);
        
    end
end

%% EMG var comparison between tasks. Fig per angle, subplot per muscle for both tasks
figure('Name','EMG var')
for j = 1:length(Aparams.angComp)
    for i = 1:length(EMGparams.channelName)-1
        subplot(length(Aparams.angComp),1,j)    
        
        iangf = find([trial_data_app_force.angle] == Aparams.angComp{j}(1));
        stem(i,var(trial_data_app_force(iangf).EMG.rect(:,i)),'b')
        hold on
        iangE = find([trial_data_app_EMG.angle] == Aparams.angComp{j}(2));
        stem(i,var(trial_data_app_EMG(iangE).EMG.rect(:,i)),'r')
        
        xlim([0 length(EMGparams.channelName)]);
        ylabel('var [-]');
        title(['Targets: Force-',num2str(rad2deg(Aparams.angComp{j}(1))),'; EMG-',num2str(rad2deg(Aparams.angComp{j}(2))),' (',EMGparams.channelName{EMGparams.channelControl(j)},')']);
        legend('ForceCO','EMGCO')
    end
    xticklabels([{''},EMGparams.channelName(1:end-1)])
end

%% EMG SNR comparison between tasks. Fig per angle, subplot per muscle for both tasks
figure('Name','EMG SNR')
for j = 1:length(Aparams.angComp)
    for i = 1:length(EMGparams.channelName)-1
        subplot(length(Aparams.angComp),1,j)    
        
        iangf = find([trial_data_app_force.angle] == Aparams.angComp{j}(1));
        stem(i,mean(trial_data_app_force(iangf).EMG.rect(:,i))./std(trial_data_app_force(iangf).EMG.rect(:,i)),'b')
        hold on
        iangE = find([trial_data_app_EMG.angle] == Aparams.angComp{j}(2));
        stem(i,mean(trial_data_app_EMG(iangE).EMG.rect(:,i))./std(trial_data_app_EMG(iangE).EMG.rect(:,i)),'r')
        
        xlim([0 length(EMGparams.channelName)]);
        ylabel('SNR [-]');
        title(['Targets: Force-',num2str(rad2deg(Aparams.angComp{j}(1))),'; EMG-',num2str(rad2deg(Aparams.angComp{j}(2))),' (',EMGparams.channelName{EMGparams.channelControl(j)},')']);
        legend('ForceCO','EMGCO')
    end
    xticklabels([{''},EMGparams.channelName(1:end-1)])
end

%% Frequency analysis
Aparams.cohparams.data = 'app';
Aparams.cohparams.tseg = 1;
Aparams.cohparams.nseg = 10;
Aparams.cohparams.my_nseg = 10;
Aparams.cohparams.window = @(N) hanning(N);
field = 'rect';

if strcmp(Aparams.cohparams.data,'avg')
    trial_data_coh_force = cohStruct(trial_data_avg_force,EMGparams.channelName,{field},Aparams.cohparams);
    trial_data_coh_EMG = cohStruct(trial_data_avg_EMG,EMGparams.channelName,{field},Aparams.cohparams);
else
    trial_data_coh_force = cohStruct(trial_data_app_force,EMGparams.channelName,{field},Aparams.cohparams);
    trial_data_coh_EMG = cohStruct(trial_data_app_EMG,EMGparams.channelName,{field},Aparams.cohparams);
end

nmusc = length(EMGparams.channelName)-1;
musccont = EMGparams.channelControl;
nmusccomb = (length(EMGparams.channelName)-1)*(length(EMGparams.channelName)-2)/2; % n*(n-1)/2
fc = 80;

data_analysis = trial_data_app_EMG;
trial_data_coh = trial_data_coh_EMG;
nangles = Aparams.targetAnglesEMG; 

%% Coherence fig for each muscle combination, subplot for each angle
for j = 1:nmusccomb
    figure;
    set(gcf,'Name','Coherence');
    for i = 1:length(nangles)
        if rem(length(nangles),2) == 0
            subplot(2,length(nangles)/2,i);
        else
            subplot(1,length(nangles),i);
        end
        plot(trial_data_coh(i).(field).fcoh(:,j),trial_data_coh(i).(field).coh(:,j));
        hold on;
        line(xlim,trial_data_coh(i).(field).CL(j)*[1 1]);
        xlim([trial_data_coh(i).(field).fcoh(2,j) fc])
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_data_coh(i).(field).muscles{j}{1},',',trial_data_coh(i).(field).muscles{j}{2},'; Target: ',num2str(rad2deg(trial_data_coh(i).angle)),' deg']);
    end
end

%% Coherence for both tasks, same fig, only for angcomp
for j = 1:nmusccomb
    figure;
    set(gcf,'Name','Coherence');
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_coh_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_coh_EMG.angle] == Aparams.angComp{i}(2));
        
        if rem(length(Aparams.angComp),2) == 0
            subplot(2,length(Aparams.angComp)/2,i);
        else
            subplot(1,length(Aparams.angComp),i);
        end
        h1 = plot(trial_data_coh_force(iangf).(field).fcoh(:,j),trial_data_coh_force(iangf).(field).coh(:,j));
        hold on;
        line(xlim,trial_data_coh_force(iangf).(field).CL(j)*[1 1],'Color','k','LineStyle','--');
        
        h2 = plot(trial_data_coh_EMG(iangE).(field).fcoh(:,j),trial_data_coh_EMG(iangE).(field).coh(:,j),'r');
        hold on;
        line(xlim,trial_data_coh_EMG(iangE).(field).CL(j)*[1 1],'Color','k','LineStyle','--');
        
        xlim([trial_data_coh_force(iangf).(field).fcoh(2,j) fc])
        xlabel('Frequency [Hz]'); ylabel('Coh [-]');
        title(['Musc: ',trial_data_coh_force(iangf).(field).muscles{j}{1},...
            ',',trial_data_coh_force(iangf).(field).muscles{j}{2},...
            '; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',EMGparams.channelName{EMGparams.channelControl(i)},') deg']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end

%% FFT for both tasks, same fig, only for angcomp
for j = 1:length(EMGparams.channelName)-1
    figure;
    set(gcf,'Name','FFT');
    for i = 1:length(Aparams.angComp)
        iangf = find([trial_data_app_force.angle] == Aparams.angComp{i}(1));
        iangE = find([trial_data_app_EMG.angle] == Aparams.angComp{i}(2));
        
        if rem(length(Aparams.angComp),2) == 0
            subplot(2,length(Aparams.angComp)/2,i);
        else
            subplot(1,length(Aparams.angComp),i);
        end
        h1 = plot(trial_data_app_force(iangf).fv,abs(trial_data_app_force(iangf).EMG.rect_fft(:,j)));
        hold on;        
        h2 = plot(trial_data_app_EMG(iangE).fv,abs(trial_data_app_EMG(iangE).EMG.rect_fft(:,j)),'r');
        
        xlim([trial_data_app_force(iangf).fv(2) fc])
        xlabel('Frequency [Hz]'); ylabel('FFT [-]');
        title(['Musc: ',EMGparams.channelName{j},...
            '; Target: ',num2str(rad2deg(Aparams.angComp{i}(1))),'-',...
            num2str(rad2deg(Aparams.angComp{i}(2))),' (',EMGparams.channelName{EMGparams.channelControl(i)},') deg']);
        legend([h1,h2],'ForceCO','EMGCO')
    end
end