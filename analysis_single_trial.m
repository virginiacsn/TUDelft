%% SINGLE-TRIAL
%% Single trial from trial_data
itrial = 213;
trial_data_ST = trial_data_EMG;
figure('Name',['Target: ', num2str(rad2deg(trial_data_ST(itrial).angle)), ' deg']);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(EMGparams.channelName)-1
    subplot(1,length(EMGparams.channelName)-1,i);
    plot(trial_data_ST(itrial).ts,trial_data_ST(itrial).EMG.rect(:,i));
    hold on;
    plot(trial_data_ST(itrial).ts,trial_data_ST(itrial).EMG.avg(:,i),'r');
    xlim([0 trial_data_ST(itrial).ts(end)]);
    %ylim([0 max(trial_data(itrial).EMG.rect(:,i))]);
    xlabel('Time [s]'); ylabel('EMG [-]');
    title([Aparams.channelNameEMG{i}]);
end

%% Compare tasks
load([filepath,'Parameters/',[date,'_s',subject,'_params_','001','.mat']])

forceAngles = [trial_data_force.angle];
EMGAngles = [trial_data_EMG.angle];

compAngles = Aparams.targetAnglesForce(ismember(Aparams.targetAnglesForce,Aparams.targetAnglesEMG));

for i = 1:length(compAngles)
    trialsAngleForce{i} = find(forceAngles == compAngles(i));
    trialsAngleEMG{i}  = find(EMGAngles == compAngles(i));
end

%% EMG
% Rectified and smoothed EMG
iAngle = 5;
for h = 8:10
    itrialForce = trialsAngleForce{iAngle}(h);
    itrialEMG = trialsAngleEMG{iAngle}(h);
    
    figure('Name',[' Target: ',num2str(rad2deg(trial_data_force(itrialForce).angle)),' deg']);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for j = 1:length(Aparams.channelNameEMG)-1
        subplot(2,length(Aparams.channelNameEMG)-1,j);
        plot(trial_data_force(itrialForce).ts,trial_data_force(itrialForce).EMG.rect(:,j));
        hold on;
        plot(trial_data_force(itrialForce).ts,trial_data_force(itrialForce).EMG.avg(:,j));
        xlabel('Time [s]'); ylabel('EMG [-]');
        xlim([0 trial_data_force(itrialForce).ts(end)]);
        ylim([0 max(trial_data_force(itrialForce).EMG.rect(:,j))]);
        title(['Musc: ',Aparams.channelNameEMG{j}]);
        
    end
    for j = 1:length(Aparams.channelNameEMG)-1
        subplot(2,length(Aparams.channelNameEMG)-1,j+length(Aparams.channelNameEMG)-1);
        plot(trial_data_EMG(itrialEMG).ts,trial_data_EMG(itrialEMG).EMG.rect(:,j));
        hold on;
        plot(trial_data_EMG(itrialEMG).ts,trial_data_EMG(itrialEMG).EMG.avg(:,j));
        xlabel('Time [s]'); ylabel('EMG [-]');
        xlim([0 trial_data_EMG(itrialEMG).ts(end)]);
        ylim([0 max(trial_data_EMG(itrialEMG).EMG.rect(:,j))]);
        title(['Musc: ',Aparams.channelNameEMG{j}]);
    end
end
