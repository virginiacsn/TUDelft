%% Check force and EMG 
% Force: subplot per target
figure('Name','Force in time');
for i = 1:length(PreAparams.targetAngles)
    if rem(length(PreAparams.targetAngles),2) == 0
        subplot(2,length(PreAparams.targetAngles)/2,i);
    else
        subplot(1,length(PreAparams.targetAngles),i);
    end
    for j = 1:2
        plot(trial_data_avg(i).ts,trial_data_avg(i).force.filt(:,j));
        hold on;
    end
    ylim(taskparams.targetForce*[-5 5]);
    xlabel('Time [s]'); ylabel('Force [N]');
    title(['Target: ',num2str(rad2deg(trial_data_avg(i).angle)),' deg']);
    legend('Fx','Fy')
end

% EMG: subplot per muscle (col) and target (row)
h = 0;
figure('Name','EMG in time');
for j = 1:length(PreAparams.targetAngles)
    for i = 1:length(EMGparams.channelName)-1
        h = h+1;
        subplot(length(PreAparams.targetAngles),length(EMGparams.channelName)-1,h);
        plot(trial_data_avg(j).ts,trial_data_avg(j).EMG.rect(:,i));
        hold on;
        plot(trial_data_avg(j).ts,trial_data_avg(j).EMG.avg(:,i));
        ylim([0 max(trial_data_avg(j).EMG.rect(:))+50]);
        xlim([0 trial_data_avg(j).ts(end)]);
        if j == length(PreAparams.targetAngles)
            xlabel('Time [s]');
        end
        if i == 1
            ylabel([num2str(rad2deg(PreAparams.targetAngles(j))) ,' deg']);
        end
        if j == 1
            title(['Musc: ',EMGparams.channelName{i}]);
        end
    end
end