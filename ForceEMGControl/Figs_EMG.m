%% EMG FIGURES
%% TIME analysis
%% Fig per muscle (all) and task. Subplot per target.
% Compare targets for each task and muscle.
for j = 1:length(EMGparams.channelName)-1
    figure('Name',['ForceCO; EMG ',EMGparams.channelName{j}]);
    for i = 1:length(Aparams.targetAnglesForce)
        if rem(length(Aparams.targetAnglesForce),2) == 0
            subplot(2,length(Aparams.targetAnglesForce) /2,i);
        else
            subplot(1,length(Aparams.targetAnglesForce) ,i);
        end
        plot(trial_avg_force(i).ts,trial_avg_force(i).EMG.rect(:,j));
        hold on;
        plot(trial_avg_force(i).ts,trial_avg_force(i).EMG.avg(:,j));
        ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_force(i).EMG.rect(:))+50]);
        xlim([0 trial_avg_force(i).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['Target: ',num2str(rad2deg(trial_avg_force(i).angle)),' deg']);
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
        plot(trial_avg_EMG(i).ts,trial_avg_EMG(i).EMG.rect(:,j));
        hold on;
        plot(trial_avg_EMG(i).ts,trial_avg_EMG(i).EMG.avg(:,j));
        ylim([0 EMGlim]); %ylim([0 EMG_lim(i,j)+10]);
        xlim([0 trial_avg_EMG(i).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['Target: ',num2str(rad2deg(trial_avg_EMG(i).angle)),' deg']);
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
        plot(trial_avg_force(j).ts,trial_avg_force(j).EMG.rect(:,i));
        hold on;
        plot(trial_avg_force(j).ts,trial_avg_force(j).EMG.avg(:,i));
        ylim([0 EMG_lim(j,i)+5]);%ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_force(j).EMG.rect(:))+50]);
        xlim([0 trial_avg_force(j).ts(end)]);
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
        plot(trial_avg_EMG(j).ts,trial_avg_EMG(j).EMG.rect(:,i));
        hold on;
        plot(trial_avg_EMG(j).ts,trial_avg_EMG(j).EMG.avg(:,i));
        ylim([0 EMG_lim(j,i)+5]);%ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_EMG(j).EMG.rect(:))+50]);
        xlim([0 trial_avg_EMG(j).ts(end)]);
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
        plot(trial_avg_force(j).ts,trial_avg_force(j).EMG.rect(:,i));
        hold on;
        plot(trial_avg_force(j).ts,trial_avg_force(j).EMG.avg(:,i));
        ylim([0 max(EMG_lim(:,i))]);%ylim([0 EMGlim])%ylim([0 max(trial_data_avg_force(j).EMG.rect(:))+10]);%
        xlim([0 trial_avg_force(j).ts(end)]);
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
        
        plot(trial_avg_EMG(j).ts,trial_avg_EMG(j).EMG.rect(:,i));
        hold on;
        plot(trial_avg_EMG(j).ts,trial_avg_EMG(j).EMG.avg(:,i));
        ylim([0 max(EMG_lim(:,i))]);%ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_EMG(j).EMG.rect(:))+50]);
        xlim([0 trial_avg_EMG(j).ts(end)]);
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
        plot(trial_avg_force(j).ts,trial_avg_force(j).EMG.rectScale(:,Aparams.chanControl(i)));
        hold on;
        plot(trial_avg_force(j).ts,trial_avg_force(j).EMG.avgScale(:,Aparams.chanControl(i)),'LineWidth',2);
        ylim([0 2]);%ylim([0 EMG_lim(j,Aparams.chanControl(i))]); %ylim([0 max(trial_data_avg_force(j).EMG.rect(:))+10]);%
        xlim([0 trial_avg_force(j).ts(end)]);
        if j == length(Aparams.targetAnglesForce)
            xlabel('Time [s]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.targetAnglesForce(j)))],'FontWeight','bold');
        end
        if j == 1
            title([Aparams.chanControlName{i}]);
        end
        set(gca,'FontSize',12);
    end
end

h = 0;
figure('Name','EMGCO');
for j = 1:length(Aparams.targetAnglesEMG)
    for i = 1:length(Aparams.chanControl)
        h = h+1;
        subplot(length(Aparams.targetAnglesEMG),length(Aparams.chanControl),h);
        
        plot(trial_avg_EMG(j).ts,trial_avg_EMG(j).EMG.rectScale(:,Aparams.chanControl(i)));
        hold on;
        plot(trial_avg_EMG(j).ts,trial_avg_EMG(j).EMG.avgScale(:,Aparams.chanControl(i)),'LineWidth',2);
        ylim([0 2]);%ylim([0 EMG_lim(j,Aparams.chanControl(i))]);%ylim([0 max(trial_data_avg_EMG(j).EMG.rect(:))+50]);
        xlim([0 trial_avg_EMG(j).ts(end)]);
        if j == length(Aparams.targetAnglesEMG)
            xlabel('Time [s]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.targetAnglesEMG(j)))],'FontWeight','bold');
        end
        if j == 1
            title([Aparams.chanControlName{i}]);
        end
        set(gca,'FontSize',12);
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
        iangf = find([trial_avg_force.angle] == Aparams.angComp{j}(1));
        plot(trial_avg_force(iangf).ts,trial_avg_force(iangf).EMG.rect(:,i))
        hold on
        plot(trial_avg_force(iangf).ts,trial_avg_force(iangf).EMG.avg(:,i),'r')
        xlim([0 trial_avg_force(iangf).ts(end)]);
        if EMG_lim(j,i)>0
            ylim([0 EMG_lim(j,i)]);
        else
            ylim([0 EMGlim])
        end
        xlabel('Time [s]'); ylabel('EMG [-]');
        title(['ForceCO; Target: ',num2str(rad2deg(Aparams.angComp{j}(1))),' deg']);
        
        subplot(length(Aparams.angComp),2,h+1)
        iangE = find([trial_avg_EMG.angle] == Aparams.angComp{j}(2));
        plot(trial_avg_EMG(iangE).ts,trial_avg_EMG(iangE).EMG.rect(:,i))
        hold on
        plot(trial_avg_EMG(iangE).ts,trial_avg_EMG(iangE).EMG.avg(:,i),'r')
        xlim([0 trial_avg_EMG(iangE).ts(end)]);
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
        iangf = find([trial_avg_force.angle] == Aparams.angComp{j}(1));
        plot(trial_avg_force(iangf).ts,trial_avg_force(iangf).EMG.rect(:,i))
        hold on
        plot(trial_avg_force(iangf).ts,trial_avg_force(iangf).EMG.avg(:,i),'r')
        ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_force(iangf).EMG.rect(:))+10]);
        xlim([0 trial_avg_force(iangf).ts(end)]);
        xlabel('Time [s]'); ylabel('EMG [-]');
        title([EMGparams.channelName{i}]);
        
        subplot(2,length(EMGparams.channelName)-1,i+length(EMGparams.channelName)-1)
        iangE = find([trial_avg_EMG.angle] == Aparams.angComp{j}(2));
        plot(trial_avg_EMG(iangE).ts,trial_avg_EMG(iangE).EMG.rect(:,i))
        hold on
        plot(trial_avg_EMG(iangE).ts,trial_avg_EMG(iangE).EMG.avg(:,i),'r')
        ylim([0 EMGlim]);%ylim([0 max(trial_data_avg_EMG(iangE).EMG.rect(:))+10]);
        xlim([0 trial_avg_EMG(iangE).ts(end)]);
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
        
        iangf = find([trial_avg_force.angle] == Aparams.angComp{j}(1));
        stem(i,mean(trial_avg_force(iangf).EMG.rect(:,i)),'b')
        hold on
        iangE = find([trial_avg_EMG.angle] == Aparams.angComp{j}(2));
        stem(i,mean(trial_avg_EMG(iangE).EMG.rect(:,i)),'r')
        
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
        
        iangf = find([trial_avg_force.angle] == Aparams.angComp{j}(1));
        stem(i,mean(trial_avg_force(iangf).EMG.rect(:,Aparams.chanControl(i))),'b')
        hold on
        iangE = find([trial_avg_EMG.angle] == Aparams.angComp{j}(2));
        stem(i,mean(trial_avg_EMG(iangE).EMG.rect(:,Aparams.chanControl(i))),'r')
        
        xlim([0 length(Aparams.chanControl)+1]);
        ylabel('Mean EMG [-]');
        title(['Targets: Force-',num2str(rad2deg(Aparams.angComp{j}(1))),'; EMG-',num2str(rad2deg(Aparams.angComp{j}(2))),' (',Aparams.muscComp{j},')']);
        
        if j == 1
            legend('Force-Control','EMG-Control')
        end
    end
    xticks([0:length(Aparams.chanControl)]);
    xticklabels([{''},Aparams.chanControlName]);  
end

%% Fig of scaled EMG mean of control muscles. Subplot per target from angComp. 
% Compare means of control muscles for each target.
for k = 1:length(Aparams.muscComp)
    xticklab{k} = {num2str(rad2deg(Aparams.angCompUni(k)));Aparams.muscComp(k)};   
end
figure('Name','EMG mean')
    for i = 1:length(Aparams.chanControl)
        subplot(length(Aparams.chanControl),1,i)
        
        for j = 1:length(Aparams.angComp)
            iangf = find([trial_avg_force.angle] == Aparams.angComp{j}(1));
            
            h1 = plot(Aparams.angComp{j}(1),mean(trial_avg_force(iangf).EMG.rectScale(:,Aparams.chanControl(i))),'bo','MarkerFaceColor','b');
            hold on
            errorbar(Aparams.angComp{j}(1),mean(trial_avg_force(iangf).EMG.rectScale(:,Aparams.chanControl(i))),...
                std(trial_avg_force(iangf).EMG.rectScale(:,Aparams.chanControl(i))),'b.')
            iangE = find([trial_avg_EMG.angle] == Aparams.angComp{j}(2));
            h2 = plot(Aparams.angComp{j}(2),mean(trial_avg_EMG(iangE).EMG.rectScale(:,Aparams.chanControl(i))),'ro','MarkerFaceColor','r');
            errorbar(Aparams.angComp{j}(2),mean(trial_avg_EMG(iangE).EMG.rectScale(:,Aparams.chanControl(i))),...
                std(trial_avg_EMG(iangE).EMG.rectScale(:,Aparams.chanControl(i))),'r.')
        end
        ylim([0 1.5]);xlim([Aparams.angCompUni(1)-0.2 Aparams.angCompUni(end)+0.2])
        ylabel('EMG [-]'); 
        xticklabels(rad2deg(Aparams.angCompUni));
        set(gca,'XTick',[Aparams.angCompUni])
        title(Aparams.chanControlName(i));
        if i == 1
            legend([h1,h2],'FC','MC')
        elseif i == length(Aparams.chanControl)
            xlabel('Target [deg]');
        end
        set(gca,'FontSize',12);
    end
%     xticklabels([{''},Aparams.chanControlName]);  

%% Fig of EMG variance of all muscles. Subplot per target from angComp.
figure('Name','EMG var')
for j = 1:length(Aparams.angComp)
    for i = 1:length(EMGparams.channelName)-1
        subplot(length(Aparams.angComp),1,j)
        
        iangf = find([trial_avg_force.angle] == Aparams.angComp{j}(1));
        stem(i,var(trial_avg_force(iangf).EMG.rect(:,i)),'b')
        hold on
        iangE = find([trial_avg_EMG.angle] == Aparams.angComp{j}(2));
        stem(i,var(trial_avg_EMG(iangE).EMG.rect(:,i)),'r')
        
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
        
        iangf = find([trial_avg_force.angle] == Aparams.angComp{j}(1));
        stem(i,mean(trial_avg_force(iangf).EMG.rect(:,i))./std(trial_avg_force(iangf).EMG.rect(:,i)),'b')
        hold on
        iangE = find([trial_avg_EMG.angle] == Aparams.angComp{j}(2));
        stem(i,mean(trial_avg_EMG(iangE).EMG.rect(:,i))./std(trial_avg_EMG(iangE).EMG.rect(:,i)),'r')
        
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
        iangf = find([trial_avg_force.angle] == Aparams.angComp{j}(1));
        plot(trial_avg_force(iangf).fv,abs(trial_avg_force(iangf).EMG.fftraw(:,i)));
        hold on;
        plot(trial_avg_force(iangf).fv,abs(trial_avg_force(iangf).EMG.fftfilt(:,i)));
        plot(trial_avg_force(iangf).fv,abs(trial_avg_force(iangf).EMG.fftrect(:,i)));
        xlim([1 fc]); ylim([0 max([fft_lim.raw(j,i),fft_lim.filt(j,i),fft_lim.rect(j,i)])]);
        if j == length(Aparams.angComp)
            xlabel('Frequency [Hz]');
        end
        ylabel('a.u. [-]');
        legend('Raw','Filt','Rect')
        title(['ForceCO; Target: ',num2str(rad2deg(Aparams.angComp{j}(1))),' deg']);
        
        subplot(length(Aparams.angComp),2,h+1)
        iangE = find([trial_avg_EMG.angle] == Aparams.angComp{j}(2));
        plot(trial_avg_EMG(iangE).fv,abs(trial_avg_EMG(iangE).EMG.fftraw(:,i)));
        hold on;
        plot(trial_avg_EMG(iangE).fv,abs(trial_avg_EMG(iangE).EMG.fftfilt(:,i)));
        plot(trial_avg_EMG(iangE).fv,abs(trial_avg_EMG(iangE).EMG.fftrect(:,i)));
        xlim([1 fc]); ylim([0 max([fft_lim.raw(j,i),fft_lim.filt(j,i),fft_lim.rect(j,i)])]);
        if j == length(Aparams.angComp)
            xlabel('Frequency [Hz]');
        end
        ylabel('a.u. [-]');
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
        subplot(length(Aparams.targetAnglesForce),length(Aparams.chanControl),h);
        plot(trial_avg_force(j).fv,abs(trial_avg_force(j).EMG.fftraw(:,Aparams.chanControl(i))));
        hold on;
        plot(trial_avg_force(j).fv,abs(trial_avg_force(j).EMG.fftfilt(:,Aparams.chanControl(i))));
        plot(trial_avg_force(j).fv,abs(trial_avg_force(j).EMG.fftrect(:,Aparams.chanControl(i))));
        xlim([1 fc]); ylim([0 max([max(fft_lim.raw(:,Aparams.chanControl(i))),max(fft_lim.filt(:,Aparams.chanControl(i))),max(fft_lim.rect(:,Aparams.chanControl(i)))])]);
        legend('Raw','Filt','Rect');
        if j == length(Aparams.targetAnglesForce)
            xlabel('Frequency [Hz]');
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
        subplot(length(Aparams.targetAnglesEMG),length(Aparams.chanControl),h);
        plot(trial_avg_EMG(j).fv,abs(trial_avg_EMG(j).EMG.fftraw(:,Aparams.chanControl(i))));
        hold on;
        plot(trial_avg_EMG(j).fv,abs(trial_avg_EMG(j).EMG.fftfilt(:,Aparams.chanControl(i))));
        plot(trial_avg_EMG(j).fv,abs(trial_avg_EMG(j).EMG.fftrect(:,Aparams.chanControl(i))));
        xlim([1 fc]); ylim([0 max([max(fft_lim.raw(:,Aparams.chanControl(i))),max(fft_lim.filt(:,Aparams.chanControl(i))),max(fft_lim.rect(:,Aparams.chanControl(i)))])]);
        legend('Raw','Filt','Rect');
        if j == length(Aparams.targetAnglesEMG)
            xlabel('Frequency [Hz]');
        end
        if i == 1
            ylabel([num2str(rad2deg(Aparams.targetAnglesEMG(j))) ,' deg']);
        end
        if j == 1
            title(['Musc: ',Aparams.chanControlName{i}]);
        end
    end
end
