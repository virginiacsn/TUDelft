function[errEMG] = meanEMGErr(trial_pp_force,trial_pp_EMG,Aparams_pp,angComp,errortype)

errEMG.angle = angComp';

for j = 1:length(angComp)
    EMGerr = [];
    for i = 1:length(trial_pp_force)
        channel = Aparams_pp(i).chanControl;
        iangf = find([trial_pp_force(i).angles] == angComp(j));
        iangE = find([trial_pp_EMG(i).angles] == angComp(j));
        if i == 1 || i == 2
            trial_force_scaleEMG = trial_pp_force(i).EMG.rect./repmat([Aparams_pp(i).EMGScale 1],[length(trial_pp_force(i).angles),1]);
            trial_EMG_scaleEMG = trial_pp_EMG(i).EMG.rect./repmat([Aparams_pp(i).EMGScale 1],[length(trial_pp_EMG(i).angles),1]);
        else
            trial_force_scaleEMG = trial_pp_force(i).EMG.rect./repmat(Aparams_pp(i).EMGScale,[length(trial_pp_force(i).angles),1]);
            trial_EMG_scaleEMG = trial_pp_EMG(i).EMG.rect./repmat(Aparams_pp(i).EMGScale,[length(trial_pp_EMG(i).angles),1]);
        end
        if ~isempty(iangf)&&~isempty(iangE)
            if strcmp(errortype,'abs')
                EMGerr = [EMGerr; abs(trial_force_scaleEMG(iangf,channel)-trial_EMG_scaleEMG(iangE,channel))];
            elseif strcmp(errortype,'rms')
                EMGerr = [EMGerr; sqrt((trial_force_scaleEMG(iangf,channel)-trial_EMG_scaleEMG(iangE,channel)).^2)];
            end
        end
    end
    
    errEMG.mean(j,:) = mean(EMGerr,1);
    errEMG.std(j,:) = std(EMGerr,1);
    errEMG.sem(j,:) = std(EMGerr,1)/sqrt(length(trial_pp_force));
    %meanEMGErr.std = std(EMGErr,2);
end

figure('Name','Mean EMG Difference');
h = plot(rad2deg(Aparams_pp(1).targetAnglesForce),errEMG.mean,'o-.');
hold on;
for i = 1:length(h)
    h(i).MarkerFaceColor = h(i).Color;
    errorbar(rad2deg(Aparams_pp(1).targetAnglesForce),errEMG.mean(:,i),errEMG.sem(:,i),'.','Color',h(i).Color)
end
xticks(rad2deg(Aparams_pp(1).targetAnglesForce));
ylim([0 1]);
xlabel('Target [deg]'); ylabel('Error [-]')
legend(Aparams_pp(1).chanControlName)
title('Mean EMG error Force-EMG tasks')
set(gca,'FontSize',12);

% % Mean EMG difference - Polar
% nmusc = 4; cols = {'r','b','g','m'};
% figure('Name','Force Magnitude Mean');
% yerr = [];
% for i = 1:nmusc
%     polar(0,1,'w'); hold on;
%     [xerr,yerr(:,i)] = getBound(errEMG.angle,errEMG.mean(:,i),errEMG.std(:,i),'cart');
%     f = polar(errEMG.angle,errEMG.mean(:,i),[cols{i},'-.o']);
%     f.MarkerFaceColor = cols{i};
%     hold on;
%     fill(xerr,yerr(:,i),cols{i},'edgecolor','none','facealpha',0.3);
%     % xticks(rad2deg(angComp)); %thetaticklabels(Aparams.muscComp);
%     % xticklabels(rad2deg(angComp));
% end
% polar(0,max(yerr(:)),'w');
% set(gca,'FontSize',12);
end