function[meanEMGerr] = meanppEMGerr(trial_pp_force,trial_pp_EMG,Aparams_pp,angComp,errortype)

meanEMGerr.angle = angComp';

for iang = 1:length(angComp)
    EMGerr = [];
    scount = 0;
    
    for isubject = 1:length(trial_pp_force)
        channel = Aparams_pp(isubject).chanControl;
        iangf = find([trial_pp_force(isubject).angle] == angComp(iang));
        iangE = find([trial_pp_EMG(isubject).angle] == angComp(iang));
        
        trial_force_scaleEMG = trial_pp_force(isubject).EMG.rect./repmat(Aparams_pp(isubject).EMGScale,[length(trial_pp_force(isubject).angle),1]);
        trial_EMG_scaleEMG = trial_pp_EMG(isubject).EMG.rect./repmat(Aparams_pp(isubject).EMGScale,[length(trial_pp_EMG(isubject).angle),1]);
        
        if ~isempty(iangf)&&~isempty(iangE)
            scount = scount+1;
            if strcmp(errortype,'abs')
                EMGerr = [EMGerr; abs(trial_force_scaleEMG(iangf,channel)-trial_EMG_scaleEMG(iangE,channel))];
            elseif strcmp(errortype,'rms')
                EMGerr = [EMGerr; sqrt((trial_force_scaleEMG(iangf,channel)-trial_EMG_scaleEMG(iangE,channel)).^2)];
            end
        end
    end
    
    meanEMGerr.mean(iang,:) = mean(EMGerr,1);
    meanEMGerr.std(iang,:) = std(EMGerr,1);
    meanEMGerr.sem(iang,:) = std(EMGerr,1)/sqrt(scount);
    %meanEMGErr.std = std(EMGErr,2);
end

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