function[meanEMGerr] = meanEMGErr(trial_pp_force,trial_pp_EMG,Aparams_pp,errortype)

EMGErr = zeros(length(trial_pp_force),4);

for i = 1:length(trial_pp_force)
    channel = Aparams_pp(i).chanControl;
    if strcmp(errortype,'abs')
        EMGErr(i,:) = abs(trial_pp_force(i).EMG.rect(i,channel)-trial_pp_EMG(i).EMG.rect(i,channel));
    elseif strcmp(errortype,'rms')
        EMGErr(i,:) = sqrt((trial_pp_force(i).EMG.rect(i,channel)-trial_pp_EMG(i).EMG.rect(i,channel)).^2);
    end
end

meanEMGerr = EMGErr;
%meanEMGErr.std = std(EMGErr,2);

figure('Name','Mean EMG Difference');
h = plot(rad2deg(Aparams_pp(1).targetAnglesForce),meanEMGerr,'o-.');
for i = 1:length(h)
    h(i).MarkerFaceColor = h(i).Color;
end
xticks(rad2deg(Aparams_pp(1).targetAnglesForce));
ylim([0 max(meanEMGerr(:))+1]);
xlabel('Target [deg]'); ylabel('Error [-]')
legend(Aparams_pp(1).chanControlName)
title('Mean EMG error Force-EMG tasks')
set(gca,'FontSize',12);
end