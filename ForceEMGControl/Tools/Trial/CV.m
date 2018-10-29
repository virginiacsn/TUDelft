function[CV] = CV(trial_pp_force,trial_pp_EMG)

for i = 1:length(trial_pp_force)
    CVforce(i,:) = trial_pp_force(i).filtmag_std./trial_pp_force(i).filtmag_mean;
    CVEMG(i,:) = trial_pp_EMG(i).filtmag_std./trial_pp_EMG(i).filtmag_mean;
end

CV.force.mean = mean(CVforce);
CV.EMG.mean = mean(CVEMG);
CV.force.std = std(CVforce);
CV.EMG.std = std(CVEMG);

end