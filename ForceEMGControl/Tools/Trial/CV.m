function[sCV] = CV(trial_pp_force,trial_pp_EMG,angComp)

sCV.angle = angComp;

for j = 1:length(angComp)
    CVforce = [];
    CVEMG = [];
    for i = 1:length(trial_pp_force)
        iangf = find([trial_pp_force(i).angles] == angComp(j));
        iangE = find([trial_pp_EMG(i).angles] == angComp(j));
        
        if ~isempty(iangf)
            CVforce = [CVforce,100*trial_pp_force(i).force.filtmag_pstd(iangf)./trial_pp_force(i).force.filtmag_mean(iangf)];
        end
        if ~isempty(iangE)
            CVEMG = [CVEMG,100*trial_pp_EMG(i).force.filtmag_pstd(iangE)./trial_pp_EMG(i).force.filtmag_mean(iangE)];
        end
    end
    sCV.force.mean(j) = mean(CVforce);
    sCV.EMG.mean(j) = mean(CVEMG);
    sCV.force.std(j) = std(CVforce);
    sCV.EMG.std(j) = std(CVEMG);
    sCV.force.sem(j) = std(CVforce)/length(trial_pp_force);
    sCV.EMG.sem(j) = std(CVEMG)/length(trial_pp_force); 
end
end