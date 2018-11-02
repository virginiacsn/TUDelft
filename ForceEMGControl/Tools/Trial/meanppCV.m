function[CV] = meanppCV(trial_pp_force,trial_pp_EMG,angComp)

CV.angle = angComp;

for j = 1:length(angComp)
    CVforce = [];
    CVEMG = [];
    scountf = 0;
    scounte = 0;
    
    for i = 1:length(trial_pp_force)
        iangf = find([trial_pp_force(i).angle] == angComp(j));
        iangE = find([trial_pp_EMG(i).angle] == angComp(j));
        
        if ~isempty(iangf)
            scountf = scountf+1;
            CVforce = [CVforce,100*trial_pp_force(i).force.mag_pstd(iangf)./trial_pp_force(i).force.mag_mean(iangf)];
        end
        if ~isempty(iangE)
            scounte = scounte+1;
            CVEMG = [CVEMG,100*trial_pp_EMG(i).force.mag_pstd(iangE)./trial_pp_EMG(i).force.mag_mean(iangE)];
        end
    end
    CV.force.mean(j) = mean(CVforce);
    CV.force.std(j) = std(CVforce);
    CV.force.sem(j) = std(CVforce)/length(scountf);
    CV.EMG.mean(j) = mean(CVEMG);
    CV.EMG.std(j) = std(CVEMG);
    CV.EMG.sem(j) = std(CVEMG)/length(scounte);
end
end