function[CV] = meanppCV(trial_pp_force,trial_pp_EMG,angComp)

CV.angle = angComp;

for iang = 1:length(angComp)
    CVforce = [];
    CVEMG = [];
    scountf = 0;
    scounte = 0;
    
    for isubject = 1:length(trial_pp_force)
        iangf = find([trial_pp_force(isubject).angle] == angComp(iang));
        iangE = find([trial_pp_EMG(isubject).angle] == angComp(iang));
        
        if ~isempty(iangf)
            scountf = scountf+1;
            CVforce = [CVforce,100*trial_pp_force(isubject).force.mag_pstd(iangf)./trial_pp_force(isubject).force.mag_mean(iangf)];
        end
        if ~isempty(iangE)
            scounte = scounte+1;
            CVEMG = [CVEMG,100*trial_pp_EMG(isubject).force.mag_pstd(iangE)./trial_pp_EMG(isubject).force.mag_mean(iangE)];
        end
    end
    CV.force.mean(iang) = mean(CVforce);
    CV.force.std(iang) = std(CVforce);
    CV.force.sem(iang) = std(CVforce)/length(scountf);
    CV.EMG.mean(iang) = mean(CVEMG);
    CV.EMG.std(iang) = std(CVEMG);
    CV.EMG.sem(iang) = std(CVEMG)/length(scounte);
end
end