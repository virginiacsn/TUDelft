function[] = data2spss(trial_pp_force,trial_pp_EMG,Aparams_pp,angComp)

tasks = {'FC','MC'};
angles = angComp;
force = {'mag','CV'};
musc = {'TLH','DA','BB','DP','ECRB','FCR'};
muscpair = {'TLH-DA','DA-BB','BB-DP'};
bands = {'alp','bet','gam'};

vars_mean = {'Subject'};
vars_coh = {'Subject'};
data_mean = [(1:length(trial_pp_force))' nan(length(trial_pp_force),2*7*(2+6))];
data_coh = [(1:length(trial_pp_force))' nan(length(trial_pp_force),2*7*3*3)];
ivar = 1;

for imusc = 1:length(musc)
    for iang = 1:length(angles)
        for itask = 1:length(tasks)
            ivar = ivar+1;
            vars_mean{ivar} = [tasks{itask},'_',musc{imusc},'_',num2str(rad2deg(angles(iang)))];
            for isubject = 1:length(trial_pp_force)
                idm = find(contains(Aparams_pp(isubject).channelNameEMG,musc{imusc}));
                
                switch tasks{itask}
                    case 'FC'
                        iangf = find([trial_pp_force(isubject).angle] == angles(iang));
                        if ~isempty(iangf)
                            data_mean(isubject,ivar) = trial_pp_force(isubject).EMG.rect(iangf,idm)./Aparams_pp(isubject).EMGScaleForce(idm);
                        end
                    case 'MC'
                        iangE = find([trial_pp_EMG(isubject).angle] == angles(iang));
                        if ~isempty(iangE)
                            data_mean(isubject,ivar) = trial_pp_EMG(isubject).EMG.rect(iangE,idm)./Aparams_pp(isubject).EMGScaleForce(idm);
                        end
                end
                
            end
        end
    end
end
for iforce = 1:length(force)
    if strcmp(force{iforce},'mag')
        field = 'mag_mean';
    else
        field = 'mag_CV';
    end
    
    for iang = 1:length(angles)
        for itask = 1:length(tasks)
            ivar = ivar+1;
            vars_mean{ivar} = [tasks{itask},'_',force{iforce},'_',num2str(rad2deg(angles(iang)))];
            for isubject = 1:length(trial_pp_force)
                
                switch tasks{itask}
                    case 'FC'
                        iangf = find([trial_pp_force(isubject).angle] == angles(iang));
                        if ~isempty(iangf)
                            data_mean(isubject,ivar) = trial_pp_force(isubject).force.(field)(iangf);
                        end
                    case 'MC'
                        iangE = find([trial_pp_EMG(isubject).angle] == angles(iang));
                        if ~isempty(iangE)
                            data_mean(isubject,ivar) = trial_pp_EMG(isubject).force.(field)(iangE);
                        end
                end
            end
        end
    end
end

save4spss(vars_mean,data_mean, 'SPSS/stat_mean');

ivar = 1;
for imuscp = 1:length(muscpair)
    mp = strsplit(muscpair{imuscp},'-');
    for iang = 1:length(angles)
        for iband = 1:length(bands)
            for itask = 1:length(tasks)
                ivar = ivar+1;
                vars_coh{ivar} = [tasks{itask},'_',muscpair{imuscp},'_',num2str(rad2deg(angles(iang))),'_',bands{iband}];
                for isubject = 1:length(trial_pp_force)
                    switch tasks{itask}
                        case 'FC'
                            iangf = find([trial_pp_force(isubject).trial_coh.angle] == angles(iang));
                            if ~isempty(iangf)
                                mp1 = sum(reshape(contains([trial_pp_force(isubject).trial_coh(iangf).rect.muscles{:}],mp{1}),[2,length(trial_pp_force(isubject).trial_coh(iangf).rect.muscles)]));
                                mp2 = sum(reshape(contains([trial_pp_force(isubject).trial_coh(iangf).rect.muscles{:}],mp{2}),[2,length(trial_pp_force(isubject).trial_coh(iangf).rect.muscles)]));
                                idmp = find(mp1&mp2);
                                data_coh(isubject,ivar) = trial_pp_force(isubject).trial_coh(iangf).rect.nasig_z(iband,idmp);
                            end
                        case 'MC'
                            iangE = find([trial_pp_EMG(isubject).trial_coh.angle] == angles(iang));
                            if ~isempty(iangE)
                                mp1 = sum(reshape(contains([trial_pp_EMG(isubject).trial_coh(iangE).rect.muscles{:}],mp{1}),[2,length(trial_pp_EMG(isubject).trial_coh(iangE).rect.muscles)]));
                                mp2 = sum(reshape(contains([trial_pp_EMG(isubject).trial_coh(iangE).rect.muscles{:}],mp{2}),[2,length(trial_pp_EMG(isubject).trial_coh(iangE).rect.muscles)]));
                                idmp = find(mp1&mp2);
                                if ~isempty(idmp)
                                    data_coh(isubject,ivar) = trial_pp_EMG(isubject).trial_coh(iangE).rect.nasig_z(iband,idmp);
                                end
                            end
                    end
                    
                end
            end
        end
    end
end

save4spss(vars_coh,data_coh, 'SPSS/stat_coh');

end