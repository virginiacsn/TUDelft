function[ppStats] = statAnalysis(trial_pp_force,trial_pp_EMG,Aparams_pp,angComp,fields)

ppStats.angle = angComp';
bands = {'alp','beta','gam'};
for iang = 1:length(angComp)
    for ifield = 1:length(fields)
        if any(strfind(fields{ifield},'.'))
            field_str = strsplit(fields{ifield},'.');
        else
            field_str = {fields{ifield}};
        end
        
        force_all = []; EMG_all = [];        
        force_alp = []; EMG_alp = [];                
        force_beta = []; EMG_beta = [];                
        force_gam = []; EMG_gam = [];
        
        scountf = 0;
        scounte = 0;
        
        for isubject = 1:length(trial_pp_force)
            if strcmp(field_str{1},'trial_coh')
                trial_force = trial_pp_force(isubject).trial_coh;
                trial_EMG = trial_pp_EMG(isubject).trial_coh;
                
                BB = sum(reshape(contains([trial_force(iangf).rect.muscles{:}],'BB'),[2,length(trial_force(iangf).rect.muscles)]));
                TLH = sum(reshape(contains([trial_force(iangf).rect.muscles{:}],'TLH'),[2,length(trial_force(iangf).rect.muscles)]));
                DA = sum(reshape(contains([trial_force(iangf).rect.muscles{:}],'DA'),[2,length(trial_force(iangf).rect.muscles)]));
                DP = sum(reshape(contains([trial_force(iangf).rect.muscles{:}],'DP'),[2,length(trial_force(iangf).rect.muscles)]));
                
                if ~isempty(iangf)
                    BB = sum(reshape(contains([trial_force(iangf).rect.muscles{:}],'BB'),[2,length(trial_force(iangf).rect.muscles)]));
                    TLH = sum(reshape(contains([trial_force(iangf).rect.muscles{:}],'TLH'),[2,length(trial_force(iangf).rect.muscles)]));
                    DA = sum(reshape(contains([trial_force(iangf).rect.muscles{:}],'DA'),[2,length(trial_force(iangf).rect.muscles)]));
                    DP = sum(reshape(contains([trial_force(iangf).rect.muscles{:}],'DP'),[2,length(trial_force(iangf).rect.muscles)]));
                elseif ~isempty(iangE)
                    BB = sum(reshape(contains([trial_EMG(iangE).rect.muscles{:}],'BB'),[2,length(trial_EMG(iangE).rect.muscles)]));
                    TLH = sum(reshape(contains([trial_EMG(iangE).rect.muscles{:}],'TLH'),[2,length(trial_EMG(iangE).rect.muscles)]));
                    DA = sum(reshape(contains([trial_EMG(iangE).rect.muscles{:}],'DA'),[2,length(trial_EMG(iangE).rect.muscles)]));
                    DP = sum(reshape(contains([trial_EMG(iangE).rect.muscles{:}],'DP'),[2,length(trial_EMG(iangE).rect.muscles)]));
                else
                    BB = [];
                    TLH = [];
                    DA = [];
                    DP = [];
                end
                
                musccomb = find((BB&DA)|(DP&BB)|(DA&TLH));
                nmusccomb = length(musccomb);
            else
                trial_force = trial_pp_force(isubject);
                trial_EMG = trial_pp_EMG(isubject);
            end
            
            iangf = find([trial_force.angle] == angComp(iang));
            iangE = find([trial_EMG.angle] == angComp(iang));
            
            if strcmp(field_str{1},'trial_coh')
                if ~isempty(iangf)
                    scountf = scountf+1;
                    force_alp = [force_alp; trial_force(iangf).(field_str{2}).(field_str{3})(1,musccomb)];
                    force_beta = [force_beta; trial_force(iangf).(field_str{2}).(field_str{3})(2,musccomb)];
                    force_gam = [force_gam; trial_force(iangf).(field_str{2}).(field_str{3})(3,musccomb)];
                end
                if ~isempty(iangE)
                    scounte = scounte+1;
                    EMG_alp = [EMG_alp; trial_EMG(iangE).(field_str{2}).(field_str{3})(1,musccomb)];
                    EMG_beta = [EMG_beta; trial_EMG(iangE).(field_str{2}).(field_str{3})(2,musccomb)];
                    EMG_gam = [EMG_gam; trial_EMG(iangE).(field_str{2}).(field_str{3})(3,musccomb)];
                end
            elseif strcmp(field_str{1},'EMG')
                channel = Aparams_pp(isubject).chanControl;
                
                trial_force_scaleEMG = trial_force.EMG.(field_str{2})./repmat(Aparams_pp(isubject).EMGScale,[length(trial_force.angle),1]);
                trial_EMG_scaleEMG = trial_EMG.EMG.(field_str{2})./repmat(Aparams_pp(isubject).EMGScale,[length(trial_EMG.angle),1]);
                
                force_all = [force_all; trial_force_scaleEMG(iangf,channel)];               
                EMG_all = [EMG_all; trial_EMG_scaleEMG(iangE,channel)];               
                
            else
                if ~isempty(iangf)
                    scountf = scountf+1;
                    force_all = [force_all; trial_force.(field_str{1}).(field_str{2})(iangf)];
                end
                if ~isempty(iangE)
                    scounte = scounte+1;
                    EMG_all = [EMG_all; trial_EMG.(field_str{1}).(field_str{2})(iangE)];
                end
                
            end
            
        end
        
        if strcmp(field_str{1},'trial_coh')
            ppStats.(field_str{3}).musc = trial_force(1).rect.muscles(musccomb);
            for n = 1:length(ppStats.(field_str{3}).musc)
                for m = 1:length(Aparams_pp(end).muscCompPair)
                    idm1 = contains(Aparams_pp(end).muscCompPair{m},ppStats.(field_str{3}).musc{n}(1));
                    idm2 = contains(Aparams_pp(end).muscCompPair{m},ppStats.(field_str{3}).musc{n}(2));
                    if sum(idm1+idm2) == 4
                        ppStats.(field_str{3}).angmusc{n} = Aparams_pp(end).angCompPair{m};
                    end
                end
            end
            for imusc = 1:nmusccomb
                [ppStats.(field_str{3}).alp.p(iang,imusc),ppStats.(field_str{3}).alp.tbl(iang,imusc)] = anova1([force_alp(:,imusc),EMG_alp(:,imusc)]);
                [ppStats.(field_str{3}).beta.p(iang,imusc),ppStats.(field_str{3}).beta.tbl(iang,imusc)] = anova1([force_beta(:,imusc),EMG_beta(:,imusc)]);
                [ppStats.(field_str{3}).gam.p(iang,imusc),ppStats.(field_str{3}).gam.tbl(iang,imusc)] = anova1([force_gam(:,imusc),EMG_gam(:,imusc)]);
                
%                 [ppStats.(field_str{3}).alp.p(iang,imusc),ppStats.(field_str{3}).alp.h(iang,imusc)] = ranksum(force_alp(:,imusc),EMG_alp(:,imusc),'alpha',0.1);
%                 [ppStats.(field_str{3}).beta.p(iang,imusc),ppStats.(field_str{3}).beta.h(iang,imusc)] = ranksum(force_beta(:,imusc),EMG_beta(:,imusc),'alpha',0.1);
%                 [ppStats.(field_str{3}).gam.p(iang,imusc),ppStats.(field_str{3}).gam.h(iang,imusc)] = ranksum(force_gam(:,imusc),EMG_gam(:,imusc),'alpha',0.1);
            end
        elseif strcmp(field_str{1},'EMG')
            for imusc = 1:length(channel)
                ppStats.(field_str{2}).musc = Aparams_pp(end).chanControlName;
                [ppStats.(field_str{2}).p(iang,imusc),ppStats.(field_str{2}).tbl(iang,imusc)] = anova1([force_all(:,imusc),EMG_all(:,imusc)]);
                
%                 [ppStats.(field_str{2}).p(iang,imusc),ppStats.(field_str{2}).h(iang,imusc)] = ranksum(force_all(:,imusc),EMG_all(:,imusc),'alpha',0.1);
            end
        else
            %[ppStats.(field_str{2}).p(iang),ppStats.(field_str{2}).tbl(iang)] = anova1([force_all,EMG_all]);
            
%             [ppStats.(field_str{2}).p(iang),ppStats.(field_str{2}).h(iang)] = ranksum(force_all,EMG_all,'alpha',0.1);
        end
    end
end
end