function[trial_pp_force,trial_pp_EMG] = extraCoh(trial_pp_force,trial_pp_EMG)

for isubject = 1:length(trial_pp_force)
    trial_force_coh = trial_pp_force(isubject).trial_coh;
    trial_EMG_coh = trial_pp_EMG(isubject).trial_coh;
    
    for iangf = 1:length(trial_force_coh)
        for imc = 1:length(trial_force_coh(iangf).rect.muscles)
            
            freq = find(trial_force_coh(iangf).rect.my_fcoh(:,imc)>=8 & trial_force_coh(iangf).rect.my_fcoh(:,imc)<=20);
            z_temp_area = zeros(size(trial_force_coh(iangf).rect.z(:,imc)));
            idx_sig_z = trial_force_coh(iangf).rect.z(:,imc) >= 1.96;
            z_temp_area(idx_sig_z) = trial_force_coh(iangf).rect.z(idx_sig_z,imc);
            
            trial_pp_force(isubject).trial_coh(iangf).rect.nasig_z_band(imc) = trapz(trial_force_coh(iangf).rect.my_fcoh(freq,imc), z_temp_area(freq));
            %trial_pp_force(isubject).trial_coh(iangf).rect.nasig_z_band =
        end
        
        BB = sum(reshape(contains([trial_force_coh(iangf).rect.muscles{:}],'BB'),[2,length(trial_force_coh(iangf).rect.muscles)]));
        TLH = sum(reshape(contains([trial_force_coh(iangf).rect.muscles{:}],'TLH'),[2,length(trial_force_coh(iangf).rect.muscles)]));
        DA = sum(reshape(contains([trial_force_coh(iangf).rect.muscles{:}],'DA'),[2,length(trial_force_coh(iangf).rect.muscles)]));
        DP = sum(reshape(contains([trial_force_coh(iangf).rect.muscles{:}],'DP'),[2,length(trial_force_coh(iangf).rect.muscles)]));
        
        muscpair = find((BB&DA)|(DP&BB)|(DA&TLH));
        
        trial_pp_force(isubject).trial_coh(iangf).rect.pool_coh = abs(mean(trial_force_coh(iangf).rect.Syx(:,muscpair),2)).^2....
            ./(mean(trial_force_coh(iangf).rect.Sxx(:,muscpair),2).*mean(trial_force_coh(iangf).rect.Syy(:,muscpair),2));
        trial_pp_force(isubject).trial_coh(iangf).rect.pool_z = sqrt(2*trial_force_coh(iangf).rect.my_nseg(:,imc))*atanh(trial_pp_force(isubject).trial_coh(iangf).rect.pool_coh);
    end
    
    for iangE = 1:length(trial_EMG_coh)
        for imc = 1:length(trial_EMG_coh(iangE).rect.muscles)
            
            freq = find(trial_EMG_coh(iangE).rect.my_fcoh(:,imc)>=8 & trial_EMG_coh(iangE).rect.my_fcoh(:,imc)<=20);
            z_temp_area = zeros(size(trial_EMG_coh(iangE).rect.z(:,imc)));
            idx_sig_z = trial_EMG_coh(iangE).rect.z(:,imc) >= 1.96;
            z_temp_area(idx_sig_z) = trial_EMG_coh(iangE).rect.z(idx_sig_z,imc);
            
            trial_pp_EMG(isubject).trial_coh(iangE).rect.nasig_z_band(imc) = trapz(trial_EMG_coh(iangE).rect.my_fcoh(freq,imc), z_temp_area(freq));
            %trial_pp_force(isubject).trial_coh(iangf).rect.nasig_z_band =
            %segm = segm+trial_EMG_coh.rect.my_nseg(:,imc);
        end
        
        BB = sum(reshape(contains([trial_EMG_coh(iangE).rect.muscles{:}],'BB'),[2,length(trial_EMG_coh(iangE).rect.muscles)]));
        TLH = sum(reshape(contains([trial_EMG_coh(iangE).rect.muscles{:}],'TLH'),[2,length(trial_EMG_coh(iangE).rect.muscles)]));
        DA = sum(reshape(contains([trial_EMG_coh(iangE).rect.muscles{:}],'DA'),[2,length(trial_EMG_coh(iangE).rect.muscles)]));
        DP = sum(reshape(contains([trial_EMG_coh(iangE).rect.muscles{:}],'DP'),[2,length(trial_EMG_coh(iangE).rect.muscles)]));
        
        muscpair = find((BB&DA)|(DP&BB)|(DA&TLH));
        
        trial_pp_EMG(isubject).trial_coh(iangE).rect.pool_coh = abs(mean(trial_EMG_coh(iangE).rect.Syx(:,muscpair),2)).^2....
            ./(mean(trial_EMG_coh(iangE).rect.Sxx(:,muscpair),2).*mean(trial_EMG_coh(iangE).rect.Syy(:,muscpair),2));
        trial_pp_EMG(isubject).trial_coh(iangE).rect.pool_z = sqrt(2*trial_EMG_coh(iangE).rect.my_nseg(:,imc)).*atanh(trial_pp_EMG(isubject).trial_coh(iangE).rect.pool_coh);
    end
end
end