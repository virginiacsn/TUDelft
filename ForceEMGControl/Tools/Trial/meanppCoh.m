function[meanCoh] = meanppCoh(trial_pp_force,trial_pp_EMG,Aparams_pp,angComp)

meanCoh.angle = angComp';

fields = {'msig_coh','asig_coh','asig_z','nasig_coh','nasig_z'};

for iang = 1:length(angComp)
    for h = 1:length(fields)
        force_coh = [];
        EMG_coh = [];
        scountf = 0;
        scounte = 0;
        
        for isubject = 1:length(trial_pp_force)
            trial_pp_force_coh = trial_pp_force(isubject).trial_coh;
            trial_pp_EMG_coh = trial_pp_EMG(isubject).trial_coh;
            
            iangf = find([trial_pp_force_coh.angle] == angComp(iang));
            iangE = find([trial_pp_EMG_coh.angle] == angComp(iang));
            
            
            BB = sum(reshape(contains([trial_pp_force_coh(iangf).rect.muscles{:}],'BB'),[2,length(trial_pp_force_coh(iangf).rect.muscles)]));
            TLH = sum(reshape(contains([trial_pp_force_coh(iangf).rect.muscles{:}],'TLH'),[2,length(trial_pp_force_coh(iangf).rect.muscles)]));
            DA = sum(reshape(contains([trial_pp_force_coh(iangf).rect.muscles{:}],'DA'),[2,length(trial_pp_force_coh(iangf).rect.muscles)]));
            DP = sum(reshape(contains([trial_pp_force_coh(iangf).rect.muscles{:}],'DP'),[2,length(trial_pp_force_coh(iangf).rect.muscles)]));
            
            % Find indexes of muscle pair combinations to plot. Combinations of control
            % muscles.
            % musccomb = find((plotBB&plotDA)|(plotDP&plotBB)|(plotDA&plotTLH)|(plotDA&plotDP)|(plotBB&plotTLH)|(plotTLH&plotDP));
            musccomb = find((BB&DA)|(DP&BB)|(DA&TLH));
            nmusccomb = length(musccomb);
            
            if ~isempty(iangf)
                scountf = scountf+1;
                force_coh = [force_coh reshape(trial_pp_force_coh(iangf).rect.(fields{h})(:,musccomb),nmusccomb*3,1)];
            end
            if ~isempty(iangE)
                scounte = scounte+1;
                EMG_coh = [EMG_coh reshape(trial_pp_EMG_coh(iangE).rect.(fields{h})(:,musccomb),nmusccomb*3,1)];
            end
        end
        
        fmean = reshape(mean(force_coh,2),3,nmusccomb);
        fstd = reshape(std(force_coh,[],2),3,nmusccomb);
        fsem = reshape(std(force_coh,[],2)/sqrt(scountf),3,nmusccomb);
        emean = reshape(mean(EMG_coh,2),3,nmusccomb);
        estd = reshape(std(EMG_coh,[],2),3,nmusccomb);
        esem = reshape(std(EMG_coh,[],2)/sqrt(scounte),3,nmusccomb);
        
        bands = {'alp','beta','gam'};
        for k = 1:3
            meanCoh.force.(fields{h}).(bands{k}).mean(iang,:) = fmean(k,:);
            meanCoh.force.(fields{h}).(bands{k}).std(iang,:) = fstd(k,:);
            meanCoh.force.(fields{h}).(bands{k}).sem(iang,:) = fsem(k,:);
            meanCoh.EMG.(fields{h}).(bands{k}).mean(iang,:) = emean(k,:);
            meanCoh.EMG.(fields{h}).(bands{k}).std(iang,:) = estd(k,:);
            meanCoh.EMG.(fields{h}).(bands{k}).sem(iang,:) = esem(k,:);
        end
    end
    meanCoh.musc = trial_pp_force_coh(1).rect.muscles(musccomb);
    
    for n = 1:length(meanCoh.musc)
        for m = 1:length(Aparams_pp(end).muscCompPair)
            idm1 = contains(Aparams_pp(end).muscCompPair{m},meanCoh.musc{n}(1));
            idm2 = contains(Aparams_pp(end).muscCompPair{m},meanCoh.musc{n}(2));
            if sum(idm1+idm2) == 4
                meanCoh.angmusc{n} = Aparams_pp(end).angCompPair{m};
            end
        end
    end
end
end