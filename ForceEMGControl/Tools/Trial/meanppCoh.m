function[meanCoh] = meanppCoh(trial_pp_force,trial_pp_EMG,Aparams_pp,angComp)

meanCoh.angle = angComp';

fields = {'my_coh','z','asig_coh','asig_z','nasig_coh','nasig_z'};

for iang = 1:length(angComp)
    
%     fprintf('Angle: %d\n\n',rad2deg(angComp(iang)));
    
    for h = 1:length(fields)
        if strcmp(fields{h},'my_coh')||strcmp(fields{h},'z')
            force_coh = zeros(1024,1);
            EMG_coh = zeros(1024,1);
            force_count = zeros(1024,1);
            EMG_count = zeros(1024,1);
            meanCoh.force.my_fcoh = trial_pp_force(1).trial_coh(1).rect.my_fcoh(:,1);
            meanCoh.EMG.my_fcoh = trial_pp_EMG(1).trial_coh(1).rect.my_fcoh(:,1);
        else
            force_coh = [];
            EMG_coh = [];
        end
        scountf = 0;
        scounte = 0;
        
        for isubject = 1:length(trial_pp_force)
            trial_pp_force_coh = trial_pp_force(isubject).trial_coh;
            trial_pp_EMG_coh = trial_pp_EMG(isubject).trial_coh;
            
            iangf = find([trial_pp_force_coh.angle] == angComp(iang));
            iangE = find([trial_pp_EMG_coh.angle] == angComp(iang));
            
            if ~isempty(iangf)
                BB = sum(reshape(contains([trial_pp_force_coh(iangf).rect.muscles{:}],'BB'),[2,length(trial_pp_force_coh(iangf).rect.muscles)]));
                TLH = sum(reshape(contains([trial_pp_force_coh(iangf).rect.muscles{:}],'TLH'),[2,length(trial_pp_force_coh(iangf).rect.muscles)]));
                DA = sum(reshape(contains([trial_pp_force_coh(iangf).rect.muscles{:}],'DA'),[2,length(trial_pp_force_coh(iangf).rect.muscles)]));
                DP = sum(reshape(contains([trial_pp_force_coh(iangf).rect.muscles{:}],'DP'),[2,length(trial_pp_force_coh(iangf).rect.muscles)]));
            elseif ~isempty(iangE)
                BB = sum(reshape(contains([trial_pp_EMG_coh(iangE).rect.muscles{:}],'BB'),[2,length(trial_pp_EMG_coh(iangE).rect.muscles)]));
                TLH = sum(reshape(contains([trial_pp_EMG_coh(iangE).rect.muscles{:}],'TLH'),[2,length(trial_pp_EMG_coh(iangE).rect.muscles)]));
                DA = sum(reshape(contains([trial_pp_EMG_coh(iangE).rect.muscles{:}],'DA'),[2,length(trial_pp_EMG_coh(iangE).rect.muscles)]));
                DP = sum(reshape(contains([trial_pp_EMG_coh(iangE).rect.muscles{:}],'DP'),[2,length(trial_pp_EMG_coh(iangE).rect.muscles)]));
            else
                BB = [];
                TLH = [];
                DA = [];
                DP = [];
            end
            % Find indexes of muscle pair combinations to plot. Combinations of control
            % muscles.
            % musccomb = find((plotBB&plotDA)|(plotDP&plotBB)|(plotDA&plotTLH)|(plotDA&plotDP)|(plotBB&plotTLH)|(plotTLH&plotDP));
            musccomb = find((BB&DA)|(DP&BB)|(DA&TLH));
            nmusccomb = length(musccomb);
            
            if ~isempty(iangf)
                scountf = scountf+1;
                if strcmp(fields{h},'my_coh')
                    force_coh = force_coh+trial_pp_force_coh(iangf).rect.(fields{h})(:,musccomb);
                    force_count = force_count+(trial_pp_force_coh(iangf).rect.(fields{h})(:,musccomb)>=trial_pp_force_coh(iangf).rect.my_CL(musccomb));
                elseif strcmp(fields{h},'z')
                     force_coh = force_coh+trial_pp_force_coh(iangf).rect.(fields{h})(:,musccomb);
                     force_count = force_count+(trial_pp_force_coh(iangf).rect.(fields{h})(:,musccomb)>=1.65);
                else
%                     if strcmp(fields{h},'nasig_z')
%                         fprintf('s%d. FC: %s-%s ',isubject,trial_pp_force_coh(iangf).rect.muscles{musccomb(1)}{1},trial_pp_force_coh(iangf).rect.muscles{musccomb(1)}{2})
%                         fprintf('%s-%s ',trial_pp_force_coh(iangf).rect.muscles{musccomb(2)}{1},trial_pp_force_coh(iangf).rect.muscles{musccomb(2)}{2})
%                         fprintf('%s-%s\n',trial_pp_force_coh(iangf).rect.muscles{musccomb(3)}{1},trial_pp_force_coh(iangf).rect.muscles{musccomb(3)}{2})
%                         
%                         fprintf('Alpha: %1.1f %1.1f %1.1f\n',trial_pp_force_coh(iangf).rect.(fields{h})(1,musccomb))
%                         fprintf('Beta: %1.1f %1.1f %1.1f\n',trial_pp_force_coh(iangf).rect.(fields{h})(2,musccomb))
%                         fprintf('Gamma: %1.1f %1.1f %1.1f\n\n',trial_pp_force_coh(iangf).rect.(fields{h})(3,musccomb))
%                     end
                    force_coh = [force_coh reshape(trial_pp_force_coh(iangf).rect.(fields{h})(:,musccomb),nmusccomb*3,1)];
                end
            end
            if ~isempty(iangE)
                scounte = scounte+1;
                if strcmp(fields{h},'my_coh')
                    EMG_coh = EMG_coh+trial_pp_EMG_coh(iangE).rect.(fields{h})(:,musccomb);
                    EMG_count = EMG_count+((trial_pp_EMG_coh(iangE).rect.(fields{h})(:,musccomb)>=trial_pp_EMG_coh(iangE).rect.my_CL(musccomb)));
                elseif strcmp(fields{h},'z')
                    EMG_coh = EMG_coh+trial_pp_EMG_coh(iangE).rect.(fields{h})(:,musccomb);
                    EMG_count = EMG_count+(trial_pp_EMG_coh(iangE).rect.(fields{h})(:,musccomb)>=1.65);
                else
%                     if strcmp(fields{h},'nasig_z')
%                         fprintf('s%d. MC: %s-%s ',isubject,trial_pp_EMG_coh(iangE).rect.muscles{musccomb(1)}{1},trial_pp_EMG_coh(iangE).rect.muscles{musccomb(1)}{2})
%                         fprintf('%s-%s ',trial_pp_EMG_coh(iangE).rect.muscles{musccomb(2)}{1},trial_pp_EMG_coh(iangE).rect.muscles{musccomb(2)}{2})
%                         fprintf('%s-%s\n',trial_pp_EMG_coh(iangE).rect.muscles{musccomb(3)}{1},trial_pp_EMG_coh(iangE).rect.muscles{musccomb(3)}{2})
%                         
%                         fprintf('Alpha: %1.1f %1.1f %1.1f\n',trial_pp_EMG_coh(iangE).rect.(fields{h})(1,musccomb))
%                         fprintf('Beta: %1.1f %1.1f %1.1f\n',trial_pp_EMG_coh(iangE).rect.(fields{h})(2,musccomb))
%                         fprintf('Gamma: %1.1f %1.1f %1.1f\n\n',trial_pp_EMG_coh(iangE).rect.(fields{h})(3,musccomb))
%                     end
                    EMG_coh = [EMG_coh reshape(trial_pp_EMG_coh(iangE).rect.(fields{h})(:,musccomb),nmusccomb*3,1)];
                end
            end
        end
        
        if strcmp(fields{h},'my_coh')||strcmp(fields{h},'z')
            meanCoh.force.(fields{h})(iang).scount = force_count;
            meanCoh.EMG.(fields{h})(iang).scount = EMG_count;

            meanCoh.force.(fields{h})(iang).mean = force_coh./scountf;
            meanCoh.EMG.(fields{h})(iang).mean = EMG_coh./scounte;
        else
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
    end
    meanCoh.musc = trial_pp_force_coh(1).rect.muscles(musccomb);
    
    for n = 1:length(meanCoh.musc)
        meanCoh.muscComp{n}{1} = meanCoh.musc{n}{1};
        meanCoh.muscComp{n}{2} = [meanCoh.musc{n}{1},',',meanCoh.musc{n}{2}];
        meanCoh.muscComp{n}{3} = meanCoh.musc{n}{2};
    end
    
    for n = 1:length(meanCoh.musc)
        for m = 1:length(Aparams_pp(end-1).muscCompPair)
            idm1 = contains(Aparams_pp(end-1).muscCompPair{m},meanCoh.musc{n}(1));
            idm2 = contains(Aparams_pp(end-1).muscCompPair{m},meanCoh.musc{n}(2));
            if sum(idm1+idm2) == 4
                meanCoh.angmusc{n} = Aparams_pp(end-1).angCompPair{m};
            end
        end
    end
    
end
end