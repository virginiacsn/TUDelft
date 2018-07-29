function[trial_data_avg] = trialAngleAvg(trial_data, epoch, fields)

angles = sort(unique(extractfield(trial_data,'angle')));
nsamp = min(extractfield(trial_data,epoch{2})-extractfield(trial_data,epoch{1}));

for iangle = 1:length(angles)
    angle_data = trial_data(find(extractfield(trial_data,'angle') == angles(iangle)));
    trial_data_avg(iangle).angle = angles(iangle);
    trial_data_avg(iangle).ntrials = length(angle_data);
    %trial_data_avg(iangle).ts = 
    
    for ifield = 1:length(fields)
        field_str = strsplit(fields{ifield},'.');
        field_col = size(angle_data(1).(field_str{1}).(field_str{2}),2);

        field_data_mean = [];
        field_data = [];
        
        for itrial = 1:length(angle_data)
            idx1 = angle_data(itrial).(epoch{1});
            idx2 = angle_data(itrial).(epoch{2});
            nsampextra = (idx2-idx1)-nsamp;
            sampv = idx1+round(nsampextra/2):idx2-round(nsampextra/2)-1;
            field_data_mean(itrial,:) = mean(angle_data(itrial).(field_str{1}).(field_str{2})(idx1:idx2,:),1);
            field_data = [field_data reshape(angle_data(itrial).(field_str{1}).(field_str{2})(sampv,:),field_col*length(sampv),1)];
        end
        trial_data_avg(iangle).ts = (0:nsamp-1)*angle_data(1).dt;
        trial_data_avg(iangle).fv = (0:nsamp-1)/trial_data_avg(iangle).ts(end);
        trial_data_avg(iangle).(field_str{1}).([field_str{2},'_mean']) = mean(field_data_mean,1);
        trial_data_avg(iangle).(field_str{1}).(field_str{2}) = reshape(mean(field_data,2),length(sampv),field_col);
        trial_data_avg(iangle).(field_str{1}).([field_str{2},'_fft']) = fft(trial_data_avg(iangle).(field_str{1}).(field_str{2}));
        
    end
end
end

        
        