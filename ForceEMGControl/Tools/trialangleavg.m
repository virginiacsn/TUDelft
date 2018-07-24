function[trial_data_avg] = trialangleavg(trial_data, epoch, fields)

angles = sort(unique(extractfield(trial_data,'angle')));
nsamp = min(extractfield(trial_data,epoch{2})-extractfield(trial_data,epoch{1}));

for iangle = 1:length(angles)
    angle_data = trial_data(find(extractfield(trial_data,'angle') == angles(iangle)));
    trial_data_avg(iangle).angle = angles(iangle);
    trial_data_avg(iangle).ntrials = length(angle_data);
    
    for ifield = 1:length(fields)
        field_data = [];
        field_data_time = [];
        field_col = size(angle_data(1).(fields{ifield}),2);
        for itrial = 1:length(angle_data)
            idx1 = angle_data(itrial).(epoch{1});
            idx2 = angle_data(itrial).(epoch{2});
            nsampextra = (idx2-idx1)-nsamp;
            sampv = idx1+round(nsampextra/2):idx2-round(nsampextra/2)-1;
            field_data(itrial,:) = mean(angle_data(itrial).(fields{ifield})(idx1:idx2,:),1);
            field_data_time = [field_data_time reshape(angle_data(itrial).(fields{ifield})(sampv,:),field_col*length(sampv),1)];
        end
        trial_data_avg(iangle).(fields{ifield}) = mean(field_data,1);
        trial_data_avg(iangle).([fields{ifield},'_time']) = reshape(mean(field_data_time,2),length(sampv),field_col);
    end
end
end

        
        