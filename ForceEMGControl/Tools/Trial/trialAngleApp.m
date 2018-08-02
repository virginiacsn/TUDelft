function[trial_data_app] = trialAngleApp(trial_data, epoch, fields)

angles = sort(unique(extractfield(trial_data,'angle')));
nsamp = min(extractfield(trial_data,epoch{2})-extractfield(trial_data,epoch{1}));

for iangle = 1:length(angles)
    angle_data = trial_data(find(extractfield(trial_data,'angle') == angles(iangle)));
    trial_data_app(iangle).angle = angles(iangle);
    trial_data_app(iangle).ntrials = length(angle_data);
    
    for ifield = 1:length(fields)
        field_str = strsplit(fields{ifield},'.');
        field_col = size(angle_data(1).(field_str{1}).(field_str{2}),2);

        field_data = [];
        
        for itrial = 1:length(angle_data)
            idx1 = angle_data(itrial).(epoch{1});
            idx2 = angle_data(itrial).(epoch{2});
            
            field_data = [field_data; angle_data(itrial).(field_str{1}).(field_str{2})(idx1:idx2,:)];
        end
        trial_data_app(iangle).ts = (0:size(field_data,1)-1)*angle_data(itrial).dt;
        trial_data_app(iangle).fv = (0:size(field_data,1)-1)/trial_data_app(iangle).ts(end);
        trial_data_app(iangle).(field_str{1}).(field_str{2}) = field_data;
        trial_data_app(iangle).(field_str{1}).([field_str{2},'_fft']) = fft(field_data);
    end
end
end