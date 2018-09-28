function[trial_data_app] = trialAngleApp(trial_data, epoch, fields, window)

angles = sort(unique(extractfield(trial_data,'angle')));

for iangle = 1:length(angles)
    angle_data = trial_data(find(extractfield(trial_data,'angle') == angles(iangle)));
    trial_data_app(iangle).angle = angles(iangle);
    trial_data_app(iangle).ntrials = length(angle_data);
    
    for ifield = 1:length(fields)
        field_str = strsplit(fields{ifield},'.');
        field_col = size(angle_data(1).(field_str{1}).(field_str{2}),2);

        field_data = [];
        app_data = [];
        
        for itrial = 1:length(angle_data)
            if length(epoch) == 2
                idx1 = angle_data(itrial).(epoch{1});
                idx2 = angle_data(itrial).(epoch{2});
            else  
                idx1 = angle_data(itrial).(epoch{1})+round(epoch{2}/angle_data(itrial).dt);
                idx2 = angle_data(itrial).(epoch{3})+round(epoch{4}/angle_data(itrial).dt);
            end
            
            if isempty(window)
                win = rectwin(length(idx1:idx2));
            else
                win = window(length(idx1:idx2));
            end
            
            field_data = [field_data; angle_data(itrial).(field_str{1}).(field_str{2})(idx1:idx2,:).*win];
            app_data = [app_data; 1; zeros(length(idx1:idx2)-1,1)];
        end
        trial_data_app(iangle).ts = (0:size(field_data,1)-1)*angle_data(itrial).dt;
        trial_data_app(iangle).fv = (0:size(field_data,1)-1)/trial_data_app(iangle).ts(end);
        trial_data_app(iangle).iapp = app_data;
        trial_data_app(iangle).(field_str{1}).(field_str{2}) = field_data;
        %trial_data_app(iangle).(field_str{1}).([field_str{2},'_fft']) = fft(field_data);
    end
end
end