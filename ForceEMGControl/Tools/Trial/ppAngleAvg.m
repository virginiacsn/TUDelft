function[trial_pp_avg] = ppAngleAvg(trial_pp, fields, Aparams_pp)

angles = sort(unique(extractfield(trial_pp,'angles')));

for iangle = 1:length(angles)
    
    trial_pp_avg(iangle).angle = angles(iangle);
    
    for ifield = 1:length(fields)
        if any(strfind(fields{ifield},'.'))
            field_str = strsplit(fields{ifield},'.');
            field_col = size(trial_pp(1).(field_str{1}).(field_str{2}),2);
        else
            field_str = {fields{ifield}};
            field_col = size(trial_pp(1).(fields{ifield}));
        end
        
        field_data_all = [];
        %field_data_var = [];
        
        for isubject = 1:length(trial_pp)
            iang = find([trial_pp(isubject).angles] == angles(iangle));
            chanControl = Aparams_pp(isubject).chanControl;

            if length(field_str) == 2
                if isrow(trial_pp(isubject).(field_str{1}).(field_str{2}))
                    field_data_all = [field_data_all; trial_pp(isubject).(field_str{1}).(field_str{2})(iang)];
                else
                    if strcmp(field_str{1},'EMG')
                        field_data_all = [field_data_all; trial_pp(isubject).(field_str{1}).(field_str{2})(iang,chanControl)];
                    else
                        field_data_all = [field_data_all; trial_pp(isubject).(field_str{1}).(field_str{2})(iang,:)];
                    end
                end
            else
                if isrow(trial_pp(isubject).(field_str{1}))
                    field_data_all = [field_data_all; trial_pp(isubject).(field_str{1})(iang)];
                else
                    field_data_all = [field_data_all; trial_pp(isubject).(field_str{1})(iang,:)];
                end
            end
            %field_data_var(isubject,:) = var(angle_data(isubject).(field_str{1}).(field_str{2})(iang,:),1);
        end
        
        if length(field_str) == 2
            trial_pp_avg(iangle).(field_str{1}).([field_str{2},'_mean']) = mean(field_data_all,1);
            trial_pp_avg(iangle).(field_str{1}).([field_str{2},'_std']) = std(field_data_all,1);
        else
            trial_pp_avg(iangle).([field_str{1},'_mean']) = mean(field_data_all,1);
            trial_pp_avg(iangle).([field_str{1},'_std']) = std(field_data_all,1);
        end
        %trial_pp_avg(iangle).(field_str{1}).([field_str{2},'_pstd']) = sqrt(sum((nsamp-1)*field_data_var,1)./((nsamp-1)*size(field_data_var,1)));
    end
end
end