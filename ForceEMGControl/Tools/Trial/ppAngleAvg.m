function[trial_pp_avg] = ppAngleAvg(trial_pp, fields, Aparams_pp)

angles = sort(unique(extractfield(trial_pp,'angle')));

for iangle = 1:length(angles)
    
    trial_pp_avg(iangle).angle = angles(iangle);
    
    for ifield = 1:length(fields)
        if any(strfind(fields{ifield},'.'))
            field_str = strsplit(fields{ifield},'.');
            if length(field_str) == 2
                field_col = size(trial_pp(1).(field_str{1}).(field_str{2}),2);
            else
                error('Only 2 levels of substructs allowed.');
            end
        else
            field_str = {fields{ifield}};
            field_col = size(trial_pp(1).(fields{ifield}));
        end
        
        field_data_all = [];
        field_data_all_scale = [];
        scount = 0;
        
        for isubject = 1:length(trial_pp)
            iang = find([trial_pp(isubject).angle] == angles(iangle));
            chanControl = Aparams_pp(isubject).chanControl;
            if isempty(iang) && (ifield == 1)
                fprintf('\nMissing target: Subject %d, angle %d',isubject,rad2deg(angles(iangle)));
            else
                scount = scount+1;
            end
            if isrow(trial_pp(isubject).(field_str{1}).(field_str{2}))
                field_data_all = [field_data_all; trial_pp(isubject).(field_str{1}).(field_str{2})(iang)];
            else
                if strcmp(field_str{1},'EMG')
                    field_data_all = [field_data_all; trial_pp(isubject).(field_str{1}).(field_str{2})(iang,chanControl)];
                    field_data_all_scale = [field_data_all_scale; trial_pp(isubject).(field_str{1}).(field_str{2})(iang,chanControl)./Aparams_pp(isubject).EMGScale(chanControl)];
                else
                    field_data_all = [field_data_all; trial_pp(isubject).(field_str{1}).(field_str{2})(iang,:)];
                end
            end
            
            %field_data_var(isubject,:) = var(angle_data(isubject).(field_str{1}).(field_str{2})(iang,:),1);
        end
        
        trial_pp_avg(iangle).(field_str{1}).([field_str{2},'_mean']) = mean(field_data_all,1);
        trial_pp_avg(iangle).(field_str{1}).([field_str{2},'_std']) = std(field_data_all,1);
        trial_pp_avg(iangle).(field_str{1}).([field_str{2},'_sem']) = std(field_data_all,1)/sqrt(scount);
        if strcmp(field_str{1},'EMG')
            trial_pp_avg(iangle).(field_str{1}).([field_str{2},'_scale_mean']) = mean(field_data_all_scale,1);
            trial_pp_avg(iangle).(field_str{1}).([field_str{2},'_scale_std']) = std(field_data_all_scale,1);
            trial_pp_avg(iangle).(field_str{1}).([field_str{2},'_scale_sem']) = std(field_data_all_scale,1)/sqrt(scount);
        end
        %trial_pp_avg(iangle).(field_str{1}).([field_str{2},'_pstd']) = sqrt(sum((nsamp-1)*field_data_var,1)./((nsamp-1)*size(field_data_var,1)));
    end
end
fprintf('\n');
end