% Tables for CV
field_tab = {'force.filtmag_CV','force.filt_CV','EMG.rect_CV'};
for iang = 1:length(Aparams.angCompUni)
    iangf = find([trial_avg_force.angle] == Aparams.angCompUni(iang));
    iangE = find([trial_avg_EMG.angle] == Aparams.angCompUni(iang));
    fprintf('\nAngle %d:\n',rad2deg(Aparams.angCompUni(iang)));
    for ifield = 1:length(field_tab)
        
        field_str = strsplit(field_tab{ifield},'.');
        if strcmp(field_str{1},'force')
            fprintf('\nField %s:\n',field_tab{ifield});
            
            for j = 1:length(trial_avg_force(iangf).(field_str{1}).(field_str{2}))
                fprintf('FC: %1.1f MC: %1.1f\n',trial_avg_force(iangf).(field_str{1}).(field_str{2})(j),trial_avg_EMG(iangE).(field_str{1}).(field_str{2})(j));
            end
        else
            fprintf('\nField %s:\n',field_tab{ifield});
            for j = 1:length(Aparams.chanControl)
                fprintf('Musc %s:\n',Aparams.channelNameEMG{Aparams.chanControl(j)});
                fprintf('FC: %1.1f MC: %1.1f\n',trial_avg_force(iangf).(field_str{1}).(field_str{2})(Aparams.chanControl(j)),trial_avg_EMG(iangE).(field_str{1}).(field_str{2})(Aparams.chanControl(j)));
            end
        end
    end
end