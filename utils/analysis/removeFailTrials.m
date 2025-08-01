function[trial_data_out] = removeFailTrials(trial_data)

i = 0;

for itrial = 1:size(trial_data,2)
    if strcmp(trial_data(itrial).outcome,'S')
        i = i+1;
        trial_data_out(i) = trial_data(itrial);
    end
end
if i == 0
    trial_data_out = [];
    disp('No successfull trials.')
else
    fprintf('Removed %d failed trials.\n',length(trial_data)-i)
end
end