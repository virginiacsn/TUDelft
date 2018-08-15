function[trial_data_trim] = trialTrim(trial_data, epoch)

%angles = sort(unique(extractfield(trial_data,'angle')));
nsamp = min(extractfield(trial_data,epoch{2})-extractfield(trial_data,epoch{1}))+1;

for itrial = 1:length(trial_data)
    trial_data_trim(itrial).block = trial_data(itrial).block;
    trial_data_trim(itrial).outcome = trial_data(itrial).outcome;
    trial_data_trim(itrial).angle = trial_data(itrial).angle;
    trial_data_trim(itrial).dt = trial_data(itrial).dt;
    trial_data_trim(itrial).(epoch{1}) = 1;
    trial_data_trim(itrial).(epoch{2}) = nsamp;
    
    idx1 = trial_data(itrial).(epoch{1});
    idx2 = trial_data(itrial).(epoch{2});
    nsampextra = (idx2-idx1)-nsamp+1;
    sampv = idx1+round(nsampextra/2):idx2-round(nsampextra/2);

    trial_data_trim(itrial).ts = (0:nsamp-1)*trial_data_trim(itrial).dt;
    trial_data_trim(itrial).fv = (0:nsamp-1)/trial_data_trim(itrial).ts(end);
    
    fieldsForce = fieldnames(trial_data(itrial).force);
    for i = 1:length(fieldsForce)
        if size(trial_data(itrial).force.(fieldsForce{i}),2)>1
            trial_data_trim(itrial).force.(fieldsForce{i}) = trial_data(itrial).force.(fieldsForce{i})(sampv,:);
        end
    end
    
    fieldsEMG = fieldnames(trial_data(itrial).EMG);
    for i = 1:length(fieldsEMG)
        if size(trial_data(itrial).EMG.(fieldsEMG{i}),2)>1
            trial_data_trim(itrial).EMG.(fieldsEMG{i}) = trial_data(itrial).EMG.(fieldsEMG{i})(sampv,:);
        end
    end
    
    fieldsTrig = fieldnames(trial_data(itrial).trigger);
    for i = 1:length(fieldsTrig)
        trial_data_trim(itrial).trigger.(fieldsTrig{i}) = trial_data(itrial).trigger.(fieldsTrig{i})(sampv,:);
    end
    
end