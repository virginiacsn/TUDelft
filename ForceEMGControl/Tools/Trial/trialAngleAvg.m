function[trial_data_avg] = trialAngleAvg(trial_data, epoch, fields)

angles = sort(unique(extractfield(trial_data,'angle')));

if ~isempty(epoch)
    if length(epoch) == 2
            nsamp = min([trial_data.(epoch{2})]-[trial_data.(epoch{1})])+1;
    else
            nsamp = min([trial_data.(epoch{3})]+round(epoch{4}/trial_data(1).dt)-([trial_data.(epoch{1})]+round(epoch{2}/trial_data(1).dt)))+1;
    end
end

for iangle = 1:length(angles)
    angle_data = trial_data(find(extractfield(trial_data,'angle') == angles(iangle)));
    trial_data_avg(iangle).angle = angles(iangle);
    trial_data_avg(iangle).ntrials = length(angle_data);
    
    for ifield = 1:length(fields)
        if any(strfind(fields{ifield},'.'))
            field_str = strsplit(fields{ifield},'.');
            field_col = size(angle_data(1).(field_str{1}).(field_str{2}),2);
        else
            field_col = size(angle_data(1).(fields{ifield}),2);
        end
        field_data_mean = [];
        field_data = [];
        %field_data_fft = [];
        
        for itrial = 1:length(angle_data)
            if ~isempty(epoch)
                if length(epoch) == 2
                    idx1 = angle_data(itrial).(epoch{1});
                    idx2 = angle_data(itrial).(epoch{2});
                else
                    idx1 = angle_data(itrial).(epoch{1})+round(epoch{2}/angle_data(itrial).dt);
                    idx2 = angle_data(itrial).(epoch{3})+round(epoch{4}/angle_data(itrial).dt);
                end
                
                nsampextra = (idx2-idx1)-nsamp+1;
                sampv = idx1+round(nsampextra/2):idx2-round(nsampextra/2);
            else
                nsamp = size(angle_data(itrial).(field_str{1}).(field_str{2}),1);
                sampv = 1:nsamp;
            end
            
            field_data_mean(itrial,:) = mean(angle_data(itrial).(field_str{1}).(field_str{2})(sampv,:),1);
            field_data = [field_data reshape(angle_data(itrial).(field_str{1}).(field_str{2})(sampv,:),field_col*length(sampv),1)];
        end
        
        if isfield(angle_data(itrial),'dt')
            trial_data_avg(iangle).ts = (0:nsamp-1)*angle_data(1).dt;
            trial_data_avg(iangle).fv = (0:nsamp-1)/trial_data_avg(iangle).ts(end);
        elseif isfield(angle_data(itrial).(field_str{1}),'fcoh')
            trial_data_avg(iangle).(field_str{1}).muscles = angle_data(iangle).(field_str{1}).muscles;
            trial_data_avg(iangle).(field_str{1}).fcoh = angle_data(iangle).(field_str{1}).fcoh;
            trial_data_avg(iangle).(field_str{1}).CL = angle_data(iangle).(field_str{1}).CL;
        end
        
        trial_data_avg(iangle).(field_str{1}).([field_str{2},'_mean']) = mean(field_data_mean,1);
        trial_data_avg(iangle).(field_str{1}).(field_str{2}) = reshape(mean(field_data,2),length(sampv),field_col);
        %trial_data_avg(iangle).(field_str{1}).(field_str{2}) = ifft(reshape(mean(field_data_fft,2),length(sampv),field_col),[],2);%reshape(mean(field_data,2),length(sampv),field_col);
        %trial_data_avg(iangle).(field_str{1}).([field_str{2},'_fft']) = fft(trial_data_avg(iangle).(field_str{1}).(field_str{2}));
        %trial_data_avg(iangle).(field_str{1}).([field_str{2},'_fft']) = reshape(mean(field_data_fft,2),length(sampv),field_col);
    end
%     for ifield = 1:length(fields)
%         field_str = strsplit(fields{ifield},'.');
%         if strcmp(field_str{2},'rect')
%             trial_data_avg(iangle).(field_str{1}).([field_str{2},'_mean']) = mean(field_data_mean,1);
%             trial_data_avg(iangle).(field_str{1}).(field_str{2}) = reshape(mean(field_data,2),length(sampv),field_col);
%             %trial_data_avg(iangle).(field_str{1}).(field_str{2}) = ifft(reshape(mean(field_data_fft,2),length(sampv),field_col),[],2);%reshape(mean(field_data,2),length(sampv),field_col);
%             trial_data_avg(iangle).(field_str{1}).([field_str{2},'_fft']) = fft(trial_data_avg(iangle).(field_str{1}).(field_str{2}));
%         elseif strcmp(field_str{2},'avg')
%             trial_data_avg(iangle).(field_str{1}).([field_str{2},'_mean']) = mean(field_data_mean,1);
%             trial_data_avg(iangle).(field_str{1}).(field_str{2}) = reshape(mean(field_data,2),length(sampv),field_col);
%             %trial_data_avg(iangle).(field_str{1}).(field_str{2}) = ifft(reshape(mean(field_data_fft,2),length(sampv),field_col),[],2);%reshape(mean(field_data,2),length(sampv),field_col);
%             trial_data_avg(iangle).(field_str{1}).([field_str{2},'_fft']) = fft(trial_data_avg(iangle).(field_str{1}).(field_str{2}));   
%         end
%     end
end
