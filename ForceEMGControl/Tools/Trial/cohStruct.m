function[trial_data_coh] = cohStruct(varargin)

trial_data = varargin{1};
channelName = varargin{2};
EMG_fields = varargin{3};


window = [];
tseg = 1;
nseg = 10;
my_nseg = 20;
fs = 1024;
alp = 0.05;

alp_band = [8 12];
beta_band = [15 30];
gam_band = [30 80];

if length(varargin)>3
    struct2vars(who,varargin{4});
end

for i = 1:length(trial_data)
    trial_data_coh(i).angle = trial_data(i).angle;
    
    if isempty(window)
        win = hanning(128);
    else
        win = window(round(length(trial_data(i).ts)/nseg));
        my_win = window(round(length(trial_data(i).ts)/my_nseg));
        %win = window(round(tseg*trial_data(i).ts(2)));
    end
    overlap = round(length(win)/2);
    my_overlap = round(length(my_win)/2);
    
    EMG_struct = trial_data(i).EMG;
    if isfield(trial_data,'iapp');
        iapp = trial_data(i).iapp;
    end
    
    for j = 1:length(EMG_fields)
        for k = 1:length(channelName)-1
            for l = k+1:length(channelName)-1
                trial_data_coh(i).(EMG_fields{j}).muscles{k+l-2} = {channelName{k},channelName{l}};
                [trial_data_coh(i).(EMG_fields{j}).coh(:,k+l-2),trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)] = mscohere(EMG_struct.(EMG_fields{j})(:,k),EMG_struct.(EMG_fields{j})(:,l),win,overlap,[],fs);

                L = round(length(EMG_struct.(EMG_fields{j})(:,k))/length(win));%*trial_data(i).ntrials;
                trial_data_coh(i).(EMG_fields{j}).CL(k+l-2) = 1-alp^(1/(L-1)); 
                
                if isfield(trial_data,'iapp')
                    [trial_data_coh(i).(EMG_fields{j}).my_coh(:,k+l-2),trial_data_coh(i).(EMG_fields{j}).my_fcoh(:,k+l-2),my_nsegtot] = coherence(EMG_struct.(EMG_fields{j})(:,k),EMG_struct.(EMG_fields{j})(:,l),fs,window,my_overlap,my_nseg,iapp);
                else
                    [trial_data_coh(i).(EMG_fields{j}).my_coh(:,k+l-2),trial_data_coh(i).(EMG_fields{j}).my_fcoh(:,k+l-2),my_nsegtot] = coherence(EMG_struct.(EMG_fields{j})(:,k),EMG_struct.(EMG_fields{j})(:,l),fs,my_win,my_overlap,my_nseg);
                end
                trial_data_coh(i).(EMG_fields{j}).my_CL(k+l-2) = 1-alp^(1/(my_nsegtot-1));
                
                alp_freq = find(trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)>=alp_band(1) & trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)<=alp_band(2));
                beta_freq = find(trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)>=beta_band(1) & trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)<=beta_band(2));
                gam_freq = find(trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)>=gam_band(1) & trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)<=gam_band(2));

                coh_temp = zeros(size(trial_data_coh(i).(EMG_fields{j}).coh(:,k+l-2)));
                coh_temp(trial_data_coh(i).(EMG_fields{j}).coh(:,k+l-2) >= trial_data_coh(i).(EMG_fields{j}).CL(k+l-2)) = 1;

                trial_data_coh(i).(EMG_fields{j}).sig_coh(:,k+l-2) = [mean(coh_temp(alp_freq)) mean(coh_temp(beta_freq)) mean(coh_temp(gam_freq)) std(coh_temp(alp_freq)) std(coh_temp(beta_freq)) std(coh_temp(gam_freq))];

                %[~,fcoh250] = min(abs(trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)-250));
                %[~,fcoh500] = min(abs(trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)-500));
                %trial_data_coh(i).(EMG_fields{j}).z(:,k+l-2) = (atanh(sqrt(trial_data_coh(i).(EMG_fields{j}).coh(:,k+l-2))/sqrt(1/(2*L)))-mean(trial_data_coh(i).(EMG_fields{j}).coh(fcoh250:fcoh500,k+l-2)));
                trial_data_coh(i).(EMG_fields{j}).z(:,k+l-2) = sqrt(2*L)*atanh(sqrt(trial_data_coh(i).(EMG_fields{j}).coh(:,k+l-2)));

                [trial_data_coh(i).(EMG_fields{j}).corr(:,k+l-2),trial_data_coh(i).(EMG_fields{j}).lags(:,k+l-2)] = xcorr(EMG_struct.(EMG_fields{j})(:,k),EMG_struct.(EMG_fields{j})(:,l),'unbiased');

            end
        end
    end
end
end