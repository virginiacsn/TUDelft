function[trial_data_coh] = cohStruct(varargin)

trial_data = varargin{1};
channelName = varargin{2};
EMG_fields = varargin{3};


window = [];
tseg = 1;
nseg = 10;
fs = 1024;
alp = 0.05;
if length(varargin)>3
    struct2vars(who,varargin{4});
end

for i = 1:length(trial_data)
    trial_data_coh(i).angle = trial_data(i).angle;
    
    if isempty(window)
        win = hanning(128);
    else
        win = window(round(length(trial_data(i).ts)/nseg));
        %win = window(round(tseg*trial_data(i).ts(2)));
    end
    overlap = round(length(win)/2);
    
    EMG_struct = trial_data(i).EMG;
    
    for j = 1:length(EMG_fields)
        for k = 1:length(channelName)-1
            for l = k+1:length(channelName)-1
                trial_data_coh(i).(EMG_fields{j}).muscles{k+l-2} = {channelName{k},channelName{l}};
                [trial_data_coh(i).(EMG_fields{j}).coh(:,k+l-2),trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)] = mscohere(EMG_struct.(EMG_fields{j})(:,k),EMG_struct.(EMG_fields{j})(:,l),win,overlap,[],fs);
                %trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2) = trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2);%/pi*fs;
                [trial_data_coh(i).(EMG_fields{j}).corr(:,k+l-2),trial_data_coh(i).(EMG_fields{j}).lags(:,k+l-2)] = xcorr(EMG_struct.(EMG_fields{j})(:,k),EMG_struct.(EMG_fields{j})(:,l),'unbiased');

                L = round(length(EMG_struct.(EMG_fields{j})(:,k))/length(win));
                trial_data_coh(i).(EMG_fields{j}).CL(k+l-2) = 1-alp^(1/(L-1)); 
                
                %[~,fcoh250] = min(abs(trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)-250));
                %[~,fcoh500] = min(abs(trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)-500));
                %trial_data_coh(i).(EMG_fields{j}).z(:,k+l-2) = (atanh(sqrt(trial_data_coh(i).(EMG_fields{j}).coh(:,k+l-2))/sqrt(1/(2*L)))-mean(trial_data_coh(i).(EMG_fields{j}).coh(fcoh250:fcoh500,k+l-2)));
                trial_data_coh(i).(EMG_fields{j}).z(:,k+l-2) = sqrt(2*L)*atanh(sqrt(trial_data_coh(i).(EMG_fields{j}).coh(:,k+l-2)));
            end
        end
    end
end
end