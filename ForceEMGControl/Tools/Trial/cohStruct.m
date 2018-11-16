function[trial_data_coh] = cohStruct(varargin)

trial_data = varargin{1};
channelName = varargin{2};
EMG_fields = varargin{3};

window = [];
tseg = 1;
fs = 1024;
alp = 0.05;
CLoverlap = 1;

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
        %win = window(round(length(trial_data(i).ts)/nseg));
        %my_win = window(round(length(trial_data(i).ts)/my_nseg));
        win = window(round(tseg/trial_data(i).ts(2)));        
        my_win = window(round(tseg/trial_data(i).ts(2)));
    end
    overlap = round(length(win)/2);
    my_overlap = round(length(my_win)/2);
    
    EMG_struct = trial_data(i).EMG;
    if isfield(trial_data,'iapp')
        iapp = trial_data(i).iapp;
    end
    
    for j = 1:length(EMG_fields)
        h = 0;
        for k = 1:length(channelName)-1
            for l = k+1:length(channelName)-1
                h = h+1;
                trial_data_coh(i).(EMG_fields{j}).muscles{h} = {channelName{k},channelName{l}};
                
%                 [trial_data_coh(i).(EMG_fields{j}).coh(:,h),trial_data_coh(i).(EMG_fields{j}).fcoh(:,h)] = mscohere(EMG_struct.(EMG_fields{j})(:,k),EMG_struct.(EMG_fields{j})(:,l),win,overlap,[],fs);
%                 L = round(length(EMG_struct.(EMG_fields{j})(:,k))/(length(win)-overlap));%*trial_data(i).ntrials;
%                 trial_data_coh(i).(EMG_fields{j}).CL(h) = 1-alp^(1/(L-1)); 
                
                if isfield(trial_data,'iapp')
                    [trial_data_coh(i).(EMG_fields{j}).my_coh(:,h),trial_data_coh(i).(EMG_fields{j}).my_fcoh(:,h),my_nsegtot] = coherence(EMG_struct.(EMG_fields{j})(:,k),EMG_struct.(EMG_fields{j})(:,l),fs,my_win,my_overlap,CLoverlap,iapp);
                else
                    [trial_data_coh(i).(EMG_fields{j}).my_coh(:,h),trial_data_coh(i).(EMG_fields{j}).my_fcoh(:,h),my_nsegtot] = coherence(EMG_struct.(EMG_fields{j})(:,k),EMG_struct.(EMG_fields{j})(:,l),fs,my_win,my_overlap,CLoverlap,[]);
                end
                trial_data_coh(i).(EMG_fields{j}).my_CL(h) = 1-alp^(1/(my_nsegtot-1));
                trial_data_coh(i).(EMG_fields{j}).my_nseg(h) = my_nsegtot;
                
                % z-score coherence
                trial_data_coh(i).(EMG_fields{j}).z(:,h) = sqrt(2*my_nsegtot)*atanh(sqrt(trial_data_coh(i).(EMG_fields{j}).my_coh(:,h)));

%                 % Significant coherence in frequency bands - MATLAB
%                 % coherence
%                 alp_freq = find(trial_data_coh(i).(EMG_fields{j}).fcoh(:,h)>=alp_band(1) & trial_data_coh(i).(EMG_fields{j}).fcoh(:,h)<=alp_band(2));
%                 beta_freq = find(trial_data_coh(i).(EMG_fields{j}).fcoh(:,h)>=beta_band(1) & trial_data_coh(i).(EMG_fields{j}).fcoh(:,h)<=beta_band(2));
%                 gam_freq = find(trial_data_coh(i).(EMG_fields{j}).fcoh(:,h)>=gam_band(1) & trial_data_coh(i).(EMG_fields{j}).fcoh(:,h)<=gam_band(2));
% 
%                 coh_temp = zeros(size(trial_data_coh(i).(EMG_fields{j}).coh(:,h)));
%                 coh_temp(trial_data_coh(i).(EMG_fields{j}).coh(:,h) >= trial_data_coh(i).(EMG_fields{j}).CL(h)) = 1;
%                 
%                 alp_sig = coh_temp(alp_freq);
%                 beta_sig = coh_temp(beta_freq);
%                 gam_sig = coh_temp(gam_freq);
%                 
%                 trial_data_coh(i).(EMG_fields{j}).mean_sig_coh(:,h) = [mean(alp_sig) mean(beta_sig) mean(gam_sig)];
%                 trial_data_coh(i).(EMG_fields{j}).SEM_sig_coh(:,h) = [std(alp_sig)/sqrt(length(alp_freq)) std(beta_sig)/sqrt(length(beta_freq)) std(gam_sig)/sqrt(length(gam_freq))];
                
                % Significant coherence in frequency bands - my coherence
                alp_freq = find(trial_data_coh(i).(EMG_fields{j}).my_fcoh(:,h)>=alp_band(1) & trial_data_coh(i).(EMG_fields{j}).my_fcoh(:,h)<=alp_band(2));
                beta_freq = find(trial_data_coh(i).(EMG_fields{j}).my_fcoh(:,h)>=beta_band(1) & trial_data_coh(i).(EMG_fields{j}).my_fcoh(:,h)<=beta_band(2));
                gam_freq = find(trial_data_coh(i).(EMG_fields{j}).my_fcoh(:,h)>=gam_band(1) & trial_data_coh(i).(EMG_fields{j}).my_fcoh(:,h)<=gam_band(2));

                coh_temp = zeros(size(trial_data_coh(i).(EMG_fields{j}).my_coh(:,h)));
                coh_temp_area = zeros(size(trial_data_coh(i).(EMG_fields{j}).my_coh(:,h)));
                z_temp_area = zeros(size(trial_data_coh(i).(EMG_fields{j}).z(:,h)));
                
                idx_sig = trial_data_coh(i).(EMG_fields{j}).my_coh(:,h) >= trial_data_coh(i).(EMG_fields{j}).my_CL(h);
                coh_temp(idx_sig) = 1;
                coh_temp_area(idx_sig) = trial_data_coh(i).(EMG_fields{j}).my_coh(idx_sig,h);
                z_temp_area(idx_sig) = trial_data_coh(i).(EMG_fields{j}).z(idx_sig,h); 

                alp_sig = coh_temp(alp_freq);
                beta_sig = coh_temp(beta_freq);
                gam_sig = coh_temp(gam_freq);
                
                tot_sig = sum(idx_sig);
                
                trial_data_coh(i).(EMG_fields{j}).msig_coh(:,h) = [mean(alp_sig), mean(beta_sig), mean(gam_sig)];
                trial_data_coh(i).(EMG_fields{j}).SEMsig_coh(:,h) = [std(alp_sig)/sqrt(length(alp_freq)), std(beta_sig)/sqrt(length(beta_freq)), std(gam_sig)/sqrt(length(gam_freq))];

                trial_data_coh(i).(EMG_fields{j}).psig_coh(:,h) = 100*[sum(alp_sig)/tot_sig, sum(beta_sig)/tot_sig, sum(gam_sig)/tot_sig];
                trial_data_coh(i).(EMG_fields{j}).pCI_sig_coh(:,h) = 100*1.96*[sqrt(sum(alp_sig)/tot_sig*(1-sum(alp_sig)/tot_sig)/tot_sig), sqrt(sum(beta_sig)/tot_sig*(1-sum(beta_sig)/tot_sig)/tot_sig), sqrt(sum(gam_sig)/tot_sig*(1-sum(gam_sig)/tot_sig)/tot_sig)];
                
                trial_data_coh(i).(EMG_fields{j}).asig_coh(:,h) = [trapz(trial_data_coh(i).(EMG_fields{j}).my_fcoh(alp_freq,h), coh_temp_area(alp_freq)),...
                    trapz(trial_data_coh(i).(EMG_fields{j}).my_fcoh(beta_freq,h), coh_temp_area(beta_freq)),...
                    trapz(trial_data_coh(i).(EMG_fields{j}).my_fcoh(gam_freq,h), coh_temp_area(gam_freq))];
                
                trial_data_coh(i).(EMG_fields{j}).asig_z(:,h) = [trapz(trial_data_coh(i).(EMG_fields{j}).my_fcoh(alp_freq,h), z_temp_area(alp_freq)),...
                    trapz(trial_data_coh(i).(EMG_fields{j}).my_fcoh(beta_freq,h), z_temp_area(beta_freq)),...
                    trapz(trial_data_coh(i).(EMG_fields{j}).my_fcoh(gam_freq,h), z_temp_area(gam_freq))];
                
                trial_data_coh(i).(EMG_fields{j}).nasig_coh(:,h) = trial_data_coh(i).(EMG_fields{j}).asig_coh(:,h)./[length(alp_freq);length(beta_freq);length(gam_freq)];  
                trial_data_coh(i).(EMG_fields{j}).nasig_z(:,h) = trial_data_coh(i).(EMG_fields{j}).asig_z(:,h)./[length(alp_freq);length(beta_freq);length(gam_freq)];
                
                %[~,fcoh250] = min(abs(trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)-250));
                %[~,fcoh500] = min(abs(trial_data_coh(i).(EMG_fields{j}).fcoh(:,k+l-2)-500));
                %trial_data_coh(i).(EMG_fields{j}).z(:,k+l-2) = (atanh(sqrt(trial_data_coh(i).(EMG_fields{j}).coh(:,k+l-2))/sqrt(1/(2*L)))-mean(trial_data_coh(i).(EMG_fields{j}).coh(fcoh250:fcoh500,k+l-2)));
                
                if isfield(trial_data,'iapp')
                    [trial_data_coh(i).(EMG_fields{j}).corr(:,h),trial_data_coh(i).(EMG_fields{j}).lags(:,h)] = correlation(EMG_struct.(EMG_fields{j})(:,k),EMG_struct.(EMG_fields{j})(:,l),iapp);
                else
                    [trial_data_coh(i).(EMG_fields{j}).corr(:,h),trial_data_coh(i).(EMG_fields{j}).lags(:,h)] = xcorr(EMG_struct.(EMG_fields{j})(:,k),EMG_struct.(EMG_fields{j})(:,l),'unbiased');
                end
            end
        end
    end
end
end