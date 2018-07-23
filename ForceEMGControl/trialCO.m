function[trial_data] = trialCO(varargin)

forceData = varargin{1}{1};
if length(varargin{1})>1
    EMGDataOut = varargin{1}{2};
    EMGData = EMGDataOut(1:2,:);
    EMGTrig = EMGDataOut(3,:);
    EMGTrig(EMGTrig ~= 6) = 0;
    EMGTrig(find(EMGTrig)) = 1;
    EMGData = EMGData(:,find(EMGTrig,1):end);
    EMGTrig = EMGTrig(:,find(EMGTrig,1):end);
else
    EMGData = [];
end
if nargin > 1
    downsample = [];
    target_angles = [0:2*pi/8:2*pi];
    struct2vars(who,varargin{2});
end
% forceData = forceDataOut;
% EMGData = EMGDataOut;

header = forceData(1,:);
header_data = header(2:end);
forceData = forceData(2:end,:);

idh_tn = ismember(header,'Trialnum');
idh_ang = ismember(header,'TargetAng');
idh_sta = ismember(header,'State');
idh_ts = ismember(header,'TimeStamp');
idh_fx = ismember(header,'Fx');
idh_fy = ismember(header,'Fy');
idh_fz = ismember(header,'Fz');
idh_trig = ismember(header,'Trigger');

tsinit = cell2mat(forceData(1,idh_ts));
tinit = tsinit(1);

trial_data = struct();
trial_nums = cell2mat(forceData(:,idh_tn));
trial_numsd = zeros(length(trial_nums),1);
trial_numsd([1;find(diff(trial_nums))]) = 1;
trial_forces = cell2mat(forceData(:,idh_fx|idh_fy|idh_fz));
trial_trig = cell2mat(forceData(:,idh_trig));
sampLength = size(forceData{1,idh_ts},1);
%trial_states = {'start','movement','hold','success','fail'};
trial_outcomes = {'S','F','I'};
%trial_data_fields = {'outcome','tstart','tmove','thold','tend','force'}; % add EMG
trials = 0:max(trial_nums);
cum_samp_trial = 0;

for itrial = 1:length(trials)
    data_trial = forceData(trial_nums==trials(itrial),2:end);
    data_outcome = data_trial(:,idh_sta(2:end));
    data_ts = (data_trial(:,idh_ts(2:end)));
    data_force = (cat(2,data_trial(:,idh_fx(2:end)),data_trial(:,idh_fy(2:end)),data_trial(:,idh_fz(2:end))));
    data_trigger = (data_trial(:,idh_trig(2:end)));
    
    if any(strcmp(data_outcome,'success'))
        trial_data(itrial).outcome = trial_outcomes{1};
        iend = find(ismember(data_outcome,'success'),1);
    elseif any(strcmp(data_outcome,'fail'))
        trial_data(itrial).outcome = trial_outcomes{2};
        iend = find(ismember(data_outcome,'fail'),1);
    else
        trial_data(itrial).outcome = trial_outcomes{3};
        iend = [];
    end
    
    trial_data(itrial).angle = target_angles(data_trial{1,idh_ang(2:end)});
    
    trial_data(itrial).dt = mean(diff(data_ts{1}));
    
    istart = find(ismember(data_outcome,'start'),1);
    imove = find(ismember(data_outcome,'movement'),1);
    ihold = find(ismember(data_outcome,'hold'),1);
    
    if ~isempty(istart) && ~isempty(imove)
        tsstart = data_ts{istart};        
        trial_data(itrial).tstart = tsstart(1);
        
        trial_data(itrial).istart = istart;
        trial_data(itrial).imove  = imove;
    elseif isempty(istart) && trials(itrial) == 0 
        trial_data(itrial).tstart = 0; % Hard coded
        
        trial_data(itrial).istart = 1;
        trial_data(itrial).imove  = 2;
    else
        trial_data(itrial).tstart = [];
    end
    
    if isempty(imove)
        trial_data(itrial).imove = [];
    end
    
    if ~isempty(ihold)
        trial_data(itrial).ihold  = ihold;
    else
        trial_data(itrial).ihold = [];
    end
    
    if ~isempty(iend) && ~isempty(istart)
        tsend = data_ts{iend};
        trial_data(itrial).iend = iend;
        trial_data(itrial).tend = tsend(1);
        
        %trial_data(itrial).ts = cell2mat(data_ts(istart:iend));
        trial_data(itrial).force = cell2mat(data_force(istart:iend,:));
        trial_data(itrial).trigger = cell2mat(data_trigger(istart:iend,:));
    elseif isempty(iend) && ~isempty(istart)
        trial_data(itrial).iend = length(data_ts);
        trial_data(itrial).tend = data_ts(end);
        
        %trial_data(itrial).ts = cell2mat(data_ts(istart:end));
        trial_data(itrial).force = cell2mat(data_force(istart:end,:));
        trial_data(itrial).trigger = cell2mat(data_trigger(istart:end,:));
    elseif ~isempty(iend) && isempty(istart)
        tsend = data_ts{iend};
        trial_data(itrial).iend = iend;
        trial_data(itrial).tend = tsend(1);
        
        %trial_data(itrial).ts = cell2mat(data_ts(1:iend));
        trial_data(itrial).force = cell2mat(data_force(1:iend,:));
        trial_data(itrial).trigger = cell2mat(data_trigger(1:iend,:));
    end
    
    trial_data(itrial).imove =  trial_data(itrial).imove*sampLength;
    trial_data(itrial).ihold = trial_data(itrial).ihold*sampLength;
    trial_data(itrial).iend =  trial_data(itrial).iend*sampLength;
    
    trial_data(itrial).trigger(find(trial_data(itrial).trigger<3)) = 0;
    trial_data(itrial).trigger(find(trial_data(itrial).trigger>3)) = 1;
    
    
    if ~isempty(downsample)
        %trial_data(itrial).ts = decimate(trial_data(itrial).ts,downsample,'fir');
        trial_data(itrial).imove =  trial_data(itrial).imove/downsample;
        trial_data(itrial).ihold =  trial_data(itrial).ihold/downsample;
        trial_data(itrial).iend =  trial_data(itrial).iend/downsample;
        
        downsamp_fx = decimate(trial_data(itrial).force(:,1),downsample,'fir');
        downsamp_fy = decimate(trial_data(itrial).force(:,2),downsample,'fir');
        downsamp_fz = decimate(trial_data(itrial).force(:,3),downsample,'fir');
        trial_data(itrial).force = [downsamp_fx downsamp_fy downsamp_fz];
        
        trial_data(itrial).trigger = decimate(trial_data(itrial).trigger,downsample,'fir');
    end
    
    samp_trial = length(trial_data(itrial).trigger);
    samp_trial_tot = round(length(cell2mat(data_ts))/2);
    trial_data(itrial).EMG = EMGData(:,1+cum_samp_trial:cum_samp_trial+samp_trial)';
    trial_data(itrial).EMGtrigger = EMGTrig(1+cum_samp_trial:cum_samp_trial+samp_trial)';

    cum_samp_trial = cum_samp_trial+samp_trial_tot;
end
end
