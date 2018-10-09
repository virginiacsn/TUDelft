function[trial_data] = trialCO(varargin)

if nargin > 1
    downsamp = 1;
    fs = 1024;
    targetAngles = [];
    block = 1;
    EMGOffset = [];
    struct2vars(who,varargin{2});
end

forceDataOut = varargin{1}{1};
if length(varargin{1}) > 1
    EMGDataOutCell = varargin{1}{2};
    
    if iscell(EMGDataOutCell)
        switch computer
            case 'MACI64'
                EMGDataOut = cell2mat(EMGDataOutCell(2:end))';
            case 'PCWIN64'
                EMGDataOut = cell2mat(EMGDataOutCell(2:end))';
            case 'PCWIN'
                EMGDataOut1 = cell2mat(EMGDataOutCell(2:round(length(EMGDataOutCell)/2)));
                EMGDataOut2 = cell2mat(EMGDataOutCell(round(length(EMGDataOutCell)/2)+1:end));
                
                EMGDataOut = [EMGDataOut1; EMGDataOut2]';
        end
    end
    
    EMG_data = EMGDataOut(1:end-1,:);
    
    EMG_trigger = EMGDataOut(end,:);
    EMG_trigger(EMG_trigger ~= 6) = 0;
    EMG_trigger(EMG_trigger > 0) = 1;
    
    EMGdiff = diff(EMG_trigger);
    EMGpulseind = find(EMGdiff>0);
    EMGindext = EMGpulseind(diff(EMGpulseind)>1024)+1;
    EMG_trigger(EMGindext) = [];
    EMG_data(:,EMGindext) = [];
    
    EMG_data = EMG_data(:,find(EMG_trigger,1):end);
    EMG_trigger = EMG_trigger(:,find(EMG_trigger,1):end);
    
else
    EMG_data = [];
end

% forceData = forceDataOut;
% EMGData = EMGDataOut;

header = forceDataOut(1,:);
%header_data = header(2:end);
force_data = forceDataOut(2:end,:);

idh_tn = ismember(header,'Trialnum');
idh_ang = ismember(header,'TargetAng');
idh_sta = ismember(header,'State');
idh_ts = ismember(header,'TimeStamp');
idh_fx = ismember(header,'Fx');
idh_fy = ismember(header,'Fy');
idh_fz = ismember(header,'Fz');
idh_trig = ismember(header,'Trigger');

% tsinit = cell2mat(force_data(1,idh_ts));
% tinit = tsinit(1);

trial_data = struct();
trial_nums = cell2mat(force_data(:,idh_tn));
% trial_numsd = zeros(length(trial_nums),1);
% trial_numsd([1;find(diff(trial_nums))]) = 1;
% trial_forces = cell2mat(force_data(:,idh_fx|idh_fy|idh_fz));
% trial_trig = cell2mat(force_data(:,idh_trig));
num_samp = size(force_data{1,idh_ts},1);
%trial_states = {'start','movement','hold','success','fail'};
trial_outcomes = {'S','F','I'};
%trial_data_fields = {'outcome','tstart','tmove','thold','tend','force'}; % add EMG
trials = 0:max(trial_nums);
cum_samp_trial = 0;

for itrial = 1:length(trials)
    trial_data(itrial).block = block;
    
    data_trial = force_data(trial_nums==trials(itrial),2:end);
    data_outcome = data_trial(:,idh_sta(2:end));
    data_ts = data_trial(:,idh_ts(2:end));
    data_force = cat(2,data_trial(:,idh_fx(2:end)),data_trial(:,idh_fy(2:end)),data_trial(:,idh_fz(2:end)));
    data_trigger = data_trial(:,idh_trig(2:end));

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

    istart = find(ismember(data_outcome,'start'),1);
    imove = find(ismember(data_outcome,'movement'),1);
    ihold = find(ismember(data_outcome,'hold'),1);
    
    if ~isempty(imove)
        trial_data(itrial).angle = targetAngles(data_trial{2,idh_ang(2:end)});
    else
        trial_data(itrial).angle = targetAngles(data_trial{1,idh_ang(2:end)});
    end
    
    trial_data(itrial).dt = 1/fs;
    
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
        trial_data(itrial).tend = tsend(end);
        
        %trial_data(itrial).ts = cell2mat(data_ts(istart:iend));
        trial_data(itrial).force.raw = cell2mat(data_force(istart:iend,:));
        trial_data(itrial).trigger.force = cell2mat(data_trigger(istart:iend,:));
    elseif isempty(iend) && ~isempty(istart)
        trial_data(itrial).iend = length(data_ts);
        tsend = data_ts{end};
        trial_data(itrial).tend = tsend(end);
        
        %trial_data(itrial).ts = cell2mat(data_ts(istart:end));
        trial_data(itrial).force.raw = cell2mat(data_force(istart:end,:));
        trial_data(itrial).trigger.force = cell2mat(data_trigger(istart:end,:));
    elseif ~isempty(iend) && isempty(istart)
        tsend = data_ts{iend};
        trial_data(itrial).iend = iend;
        trial_data(itrial).tend = tsend(end);
        
        %trial_data(itrial).ts = cell2mat(data_ts(1:iend));
        trial_data(itrial).force.raw = cell2mat(data_force(1:iend,:));
        trial_data(itrial).trigger.force = cell2mat(data_trigger(1:iend,:));
    end
    
    trial_data(itrial).imove =  trial_data(itrial).imove*num_samp;
    trial_data(itrial).ihold = trial_data(itrial).ihold*num_samp;
    trial_data(itrial).iend =  trial_data(itrial).iend*num_samp;
    
    trial_data(itrial).trigger.force(find(trial_data(itrial).trigger.force < 3)) = 0;
    trial_data(itrial).trigger.force(find(trial_data(itrial).trigger.force >= 3)) = 1;
    
    
    if ~isempty(downsamp)
        %trial_data(itrial).ts = decimate(trial_data(itrial).ts,downsample,'fir');
        trial_data(itrial).imove =  trial_data(itrial).imove/downsamp;
        trial_data(itrial).ihold =  trial_data(itrial).ihold/downsamp;
        trial_data(itrial).iend =  trial_data(itrial).iend/downsamp;
        
        downsamp_fx = decimate(trial_data(itrial).force.raw(:,1),downsamp,'fir');
        downsamp_fy = decimate(trial_data(itrial).force.raw(:,2),downsamp,'fir');
        downsamp_fz = decimate(trial_data(itrial).force.raw(:,3),downsamp,'fir');
        trial_data(itrial).force.raw = [downsamp_fx downsamp_fy downsamp_fz];
        
        trial_data(itrial).trigger.force = downsample(trial_data(itrial).trigger.force,downsamp);
    end
    
    samp_trial = length(trial_data(itrial).trigger.force);
    samp_trial_tot = ceil(length(cell2mat(data_ts))/downsamp);
    
    trial_data(itrial).EMG.offset = EMGOffset';
    
    if size(EMG_data,2) < cum_samp_trial+samp_trial
        trial_data(itrial).outcome = trial_outcomes{3};
    else
        trial_data(itrial).EMG.raw = EMG_data(:,1+cum_samp_trial:cum_samp_trial+samp_trial)';
        trial_data(itrial).trigger.EMG = EMG_trigger(1+cum_samp_trial:cum_samp_trial+samp_trial)';
        
        cum_samp_trial = cum_samp_trial+samp_trial_tot;
        
        trial_data(itrial).ts = (0:length(trial_data(itrial).EMG.raw)-1)*trial_data(itrial).dt;
        trial_data(itrial).fv = (0:length(trial_data(itrial).EMG.raw)-1)/trial_data(itrial).ts(end);
    end
end

nfields = length(fieldnames(trial_data));
trial_data = orderfields(trial_data,[1:nfields-5 nfields-1 nfields nfields-4 nfields-2 nfields-3]);
end