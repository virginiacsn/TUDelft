function[trial_data] = procForce(trial_data,Aparams)
% Initialize variables
fclF = 5;
fs = 2048;
struct2vars(who,Aparams)

% Filtering parameters
wn = (2/fs)*fclF;
[b,a] = butter(2,wn,'low');

for i = 1:length(trial_data)
    trial_data(i).force.rawmag = sqrt(trial_data(i).force.raw(:,1).^2+trial_data(i).force.raw(:,2).^2);
    trial_data(i).force.filt = filtfilt(b,a,trial_data(i).force.raw);
    trial_data(i).force.filtmag = sqrt(trial_data(i).force.filt(:,1).^2+trial_data(i).force.filt(:,2).^2);
end
end