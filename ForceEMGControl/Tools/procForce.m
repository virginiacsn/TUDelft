function[trial_data] = procForce(trial_data,Aparams)
% Initialize variables
fclF = 5;
sampleRateForce = 2048;
struct2vars(who,Aparams)

% Filtering parameters
wn = (2/sampleRateForce)*fclF;
[b,a] = butter(2,wn,'low');

for i = 1:length(trial_data)
    trial_data(i).forcefilt = filtfilt(b,a,trial_data(i).force);
    trial_data(i).forcemag = sqrt(trial_data(i).forcefilt(:,1).^2+trial_data(i).forcefilt(:,2).^2);
end
end