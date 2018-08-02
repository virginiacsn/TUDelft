load([fileparams.filepath,fileparams.filenameEMG])
% load('emg_proc.mat')
EMGdata = (EMGDataOut);
fchEMG = 5;
sampleRateEMG = 1024;
wn = [(2/sampleRateEMG)*fchEMG];
[b,a] = butter(2,wn,'high');
filtEMG = filtfilt(b,a,(EMGdata)')';
 
y = movingavg(abs(filtEMG),200);

figure
plot((filtEMG(1,:)))
hold on
plot(EMGdata(1,:),'r')
plot(y(1,:),'g')
legend('FiltEMG','EMGdata','movavgFiltEMG')

figure
plot(emg_save(1,:))

%%

x = [ sin(0.1:0.002:0.9) 0.9:-0.01:0 zeros(1,200)];
filtx = filter(b,a,x,[],2);

figure
plot(x)
hold on
plot(filtx)


