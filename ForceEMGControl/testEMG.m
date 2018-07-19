load([fileparams.filepath,fileparams.filenameEMG])
fchEMG = 10;
sampleRateEMG = 1024;
wn = [(2/sampleRateEMG)*fchEMG];
[b,a] = butter(2,wn,'high');
filtEMG = filter(b,a,(EMGDataOut),[],2);
 
y = movingavg(abs(filtEMG),200);

figure
plot((filtEMG(1,:)))
hold on
plot(EMGDataOut(1,:),'r')
plot(y(1,:),'g')
