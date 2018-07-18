load([fileparams.filepath,fileparams.filenameEMG])
fchEMG = 20;
sampleRateEMG = 1024;
wn = [(2/sampleRateEMG)*fchEMG];
[b,a] = butter(2,wn,'high');
filtEMG = filter(b,a,EMGDataOut,[],2);

figure
plot(abs(filtEMG(1,:)))

