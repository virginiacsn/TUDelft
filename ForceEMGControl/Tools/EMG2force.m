function[H] = EMG2force(forceData,EMGData)

H = regress(forceData',[ones(size(EMGData,1),1) EMGData]'); 
end