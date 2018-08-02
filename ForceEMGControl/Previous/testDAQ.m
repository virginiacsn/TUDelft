function[storedata,storetime]= testDAQ()
device = daq.getDevices;
s = daq.createSession('ni');

addAnalogInputChannel(s,device.ID,0:6,'Voltage');
addAnalogOutputChannel(s,device.ID,'ao0','Voltage');

% outputSingleScan(s,1);
% wait(s,0.1);
% outputSingleScan(s,0);
s.Rate = 2000;

availsam = 200;
outputData = [5*ones(1,availsam),zeros(1,s.Rate-availsam)]';
queueOutputData(s,outputData);

data = s.startForeground();
outputData = [5*ones(1,availsam),zeros(1,s.Rate-availsam)]';
queueOutputData(s,outputData);
lh1 = addlistener(s,'DataRequired',@(src,event) src.queueOutputData(outputData));
s.NotifyWhenDataAvailableExceeds = 200;

storedata = [];
storetime = [];

lh = addlistener(s,'DataAvailable',@saveData);
s.IsContinuous = true;

s.startBackground();

input('\Press enter to stop acquisition.')


s.stop()
delete(lh);
delete(lh1);

function saveData(src,event)
        storedata = [storedata event.Data(:,7)'];
        storetime = [storetime event.TimeStamps'];
end
end
