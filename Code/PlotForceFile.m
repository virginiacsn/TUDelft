date = '20180710';
task = 'CO';
code = '003';
filename = [date,'_',task,'_',code,'.txt'];
filepath = ['D:\Student_experiments\Virginia\FTSensor\Data\',filename];

fileData = importdata(filepath,'\t');
time = fileData.data(:,1);
Fx = fileData.data(:,2);
Fy = fileData.data(:,3);

figure(2)
plot(time,Fx)
hold on
plot(time,Fy,'r')
legend('Fx','Fy')
xlabel('Time [s]'); ylabel('Force [N]');

