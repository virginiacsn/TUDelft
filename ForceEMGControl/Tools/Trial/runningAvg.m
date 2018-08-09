function[y] = runningAvg(x,window,overlap)
sampleRateEMG = 1024;
if size(x,2)<size(x,1)
    xtemp = x'; % Row is channel
else
    xtemp = x;
end

buffer = zeros(size(xtemp,1),window);

j = 0;
ytemp = [];

for n = 0:floor(size(xtemp,2)/(overlap))-1
    j = j+1;
    bufferTemp = buffer;
    bufferTemp(:,1:overlap) = xtemp(:,1+(overlap)*n:overlap*(n+1));
    bufferTemp(:,overlap+1:end) = buffer(:,1:window-overlap);
    buffer = bufferTemp;
    
    wnh = (2/sampleRateEMG)*10;
    wnl = (2/sampleRateEMG)*30;
    [b,a] = butter(2,wnh,'high');
    [d,c] = butter(2,wnl,'low');
    
    filtEMGBuffer = filtfilt(b,a,buffer')';
    filtEMGBuffer = filter(d,c,filtEMGBuffer,[],1);
    rectEMG = abs(filtEMGBuffer);

    ytemp(:,j) = mean(rectEMG,2);
end
        
% nytemp = floor(size(xtemp,2)/window)*window;
% if nytemp ~= size(xtemp,2)
%     ytemp(:,j+1) = mean(xtemp(:,window*(1+n)+1:end),2);
% end

% y = zeros(size(xtemp));
% for i = 1:size(y,1)
%     if nytemp ~= size(xtemp,2)
%         y(i,1:nytemp) = interp(ytemp(i,1:end-1),window);
%         y(i,nytemp:end) = mean([y(i,nytemp-10:nytemp),ytemp(i,end)]);
%     else
%         y(i,:) = interp(ytemp(i,:),window);
%     end
% end
y = ytemp;
if size(x,2)<size(x,1)
    y = y';
end
end