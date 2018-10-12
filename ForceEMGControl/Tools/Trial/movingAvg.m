function[y] = movingAvg(x,window)

if size(x,2)<size(x,1)
    xtemp = x';
else
    xtemp = x;
end

j = 0;
ytemp = zeros(size(xtemp,1),floor(size(xtemp,2)/window));
for n = 0:floor(size(xtemp,2)/window)-1
    j = j+1;
    ytemp(:,j) = mean(xtemp(:,1+window*n:window*(1+n)),2);
end

nytemp = floor(size(xtemp,2)/window)*window;
if nytemp ~= size(xtemp,2)
    ytemp(:,j+1) = mean(xtemp(:,window*(1+n)+1:end),2);
end

y = zeros(size(xtemp));
for i = 1:size(y,1)
    if length(ytemp(i,:)) > 9
        if nytemp ~= size(xtemp,2)
            y(i,1:nytemp) = interp(ytemp(i,1:end-1),window);
            y(i,nytemp:end) = mean([y(i,nytemp-10:nytemp),ytemp(i,end)]);
        else
            y(i,:) = interp(ytemp(i,:),window);
        end
    end
end

if size(x,2)<size(x,1)
    y = y';
end
end

