function[y] = movingavg(x,window)

if size(x,2)<size(x,1)
    xtemp = x';
else
    xtemp = x;
end

j = 0;
ytemp = zeros(size(xtemp,1),floor(size(xtemp,2)/window)-1);
for n = 0:floor(size(xtemp,2)/window)-1
    j = j+1;
    ytemp(:,j) = mean(xtemp(:,1+window*n:window*(1+n)),2);
end

y = zeros(size(xtemp,1),window*(1+n));
for i = 1:size(y,1)
    y(i,:) = interp(ytemp(i,:),window);
end

if size(x,2)<size(x,1)
    y = y';
end
end

