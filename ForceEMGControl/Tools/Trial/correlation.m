function[cor,lags] = correlation(x,y,iapp)

if length(find(iapp))>1
    N = min(diff(find(iapp)));
else
    N = length(x);
end

idxapp = find(iapp);
napp = sum(iapp);

%totseg = floor(N/napp);

X = zeros(N*2-1,napp);

for i = 1:napp-1
    [X(:,i),lags] = xcorr(x(idxapp(i):idxapp(i+1)-1),y(idxapp(i):idxapp(i+1)-1),'unbiased');
end

cor = mean(X,2);
end