function[cor,lags] = correlation(x,y,iapp)

if length(find(iapp))>1
    N = min(diff(find(iapp)));
    
    idxapp = find(iapp);
    napp = sum(iapp);
    
    %totseg = floor(N/napp);
    
    % X = zeros(N,napp);
    % Y = zeros(N,napp);
    X = zeros(N*2-1,napp);
    
    for i = 1:napp-1
        %     X(:,i) = x(idxapp(i):idxapp(i+1)-1);
        %     Y(:,i) = y(idxapp(i):idxapp(i+1)-1);
        [X(:,i),lags] = xcorr(x(idxapp(i):idxapp(i+1)-1),y(idxapp(i):idxapp(i+1)-1),'unbiased');
    end
    cor = mean(X,2)./mean(X(N,:));
else
    [cor,lags] = xcorr(x,y,'unbiased');
end

% xcor = xcorr(mean(X,2),'unbiased');
% ycor = xcorr(mean(Y,2),'unbiased');

% [cor,lags] = xcorr(mean(X,2),mean(Y,2),'unbiased');
% cor = cor./(xcor(N)*ycor(N));
end
