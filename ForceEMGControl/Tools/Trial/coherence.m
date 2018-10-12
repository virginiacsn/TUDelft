function[coh,fcoh,CLseg] = coherence(x,y,fs,win,overlap,iapp)

if isempty(iapp)
    
    N = length(x);
    sampseg = length(win);
    
    if sampseg <= overlap
        error('Overlap should be less than samples per segement. Reduce number of segments.')
    else
        totseg = floor(N/(sampseg-overlap));
    end
    
    X = zeros(sampseg,totseg-1);
    Y = zeros(sampseg,totseg-1);
    
    for i = 0:totseg-2
        X(:,i+1) = fft(x(1+i*(sampseg-overlap):i*(sampseg-overlap)+sampseg).*win);
        Y(:,i+1) = fft(y(1+i*(sampseg-overlap):i*(sampseg-overlap)+sampseg).*win);
    end
    
    %     X = fft(reshape(x,[sampseg,nseg]));
    %     Y = fft(reshape(x,[sampseg,nseg]));
   
    % WHAT?
    Sxx = mean(1/sampseg*real(conj(X).*X),2);
    Syy = mean(1/sampseg*real(conj(Y).*Y),2);
    Sxy = mean(1/sampseg*conj(X).*Y,2);
    
    CLseg = totseg;
    
elseif ~isempty(iapp)
    
    N = min(diff(find(iapp)));
    idxapp = find(iapp);
    napp = sum(iapp);
    sampseg = length(win);%floor(N/nseg);
    overlap = floor(sampseg/2);
    
    if sampseg <= overlap
        error('Overlap should be less than samples per segement. Reduce number of segments.')
    else
        totseg = floor(N/(sampseg-overlap));
    end
        
    X = zeros(sampseg,totseg*napp);
    Y = zeros(sampseg,totseg*napp);
    count = 1;
    
        for j = 1:napp
            for i = 0:totseg-2
                X(:,count) = fft(x(idxapp(j)+(sampseg-overlap)*i:idxapp(j)+(sampseg-overlap)*i+sampseg-1).*win);
                Y(:,count) = fft(y(idxapp(j)+(sampseg-overlap)*i:idxapp(j)+(sampseg-overlap)*i+sampseg-1).*win);
                count = count+1;
            end
        end
    
    Sxx = mean(1/sampseg*real(conj(X).*X),2);
    Syy = mean(1/sampseg*real(conj(Y).*Y),2);
    Sxy = mean(1/sampseg*conj(X).*Y,2);
    
    CLseg = totseg*napp;
    
elseif isempty(win)
    N = length(x);
    
    X = fft(x);
    Y = fft(y);
    
    Sxx = 1/N*X.*conj(X);
    Syy = 1/N*Y.*conj(Y);
    Sxy = 1/N*X.*conj(Y);
    
    CLseg = 1;
end

coh = abs(Sxy).^2./(Sxx.*Syy);
coh_avg = movingAvg(coh,3);
fcoh = (0:length(coh_avg)-1)/(length(coh_avg)/fs);

end

    %     if napp == 1
    %         X(:,1) = fft(x(idxapp(1):idxapp(1)+N-1));
    %         Y(:,1) = fft(y(idxapp(1):idxapp(1)+N-1));
    %     else
    %         for j = 1:napp-1
    %             %         for i = 0:totseg-1
    %             %             X(:,count) = fft(x(idxapp(j)+(sampseg-overlap)*i:idxapp(j)+(sampseg-overlap)*i+sampseg-1)*window(sampseg));
    %             %             Y(:,count) = fft(y(idxapp(j)+(sampseg-overlap)*i:idxapp(j)+(sampseg-overlap)*i+sampseg-1)*window(sampseg));
    %             %             count = count+1;
    %             %         end
    %
    %             X(:,j) = fft(x(idxapp(j):idxapp(j)+N-1));
    %             Y(:,j) = fft(y(idxapp(j):idxapp(j)+N-1));
    %         end
    %     end
%     if napp == 1
%         for i = 0:totseg-2
%             X(:,count) = fft(x(idxapp(1)+(sampseg-overlap)*i:idxapp(1)+(sampseg-overlap)*i+sampseg-1).*window(sampseg));
%             Y(:,count) = fft(y(idxapp(1)+(sampseg-overlap)*i:idxapp(1)+(sampseg-overlap)*i+sampseg-1).*window(sampseg));
%             count = count+1;
%         end
%     else