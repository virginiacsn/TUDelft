function[coh,fcoh,totseg] = coherence(varargin)

x = varargin{1};
y = varargin{2};
fs = varargin{3};

if length(varargin) > 2
    window = varargin{4};
    overlap = varargin{5};
    nseg = varargin{6};
    if length(varargin) > 6
        iapp = varargin{7};
    end
end

if length(varargin) > 2 && length(varargin) < 6
    
    N = length(x);
    sampseg = floor(N/nseg);
    
    if sampseg <= overlap
        error('Overlap should be less than samples per segement. Reduce number of segments.')
    else
        totseg = floor(N/(sampseg-overlap));
    end
    
    X = zeros(sampseg,totseg-1);
    Y = zeros(sampseg,totseg-1);
    
    for i = 0:totseg-2
        X(:,i+1) = fft(x(1+i*(sampseg-overlap):i*(sampseg-overlap)+sampseg)*window(sampseg));
        Y(:,i+1) = fft(y(1+i*(sampseg-overlap):i*(sampseg-overlap)+sampseg)*window(sampseg));
    end
    
%     X = fft(reshape(x,[sampseg,nseg]));
%     Y = fft(reshape(x,[sampseg,nseg]));
    
    Sxx = mean(nseg/sampseg*real(conj(X).*X),2);
    Syy = mean(nseg/sampseg*real(conj(Y).*Y),2);
    Sxy = mean(nseg/sampseg*conj(X).*Y,2);
  
elseif length(varargin) > 6
   
   N = min(diff([find(iapp); length(iapp)])); 
   napp = sum(iapp);
   sampseg = floor(N/nseg);
   overlap = floor(sampseg/2);
   
    if sampseg <= overlap
        error('Overlap should be less than samples per segement. Reduce number of segments.')
    else
        totseg = floor(N/(sampseg-overlap));
    end
    
    count = 1;
    for j = 0:napp-1
        X = zeros(sampseg,totseg-1);
        Y = zeros(sampseg,totseg-1);
        
        for i = 0:totseg-2
            X(:,count) = fft(x(1+i*(sampseg-overlap):i*(sampseg-overlap)+sampseg)*window(sampseg));
            Y(:,count) = fft(y(1+i*(sampseg-overlap):i*(sampseg-overlap)+sampseg)*window(sampseg));
            count = count+1;
        end
    end
    
    Sxx = mean(nseg/sampseg*real(conj(X).*X),2);
    Syy = mean(nseg/sampseg*real(conj(Y).*Y),2);
    Sxy = mean(nseg/sampseg*conj(X).*Y,2);
    
else
    X = fft(x);
    Y = fft(y);
    
    Sxx = 1/N*X.*conj(X);
    Syy = 1/N*Y.*conj(Y);
    Sxy = 1/N*X.*conj(Y);
end

coh = abs(Sxy).^2./(Sxx.*Syy);
fcoh = (0:length(coh)-1)/(length(coh)/fs);

end