function[coh,fcoh] = coherence(varargin)

x = varargin{1};
y = varargin{2};
fs = varargin{3};

if length(varargin) > 2
    window = varargin{4};
    overlap = varargin{5};
    nseg = varargin{6};
end

N = length(x);

if length(varargin) > 2
    
    sampseg = floor(length(x)/nseg);
    totseg = floor(length(x)/(sampseg-overlap));
    
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