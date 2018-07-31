function[coh] = coherence(varargin)

x = varargin{1};
y = varargin{2};

if length(varargin) > 2
    window = varargin{3};
    overlap = varargin{4};
    nseg = varargin{5};
end

N = length(x);

X = fft(x);
Y = fft(y);

if length(varargin) > 2
    
end
Sxx = 1/N*X.*conj(X);
Syy = 1/N*Y.*conj(Y);
Sxy = 1/N*X.*conj(Y);

coh = (abs(Sxy).^2)./(Sxx.*Syy);

end