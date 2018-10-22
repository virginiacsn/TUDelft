x = 8*ones(1,800);
xlong = 10*ones(1,5000);
[b,a] = butter(2,5/(1024/2),'low');
y = filter(b,a,x);
ylong = filter(b,a,xlong);
mean(y)
mean(ylong)

figure;
subplot(1,2,1)
plot(y)
subplot(1,2,2)
plot(ylong)