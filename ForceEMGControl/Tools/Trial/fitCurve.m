function[xfit,yfit] = fitCurve(x,y,order)
coef = polyfit(x,y,order);
xfit = min(x):5:max(x);
yfit = polyval(coef,xfit);
end