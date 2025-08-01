function[xBound,yBound] = getBound(x,y,CI,outType)

upR = y+CI;
lowR = y-CI;

if strcmp(outType,'cart')
    [upX,upY] = pol2cart(x,upR);
    [lowX,lowY] = pol2cart(x,lowR);
    
    if isrow(upX)
        xBound = [upX fliplr(lowX)];
        yBound = [upY fliplr(lowY)];
    else
        xBound = [upX; flipud(lowX)];
        yBound = [upY; flipud(lowY)];
    end
elseif strcmp(outType,'polar')
    xBound = x;
    if isrow(upR)
        yBound = [upR; lowR];
    else
        yBound = [upR lowR]';
    end
end