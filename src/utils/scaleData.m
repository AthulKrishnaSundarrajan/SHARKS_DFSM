%
function [X,Y,maxX,maxY] = scaleData(X,Y)

% maximum column values
maxX = max(abs(X),[],1);
maxY = max(abs(Y),[],1);

% avoid zeros
maxX(maxX==0) = 1;
maxY(maxY==0) = 1;

% scale
X = X./maxX;
Y = Y./maxY;

end