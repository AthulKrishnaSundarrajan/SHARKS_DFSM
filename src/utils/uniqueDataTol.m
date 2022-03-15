%
function [X,Y] = uniqueDataTol(X,Y,tol)

% combine
xy = [X,Y];

% determine unique values (lowest first)
[~,IA_lowest,~] = uniquetol(xy,tol,'lowest','ByRows',true);
I_interior_lowest = false(size(xy,1),1);
I_interior_lowest(IA_lowest) = true;

% determine unique values (highest first)
[~,IA_highest,~] = uniquetol(xy,tol,'highest','ByRows',true);
I_interior_highest = false(size(xy,1),1);
I_interior_highest(IA_highest) = true;

% combine indices
I = I_interior_highest | I_interior_lowest;

% extract unique data points
X = X(I,:);
Y = Y(I,:);

end