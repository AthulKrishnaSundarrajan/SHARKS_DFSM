%
function [X,Y] = uniqueDataTol(X,Y,tol)

xy = [X,Y];

% tess = convhulln(xy);
% in = inhull(xy,xy,tess);

% I_convex = false(size(xy,1),1);
% I_convex(unique(tess)) = true;


% determine unique values (lowest first)
[~,IA_lowest,~] = uniquetol(xy,tol,'lowest','ByRows',true);
I_interior_lowest = false(size(xy,1),1);
I_interior_lowest(IA_lowest) = true;

% determine unique values (highest first)
[~,IA_highest,~] = uniquetol(xy,tol,'highest','ByRows',true);
I_interior_highest = false(size(xy,1),1);
I_interior_highest(IA_highest) = true;

I = I_interior_highest | I_interior_lowest;

% extract unique indices
X = X(I,:);
Y = Y(I,:);

end
