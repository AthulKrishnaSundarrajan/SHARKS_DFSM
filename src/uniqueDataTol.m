%
function [X,Y] = uniqueDataTol(X,Y,tol)

% determine unique values
[~,IA,~] = uniquetol([X,Y],tol,'ByRows',true);

% extract unique indices
X = X(IA,:);
Y = Y(IA,:);

end