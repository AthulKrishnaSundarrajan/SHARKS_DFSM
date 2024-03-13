function V = estimate_loss(X_m,X_e,W)

% get the length
N = length(X_m);

% initialize
V = zeros(N,1);

% loop through and estimate error
for i = 1:N
    
    % error
    e = (X_m(i,:) - X_e(i,:))';
    
    % loss
    V(i) = e'*W*e;
end

% mean loss
V = mean(V);


end