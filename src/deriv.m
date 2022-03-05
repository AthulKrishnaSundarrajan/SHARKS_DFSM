function dx = deriv(t,x,U,net,maxX,maxY,W_centers)

% get wind speed
u = U(t);

% get wind speed
w = u(1);

% determine which neural net must be used
net_num = 1:length(W_centers);
n = ceil(interp1(W_centers,net_num,w));

% get scaled input values
X = [x',u]./maxX;

% use net
net = net{n};

% evaluate derivative value from neural network
dx = net(X');

% unscale
dx = dx.*maxY';

end