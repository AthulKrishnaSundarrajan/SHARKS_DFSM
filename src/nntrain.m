%
function net = nntrain(X,Y)

% transpose input data
x = X';
t = Y';

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
trainFcn = 'trainbr';  % Levenberg-Marquardt backpropagation.

% Create a Fitting Network
% hiddenLayerSize = 10;
% net = fitnet(hiddenLayerSize,trainFcn);

hiddenLayerSize1 = 10;
hiddenLayerSize2 = 10;
net = fitnet([hiddenLayerSize1,hiddenLayerSize2],trainFcn);

net.trainParam.max_fail = 10000;
net.trainParam.mu_max = 1e20;
net.trainParam.min_grad = 1e-8;
net.trainParam.epochs = 5000;

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 20/100;
net.divideParam.valRatio = 10/100;
net.divideParam.testRatio = 70/100;
net.divideFcn = 'divideint';

% Train the Network
[net,tr] = train(net,x,t,'UseParallel','yes');

% Test the Network
% y = net(x);
% e = gsubtract(t,y);
% performance = perform(net,t,y);

% View the Network
% view(net)

end