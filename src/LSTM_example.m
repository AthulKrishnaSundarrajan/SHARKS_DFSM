clc; clear; close all;


% load waveform data
waveform_data = load('WaveformData.mat');
data = waveform_data.data;
data1 = data{1};

% get number of channels
n_channels = size(data{1},1);

% visualize the first few seqences
hf = figure;
hf.Color = 'w';

subplot(3,1,1)
plot(data1(1,:))
xlabel("Time Step")

subplot(3,1,2)
plot(data1(2,:))
xlabel("Time Step")

subplot(3,1,3)
plot(data1(3,:))
xlabel("Time Step")

% partition the data
numObservations = numel(data);
idxTrain = 1:floor(0.9*numObservations);
idxTest = floor(0.9*numObservations)+1:numObservations;

dataTrain = data(idxTrain);
dataTest = data(idxTest);

% prepare data for training and testing

for n = 1:numel(dataTrain)

    X = dataTrain{n};
    XTrain{n} = X(:,1:end-1);
    TTrain{n} = X(:,2:end);

end

% normalize

muX = mean(cat(2,XTrain{:}),2);
sigmaX = std(cat(2,XTrain{:}),0,2);

muT = mean(cat(2,TTrain{:}),2);
sigmaT = std(cat(2,TTrain{:}),0,2);

for n = 1:numel(XTrain)
    XTrain{n} = (XTrain{n} - muX) ./ sigmaX;
    TTrain{n} = (TTrain{n} - muT) ./ sigmaT;
end

% define architecture
layers = [
    sequenceInputLayer(n_channels)
    lstmLayer(120)
    fullyConnectedLayer(n_channels)
    regressionLayer];

% training options
options = trainingOptions("adam", ...
    MaxEpochs=200, ...
    SequencePaddingDirection="left", ...
    Shuffle="every-epoch", ...
    Plots="training-progress", ...
    Verbose=1);

% train
net = trainNetwork(XTrain,TTrain,layers,options);

% normalize test dataset
for n = 1:size(dataTest,1)
    X = dataTest{n};
    XTest{n} = (X(:,1:end-1) - muX) ./ sigmaX;
    TTest{n} = (X(:,2:end) - muT) ./ sigmaT;
end

idx = 2;
X = XTest{idx};
T = TTest{idx};

figure
stackedplot(X',DisplayLabels="Channel " + (1:n_channels))
xlabel("Time Step")
title("Test Observation " + idx)

%net = resetState(net);
offset = 75;
[net,~] = predictAndUpdateState(net,X(:,1:offset));

numTimeSteps = size(X,2);
numPredictionTimeSteps = numTimeSteps - offset;
Y = zeros(n_channels,numPredictionTimeSteps);

for t = 1:numPredictionTimeSteps
    Xt = X(:,offset+t);
    [net,Y(:,t)] = predictAndUpdateState(net,Xt);
end


hf = figure;
hf.Color = 'w';

ic = 1;
subplot(3,1,ic); hold on;
plot(T(ic,:))
plot(offset:numTimeSteps,[T(ic,offset) Y(1,:)],'--');
ylabel(['Channel ',num2str(ic)])

ic = ic+1;
subplot(3,1,ic); hold on;
plot(T(ic,:))
plot(offset:numTimeSteps,[T(ic,offset) Y(1,:)],'--');
ylabel(['Channel ',num2str(ic)])

ic = ic+1;
subplot(3,1,ic); hold on;
plot(T(ic,:))
plot(offset:numTimeSteps,[T(ic,offset) Y(1,:)],'--');
ylabel(['Channel ',num2str(ic)])