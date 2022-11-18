function [Y,Xf,Af] = DFSM_TEST_net(X,~,~)
%DFSM_TEST_NET neural network simulation function.
%
% Auto-generated by MATLAB, 28-Oct-2022 12:01:54.
% 
% [Y] = DFSM_TEST_net(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 3xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 2xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-4.36139346181444;-6.40734442251825;-4.99999994854441];
x1_step1.gain = [0.243040300991632;0.142831682022955;0.200000001659478];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.7995295872486094524;2.419991744318879423;1.6347951139213678928;-3.017308764160606227;-1.0853218033353755345;-0.0078028471949377602862;-0.62170073318032592624;0.19189881586639956357;1.907464720224855359;1.2361062888448770547;-0.092339567235242864252;-1.5687728714119224538;1.8237331081062511107;-1.8257577796264687908;4.5024757017773806567];
IW1_1 = [0.78691217839765659203 -0.61030168178059851503 -0.00013093862337849860016;-1.3615520348394536398 -1.1399482156128351029 0.0020552081533085227692;-0.51019390225530714122 -1.5408015223591371434 -1.9148762328994757276;3.4590555870044830122 6.1148012903847872579 0.077691290703941151552;1.1609269814935303522 1.0107974386325115912 -0.0050868263453261531107;0.1221580004066792069 -0.71655133149187555386 -0.020534015209143736264;3.2800191308079815755 2.8673149640210602307 -0.41243494570569949609;0.11344067369202995299 0.62046755276877330054 -0.018733932625370645908;0.055791738546154094613 1.1129977254043670332 0.013806699642607353148;-0.2021391674316619369 -2.508478274994939472 0.14586193271157579798;0.22885797057092857609 -1.4579758406691476491 -0.04979431739883890401;-0.7462123386138738157 0.70976719379090824624 -0.00071209727595503093025;0.81554893350043544675 0.79030103978564747447 0.0017809063173774108794;-0.15405673053684521001 1.2742375672797801034 0.012603669009005396759;0.87588566823347357992 1.6729453913186307368 -1.8301506364091943979];

% Layer 2
b2 = [-0.25593427332594309576;2.6327605272645060452];
LW2_1 = [-0.10101711425792052335 -0.029547030465074272665 -2.7098139390856987871e-05 -0.00080525443638853125283 0.017639015487076133842 -0.6610602355378031314 -0.00041452571978016413421 0.71219772783746815303 0.53474239096367059965 -0.020423713325113616174 0.019418433859820286719 -0.042535642597429736134 -0.043139235707379992857 0.45811673158614363288 0.006432284895188370398;10.997543722222726359 4.5087682211483626205 0.00043410289900961963579 0.00047072554688554982955 -1.062506455703306818 -5.6033425174160500504 -0.00068650794315365591507 -2.3055028527502430258 7.0675723331118032178 -0.0062160686743665009499 0.38395000400526313467 -11.123879495112552362 -9.4418443782209937609 5.5929149053723135054 -0.0042274458601894073645];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = [0.142831243967713;0.0308924061633254];
y1_step1.xoffset = [-6.40736420566931;-23.639246214651];

% ===== SIMULATION ========

% Format Input Arguments
isCellX = iscell(X);
if ~isCellX
  X = {X};
end

% Dimensions
TS = size(X,2); % timesteps
if ~isempty(X)
  Q = size(X{1},2); % samples/series
else
  Q = 0;
end

% Allocate Outputs
Y = cell(1,TS);

% Time loop
for ts=1:TS

    % Input 1
    Xp1 = mapminmax_apply(X{1,ts},x1_step1);
    
    % Layer 1
    a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*Xp1);
    
    % Layer 2
    a2 = repmat(b2,1,Q) + LW2_1*a1;
    
    % Output 1
    Y{1,ts} = mapminmax_reverse(a2,y1_step1);
end

% Final Delay States
Xf = cell(1,0);
Af = cell(2,0);

% Format Output Arguments
if ~isCellX
  Y = cell2mat(Y);
end
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
  x = bsxfun(@minus,y,settings.ymin);
  x = bsxfun(@rdivide,x,settings.gain);
  x = bsxfun(@plus,x,settings.xoffset);
end