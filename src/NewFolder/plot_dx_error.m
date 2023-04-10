clc; clear; close all;

% load results
load('TLR_dx_error.mat');

dx_error_nl = vertcat(dx_error_nl{:});
dx_error_mf = vertcat(dx_error_mf{:});

hf = figure;
hf.Color = 'w';
hold on;
commonFigureProperties

idx = 1;
histogram(dx_error_nl(:,idx),'FaceAlpha',0.7)
histogram(dx_error_mf(:,idx),'FaceAlpha',0.7)
legend('nonlinear','multifidelity')
%----------------------------------------------------
hf = figure;
hf.Color = 'w';
hold on;
commonFigureProperties

idx = idx+1;
histogram(dx_error_nl(:,idx),'FaceAlpha',0.7)
histogram(dx_error_mf(:,idx),'FaceAlpha',0.7)
legend('nonlinear','multifidelity')
%---------------------------------------------------
hf = figure;
hf.Color = 'w';
hold on;
commonFigureProperties

idx = idx+1;
histogram(dx_error_nl(:,idx),'FaceAlpha',0.7)
histogram(dx_error_mf(:,idx),'FaceAlpha',0.7)
legend('nonlinear','multifidelity')
%----------------------------------------------------
hf = figure;
hf.Color = 'w';
hold on;
commonFigureProperties

idx = idx+1;
histogram(dx_error_nl(:,idx),'FaceAlpha',0.7)
histogram(dx_error_mf(:,idx),'FaceAlpha',0.7)
legend('nonlinear','multifidelity')
