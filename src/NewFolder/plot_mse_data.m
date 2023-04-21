clc;clear;close all;

% load data file
load('mse-data.mat');

% find average
mean_dx_linfit = squeeze(mean(mse_dx_linfit,2));
max_dx_linfit = squeeze(max(mse_dx_linfit,[],2));
min_dx_linfit = squeeze(min(mse_dx_linfit,[],2));

mean_dx_taylor = squeeze(mean(mse_dx_taylor,2));
max_dx_taylor = squeeze(max(mse_dx_taylor,[],2));
min_dx_taylor = squeeze(min(mse_dx_taylor,[],2));

hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

markersize = 15;
linewidth = 1.5;

plot(fac,mean_dx_taylor(:,1),'r.-','linewidth',linewidth,'markersize',markersize)
plot(fac,mean_dx_linfit(:,1),'k.-','linewidth',linewidth,'markersize',markersize)

xlabel('k');ylabel('MSE')

ha = gca;
ha.XScale = "log";
%ha.YScale = "log";
grid on

taskflag = 'axes';commonFigureTasks;


legend('taylor','linfit')
%-----------------------------------
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];


for i = 1:length(fac)

    input = inputs_cell{i};

    max_input(i,:) = max(input,[],1) ;

    min_input(i,:) = min(input,[],1);

end

plot(fac,max_input,'k.-','linewidth',linewidth,'markersize',markersize)
plot(fac,min_input,'k.-','linewidth',linewidth,'markersize',markersize)
xlabel('k');ylabel('Range')

ha = gca;
ha.XScale = "log";
grid on
taskflag = 'axes';commonFigureTasks;


return