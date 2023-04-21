clc; clear; close all;

% set seed
rng(4357)

load('TLR_simulations.mat');

nsamples = length(dx_cell_lin);

ind = randsample(1:nsamples,1);


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

plot(time,controls(:,1),'-','linewidth',1.5)
plot(time,controls(:,2),'-','linewidth',1.5)
xlabel('Time [s]')
ylabel('$u$')
taskflag = 'axes';commonFigureTasks

%-----------------------------------------
t0 = time(1); tf = time(end);

nt = length(dx_cell_nl{1,1});
time_ = linspace(t0,tf,nt);

dx_act = dx_cell_lin{1,1};
dx_lin = dx_cell_lin{1,2};
dx_mf = dx_cell_mf{1,2};
dx_nl = dx_cell_nl{1,2};

X_act = X_cell_lin{1,1};
X_mf = X_cell_mf{1,2};
X_nl = X_cell_nl{1,2};

time_x = linspace(t0,tf,length(X_mf));

blue = C.blue(9,:);
yellow = C.amber(8,:);
red = C.red(9,:);
green = C.green(9,:);
orange = C.deeporange(6,:);



for idx = 1:4

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [1314 469 570 413];
    commonFigureProperties

    plot(time_,dx_act(:,idx),'linewidth',1.5,'color',blue)
    plot(time_,dx_lin(:,idx),'linewidth',1.5,'color',yellow)
    xlabel('Time [s]');ylabel(['$\dot{\xi}_',num2str(idx),'$']);

    taskflag = 'axes';commonFigureTasks



end

for idx = 1:4

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [1314 469 570 413];
    commonFigureProperties

    plot(time_,dx_act(:,idx),'linewidth',1.5,'color',blue)
    plot(time_,dx_mf(:,idx),'linewidth',1.5,'color',red)
    xlabel('Time [s]');ylabel(['$\dot{\xi}_',num2str(idx),'$']);

    taskflag = 'axes';commonFigureTasks



end
%close all;
for idx = 1:4

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [1314 469 570 413];
    commonFigureProperties

    plot(time_x,X_act(:,idx),'linewidth',1.5,'color',blue)
    plot(time_x,X_mf(:,idx),'linewidth',1.5,'color',red)
    plot(time_x,X_nl(:,idx),'linewidth',1.5,'color',yellow)
    xlabel('Time [s]');ylabel(['$\xi_',num2str(idx),'$']);

    taskflag = 'axes';commonFigureTasks



end
