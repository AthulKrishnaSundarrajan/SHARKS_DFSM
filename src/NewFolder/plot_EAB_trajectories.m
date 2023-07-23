clc;clear;close all;

% load results
load('DFSM_FOWT_validation_EAB.mat')
saveflag = false;

C = materialColors;

blue = C.blue(9,:);
yellow = C.yellow(10,:);
red = C.red(7,:);


% plot ptfm pitch
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [134 69 570 413];
commonFigureProperties



t = results_cell{1}{1};
X_cell = results_cell{2};

X_act = X_cell{1};
X_dfsm = X_cell{2};

plot(t,X_dfsm(:,1),'linewidth',linewidth,'color',red)
plot(t,X_act(:,1),'linewidth',linewidth,'color',blue)
plot(t,X_dfsm(:,1)-X_act(:,1),'linewidth',linewidth,'color','k')


xlabel('Time [s]');ylabel('Platform Pitch [deg]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['platform_pitch_',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

% gen speed
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [131 49 570 413];
commonFigureProperties



plot(t,X_dfsm(:,2),'linewidth',linewidth,'color',red)
plot(t,X_act(:,2),'linewidth',linewidth,'color',blue)
plot(t,X_dfsm(:,2)-X_act(:,2),'linewidth',linewidth,'color','k')

xlabel('Time [s]');ylabel('Generator Speed [rad/s]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['gen_speed_',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end


% outputs
Y_cell = results_cell{5};