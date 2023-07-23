clc;clear;close all;

% load results

load("DFSM_MHK_validation_EAB.mat")
saveflag = ~false;
fol_name = 'plot_oloc_results/validation/EAB_MHK2';

C = materialColors;

blue = C.blue(9,:);
yellow = C.yellow(10,:);
red = C.red(7,:);
grey = C.grey(7,:);

% extract

time = results_cell{1};
input_test = results_cell{2};
X_cell = results_cell{3};

% plot controls
nu = size(input_test,2);

ylabel_cell = {'Wave Elevation [m]','Current Speed [m/s]'};
savename_cell = {'wave','current'};

for i = 1:nu

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [131 49 570 413];
    commonFigureProperties
    
    
    plot(time,input_test(:,i),'linewidth',linewidth,'color',blue)
    
    xlabel('Time [s]');ylabel(ylabel_cell{i});
    taskflag = 'axes'; commonFigureTasks;
    
    if saveflag
        savename = [savename_cell{i},'test'];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -svg");
        eval(str)
    end


end

for i = 1:nu

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [131 49 570 413];
    commonFigureProperties
    
    
    plot(time,input_train(:,i),'linewidth',linewidth,'color',blue)
    
    xlabel('Time [s]');ylabel(ylabel_cell{i});
    taskflag = 'axes'; commonFigureTasks;
    
    if saveflag
        savename = [savename_cell{i},'train'];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -svg");
        eval(str)
end


end

% states
nx = size(X_cell{1},2);
ylabel_cell =  {'Surge [m]','Sway [m]','Heave [m]','Roll [deg]','Pitch [deg]','Yaw [deg]'};
savename_cell =  {'PtfmSurge','PtfmSway','PtfmHeave','PtfmRoll','PtfmPitch','PtfmYaw'};

for i = 1:nx

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [131 49 570 413];
    commonFigureProperties

    X_act = X_cell{1};
    X_dfsm = X_cell{2};
    
    
    plot(time,X_dfsm(:,i),'linewidth',linewidth,'color',red)
    plot(time,X_act(:,i),'linewidth',linewidth,'color',blue)
    plot(time,X_dfsm(:,i)-X_act(:,i),'-','linewidth',linewidth-0.5,'color',grey)

    
    xlabel('Time [s]');ylabel(ylabel_cell{i});
    taskflag = 'axes'; commonFigureTasks;
    
    if saveflag
        savename = [savename_cell{i}];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -svg");
        eval(str)
    end


end

%------------------------------------------------------------------
hf2 = copyfig(hf);
hf2.Position = [263.4000 417.8000 500.2000 70.4000];
h = findall(hf2,'type','line');
close(hf)

axis off; 
for i = 1:length(h)

    h(i).XData = NaN; %ignore warnings

end
legend_entries = {'DFSM','OpenFAST','Error'};
% legend
legend(legend_entries);
fontlegend = 22; nCol = 3;lcn = 'best';
taskflag = 'legend';commonFigureTasks;

% export
if saveflag
    savename = 'legend_common';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -svg");
    eval(str)
end