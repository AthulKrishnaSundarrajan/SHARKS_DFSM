clc; clear; close all;

% load results
load('DFSM_oloc_results_belowrated.mat')

fol_name = 'plot_oloc_results/transition';
saveflag = false;


ntests = length(time_cell);
linewidth = 0.75;

plot_inds = [1,2];

C = materialColors;

blue = C.blue(7,:);
yellow = C.yellow(10,:);
red = C.red(7,:);

color_map = cell(ntests,1);
color_map{plot_inds(1)} = blue;
color_map{plot_inds(2)} = yellow;


%---------------------------------------------------
% plot wind speed
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

for i = plot_inds

U = U_cell{i};
plot(time_cell{i},U(:,1),'linewidth',linewidth,'color',color_map{i})

end

xlabel('Time [s]');ylabel('Wind Speed [m/s]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = 'wind_speed';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

%-----------------------------------------------------------------

% plot generator torque
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

for i = plot_inds

U = U_cell{i};
plot(time_cell{i},U(:,2),'linewidth',linewidth,'color',color_map{i})
ylim([19,20])

end

xlabel('Time [s]');ylabel('Generator Torque [MWm]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = 'gen_torque';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end
%--------------------------------------------------------------
% plot blade pitch
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

for i = plot_inds

U = U_cell{i};
plot(time_cell{i},U(:,3),'linewidth',linewidth,'color',color_map{i})

end

xlabel('Time [s]');ylabel('Blade Pitch [deg]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = 'blade_pitch';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end
%---------------------------------------------------------
% plot platform pitch
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

for i = plot_inds

X = X_cell{i};
plot(time_cell{i},X(:,1),'linewidth',linewidth,'color',color_map{i})

end

xlabel('Time [s]');ylabel('Platform Pitch [deg]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = 'ptfm_pitch';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end
%-------------------------------------------------------
% plot generator speed
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

for i = plot_inds

X = X_cell{i};
plot(time_cell{i},X(:,2),'linewidth',linewidth,'color',color_map{i})
ylim([6.5,7.5])

end

xlabel('Time [s]');ylabel('Generator Speed [rpm]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = 'gen_speed';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end
%------------------------------------------------------
% plot generator power
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

for i = plot_inds

X = X_cell{i};

U = U_cell{i};

power = U(:,2).*(X(:,2).*0.1047198*0.9941);
plot(time_cell{i},power,'linewidth',linewidth,'color',color_map{i})
ylim([14.5,15.5])

end

xlabel('Time [s]');ylabel('Generator Power [MW]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = 'gen_power';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
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
legend_entries = {'Case 1','Case 2'};
% legend
legend(legend_entries);
fontlegend = 22; nCol = 3;lcn = 'best';
taskflag = 'legend';commonFigureTasks;

% export
if saveflag
    savename = 'legend_common';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end