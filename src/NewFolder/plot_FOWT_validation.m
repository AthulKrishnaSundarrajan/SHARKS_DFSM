clc;clear; close all;

% load result
load('DFSM_FOWT_validation_QR.mat')

C = materialColors;

blue = C.blue(9,:);
yellow = C.yellow(10,:);
red = C.red(7,:);

color_map = cell(2,1);
color_map{1} = blue;
color_map{2} = red;
color_map{3} = yellow;


fol_name = 'plot_oloc_results/validation/QR';
saveflag = ~false;

%% plot wind speed



for i = 1
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties


t = results_cell{1}{1};
U = results_cell{1}{2};

plot(t,U(:,i),'linewidth',linewidth,'color',blue)

xlabel('Time [s]');ylabel('Wind Speed [m/s]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['wind_speed_',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

end



%% plot platform pitch


for i = 1

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [1314 469 570 413];
    commonFigureProperties
    
    
    
    t = results_cell{i}{1};
    X_cell = results_cell{i}{3};
    
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


end


%% plot gen speed


for i = 1


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties


t = results_cell{i}{1};
X_cell = results_cell{i}{3};

X_act = X_cell{1};
X_dfsm = X_cell{2};


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


end

%%
for i = 1


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties


t = results_cell{i}{1};
Y_cell = results_cell{i}{4};

Y_act = Y_cell{1};
Y_dfsm = Y_cell{2};


    plot(t,Y_dfsm(:,1),'linewidth',linewidth,'color',red)
    plot(t,Y_act(:,1),'linewidth',linewidth,'color',blue)
    plot(t,Y_dfsm(:,1)-Y_act(:,1),'linewidth',linewidth,'color','k')

xlabel('Time [s]');ylabel('TwrBsFxt [kN]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['TwrBsFxt',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

%-----------------------------------------------------

% hf = figure;
% hf.Color = 'w';
% hold on;
% hf.Position = [1314 469 570 413];
% commonFigureProperties
% 
% plot(t,Y_dfsm(:,1)-Y_act(:,1),'linewidth',linewidth)
% 
% xlabel('Time [s]');ylabel('Error [kN]');
% 
% %xlabel('Time [s]');ylabel('TwrBsFxt [kN]');
% taskflag = 'axes'; commonFigureTasks;
% 
% if saveflag
%     savename = ['TwrBsFxt_hist',num2str(i)];
%     pathpdf = mfoldername(mfilename('fullpath'),fol_name);
%     filename = fullfile(pathpdf,savename);
%     str = strcat("export_fig '",filename,"' -pdf");
%     eval(str)
% end
% 

end
%%
for i = 1


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties


t = results_cell{i}{1};
Y_cell = results_cell{i}{4};

Y_act = Y_cell{1};
Y_dfsm = Y_cell{2};


plot(t,Y_dfsm(:,4),'linewidth',linewidth,'color',red)
plot(t,Y_act(:,4),'linewidth',linewidth,'color',blue)
plot(t,Y_dfsm(:,4)-Y_act(:,4),'linewidth',linewidth,'color','k')

xlabel('Time [s]');ylabel('TwrBsMxt [kNm]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['TwrBsMxt',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

%-----------------------------------------------------
% hf = figure;
% hf.Color = 'w';
% hold on;
% hf.Position = [1314 469 570 413];
% commonFigureProperties
% 
% plot(t,Y_dfsm(:,4)-Y_act(:,4),'linewidth',linewidth)
% 
% xlabel('Time [s]');ylabel('Error [kNm]');
% taskflag = 'axes'; commonFigureTasks;
% 
% if saveflag
%     savename = ['TwrBsMxt_hist',num2str(i)];
%     pathpdf = mfoldername(mfilename('fullpath'),fol_name);
%     filename = fullfile(pathpdf,savename);
%     str = strcat("export_fig '",filename,"' -pdf");
%     eval(str)
% end


end

%%
for i = 1


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties


t = results_cell{i}{1};
Y_cell = results_cell{i}{4};

Y_act = Y_cell{1};
Y_dfsm = Y_cell{2};


plot(t,Y_dfsm(:,3),'linewidth',linewidth,'color',red)
plot(t,Y_act(:,3),'linewidth',linewidth,'color',blue)
plot(t,Y_dfsm(:,3)-Y_act(:,3),'linewidth',linewidth,'color','k')

xlabel('Time [s]');ylabel('RtAeroCt [-]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['rt_aero_ct_',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end
%--------------------------------------------------
% hf = figure;
% hf.Color = 'w';
% hold on;
% hf.Position = [1314 469 570 413];
% commonFigureProperties
% 
% plot(t,Y_dfsm(:,3)-Y_act(:,3),'linewidth',linewidth)
% 
% xlabel('Time [s]');ylabel('Error');
% taskflag = 'axes'; commonFigureTasks;
% 
% if saveflag
%     savename = ['rt_aero_cp_hist',num2str(i)];
%     pathpdf = mfoldername(mfilename('fullpath'),fol_name);
%     filename = fullfile(pathpdf,savename);
%     str = strcat("export_fig '",filename,"' -pdf");
%     eval(str)
% end



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
legend_entries = {'DFSM','WEIS','Error'};
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
