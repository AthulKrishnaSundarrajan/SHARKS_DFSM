clc;clear; close all;

% load result
load('DFSM_MHK_TR_validation_Final.mat')

C = materialColors;
linewidth  = 1;

n = 4;

blue = C.blue(9,:);
yellow = C.yellow(10,:);
red = C.red(7,:);
grey = C.grey(7,:);

color_map = cell(2,1);
color_map{1} = blue;
color_map{2} = red;
color_map{3} = yellow;


fol_name = 'plot_oloc_results_MHK/TR/Final';
saveflag = ~false;

%% plot wind speed

ww = 0.087;

for i = n
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [131 49 570 413];
commonFigureProperties

results_cell_ = results_cell{i};


t = results_cell_{1};
t = t - 100;
U = results_cell_{2};

plot(t,U(:,1),'linewidth',linewidth,'color',blue)
xlim([0,600])
xlabel('Time [s]');ylabel('Current Speed [m/s]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['curr_speed_',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

end

for i = n
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [131 49 570 413];
commonFigureProperties

results_cell_ = results_cell{i};


t = results_cell_{1};
t = t - 100;
U = results_cell_{2};

plot(t,U(:,4),'linewidth',linewidth,'color',blue)
xlim([0,600])
xlabel('Time [s]');ylabel('Wave Elevation [m]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['wave_elev_',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end


end




%% plot platform pitch


for i = n

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [314 69 570 413];
    commonFigureProperties
    
    results_cell_ = results_cell{i};
    
    t = results_cell_{1};
    t = t - 100;
    X_cell = results_cell_{3};
    
    X_act = X_cell{1};
    X_dfsm = X_cell{2};
    
    plot(t,X_dfsm(:,1),'linewidth',linewidth,'color',red)
    plot(t,X_act(:,1),'linewidth',linewidth,'color',blue)
    plot(t,X_dfsm(:,1)-X_act(:,1),'linewidth',linewidth/2,'color',grey)
    xlim([0,600])

    xlabel('Time [s]');ylabel('Platform Pitch [deg]');
    taskflag = 'axes'; commonFigureTasks;

    if saveflag
        savename = ['platform_pitch_',num2str(i)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)
    end
 %--------------------------------------
    time = t;
    [fs,FFT_act] = perform_FFT(time,X_act(:,1));
    [fs,FFT_dfsm] = perform_FFT(time,X_dfsm(:,1));
    
    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [131 49 570 413];
    commonFigureProperties
    
    plot(fs,FFT_dfsm,'linewidth',linewidth,'color',red)
    plot(fs,FFT_act,'linewidth',linewidth,'color',blue)
    xlim([min(fs),max(fs)])
    xline([0.088,0.25],'linewidth',1)

    text(ww,1e-5,'Wave','FontSize',16,'Interpreter','latex')
    
    
    taskflag = 'axes'; commonFigureTasks;
    ha.XScale = 'log';
    ha.YScale = 'log';
    xlabel('Freq. [Hz]')
    ylabel(['PtfmPitch PSD',' ','[$','deg','^2$/Hz]'])
    %sgtitle('TwrBsFxt','fontsize',fonttick)
    
     if saveflag
        savename = ['ptfmpitch_',num2str(i),'_FFT'];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)
     end
end


%% plot gen speed


for i = n


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [114 46 570 413];
commonFigureProperties

results_cell_ = results_cell{i};


t = results_cell_{1};
t = t - 100;
X_cell = results_cell_{3};

X_act = X_cell{1};
X_dfsm = X_cell{2};


    plot(t,X_dfsm(:,2)*100*0.1047,'linewidth',linewidth,'color',red)
    plot(t,X_act(:,2)*100*0.1047,'linewidth',linewidth,'color',blue)
    plot(t,X_dfsm(:,2)*100*0.1047-X_act(:,2)*100*0.1047,'linewidth',linewidth/2,'color',grey)
xlim([0,600])
xlabel('Time [s]');ylabel('Generator Speed [rad/s]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['gen_speed_',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

 %--------------------------------------
    time = t;
    [fs,FFT_act] = perform_FFT(time,X_act(:,2));
    [fs,FFT_dfsm] = perform_FFT(time,X_dfsm(:,2));
    
    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [131 49 570 413];
    commonFigureProperties
    
    plot(fs,FFT_dfsm,'linewidth',linewidth,'color',red)
    plot(fs,FFT_act,'linewidth',linewidth,'color',blue)
    xlim([min(fs),max(fs)])
   xline([0.088,0.25,2.42/6.28],'linewidth',1)

     text(ww,1e-5,'Wave','FontSize',16,'Interpreter','latex')
    text(0.4,1e-7,'2P','FontSize',16,'Interpreter','latex')
    
    taskflag = 'axes'; commonFigureTasks;
    ha.XScale = 'log';
    ha.YScale = 'log';
    xlabel('Freq. [Hz]')
    ylabel(['GenSpeed PSD',' ','[$','deg','^2$/Hz]'])
    %sgtitle('TwrBsFxt','fontsize',fonttick)
    
     if saveflag
        savename = ['genspeed_',num2str(i),'_FFT'];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)
     end
end

%%
% for i = 1:6
% 
% 
% hf = figure;
% hf.Color = 'w';
% hold on;
% hf.Position = [114 69 570 413];
% commonFigureProperties
% 
% results_cell_ = results_cell{i};
% 
% t = results_cell_{1};
% t = t - 100;
% Y_cell = results_cell_{5};
% 
% Y_act = Y_cell{1};
% Y_dfsm = Y_cell{2};
% 
% 
%     plot(t,Y_dfsm(:,1),'linewidth',linewidth,'color',red)
%     plot(t,Y_act(:,1),'linewidth',linewidth,'color',blue)
%     plot(t,Y_dfsm(:,1)-Y_act(:,1),'linewidth',linewidth/2,'color',grey)
% xlim([0,600])
% xlabel('Time [s]');ylabel('TwrBsFxt [kN]');
% taskflag = 'axes'; commonFigureTasks;
% 
% if saveflag
%     savename = ['TwrBsFxt_',num2str(i)];
%     pathpdf = mfoldername(mfilename('fullpath'),fol_name);
%     filename = fullfile(pathpdf,savename);
%     str = strcat("export_fig '",filename,"' -pdf");
%     eval(str)
% end
% 
% %--------------------------------------
%     %--------------------------------------
%     time = t;
%     [fs,FFT_act] = perform_FFT(time,X_act(:,1));
%     [fs,FFT_dfsm] = perform_FFT(time,X_dfsm(:,1));
%     
%     hf = figure;
%     hf.Color = 'w';
%     hold on;
%     hf.Position = [131 49 570 413];
%     commonFigureProperties
%     
%     plot(fs,FFT_dfsm,'linewidth',linewidth,'color',red)
%     plot(fs,FFT_act,'linewidth',linewidth,'color',blue)
%     xlim([min(fs),max(fs)])
%     xline([0.022],'linewidth',1)
% 
%     text(0.025,1e-5,'Platform','FontSize',16,'Interpreter','latex')
% 
%     taskflag = 'axes'; commonFigureTasks;
%     ha.XScale = 'log';
%     ha.YScale = 'log';
%     xlabel('Freq. [Hz]')
%     ylabel(['PtfmPitch PSD',' ','[$','deg','^2$/Hz]'])
%     %sgtitle('TwrBsFxt','fontsize',fonttick)
%      if saveflag
%         savename = ['twrbsfxt_',num2str(i),'_FFT'];
%         pathpdf = mfoldername(mfilename('fullpath'),fol_name);
%         filename = fullfile(pathpdf,savename);
%         str = strcat("export_fig '",filename,"' -pdf");
%         eval(str)
%      end​
% end
% 
% 

%----------
for i = n


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [131 49 570 413];
commonFigureProperties

results_cell_ = results_cell{i};

t = results_cell_{1};
t = t - 100;
Y_cell = results_cell_{5};

Y_act = Y_cell{1};
Y_dfsm = Y_cell{2};


plot(t,Y_dfsm(:,1),'linewidth',linewidth,'color',red)
plot(t,Y_act(:,1),'linewidth',linewidth,'color',blue)
plot(t,Y_dfsm(:,1)-Y_act(:,1),'linewidth',linewidth/2,'color',grey)
xlim([0,600])
xlabel('Time [s]');ylabel('TwrBsFxt [kNm]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['TwrBsFxt_',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

%--------------------------------------
     %--------------------------------------
    time = t;
    [fs,FFT_act] = perform_FFT(time,Y_act(:,1));
    [fs,FFT_dfsm] = perform_FFT(time,Y_dfsm(:,1));
    
    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [131 49 570 413];
    commonFigureProperties
    
    plot(fs,FFT_dfsm,'linewidth',linewidth,'color',red)
    plot(fs,FFT_act,'linewidth',linewidth,'color',blue)
    xlim([min(fs),max(fs)])
    xline([0.088,0.25,2.42/6.28],'linewidth',1)

    text(ww,1e-5,'Wave','FontSize',16,'Interpreter','latex')
text(0.4,1e-5,'2P','FontSize',16,'Interpreter','latex')
    
    taskflag = 'axes'; commonFigureTasks;
    ha.XScale = 'log';
    ha.YScale = 'log';
    xlabel('Freq. [Hz]')
    ylabel(['TwrBsFxt PSD',' ','[$','deg','^2$/Hz]'])
    %sgtitle('TwrBsFxt','fontsize',fonttick)
    
     if saveflag
        savename = ['twrbsfxt_',num2str(i),'_FFT'];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)
     end
end


%----------


for i = n


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [131 49 570 413];
commonFigureProperties

results_cell_ = results_cell{i};

t = results_cell_{1};
t = t - 100;
Y_cell = results_cell_{5};

Y_act = Y_cell{1};
Y_dfsm = Y_cell{2};


plot(t,Y_dfsm(:,2),'linewidth',linewidth,'color',red)
plot(t,Y_act(:,2),'linewidth',linewidth,'color',blue)
plot(t,Y_dfsm(:,2)-Y_act(:,2),'linewidth',linewidth/2,'color',grey)
xlim([0,600])
xlabel('Time [s]');ylabel('TwrBsMyt [kNm]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['TwrBsMxt_',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

%--------------------------------------
     %--------------------------------------
    time = t;
    [fs,FFT_act] = perform_FFT(time,Y_act(:,2));
    [fs,FFT_dfsm] = perform_FFT(time,Y_dfsm(:,2));
    
    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [131 49 570 413];
    commonFigureProperties
    
    plot(fs,FFT_dfsm,'linewidth',linewidth,'color',red)
    plot(fs,FFT_act,'linewidth',linewidth,'color',blue)
    xlim([min(fs),max(fs)])
    xline([0.088,0.25,2.42/6.28],'linewidth',1)

    text(ww,1e-3,'Wave','FontSize',16,'Interpreter','latex')
     text(0.38,1e-4,'2P','FontSize',16,'Interpreter','latex')
    
    taskflag = 'axes'; commonFigureTasks;
    ha.XScale = 'log';
    ha.YScale = 'log';
    xlabel('Freq. [Hz]')
    ylabel(['TwrBsMyt PSD',' ','[$','deg','^2$/Hz]'])
    %sgtitle('TwrBsFxt','fontsize',fonttick)
    
     if saveflag
        savename = ['twrbsmyt_',num2str(i),'_FFT'];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)
     end
end

%%
% for i = 1:3
% 
% 
% hf = figure;
% hf.Color = 'w';
% hold on;
% hf.Position = [114 49 570 413];
% commonFigureProperties
% 
% results_cell_ = results_cell{i};
% 
% 
% t = results_cell_{1};
% Y_cell = results_cell_{5};
% 
% Y_act = Y_cell{1};
% Y_dfsm = Y_cell{2};
% 
% 
% plot(t,Y_dfsm(:,2),'linewidth',linewidth,'color',red)
% plot(t,Y_act(:,2),'linewidth',linewidth,'color',blue)
% plot(t,Y_dfsm(:,2)-Y_act(:,2),'linewidth',linewidth/2,'color',grey)
% 
% xlabel('Time [s]');ylabel('TwrBsMyt [kNm]');
% taskflag = 'axes'; commonFigureTasks;
% 
% if saveflag
%     savename = ['twrbsmyt',num2str(i)];
%     pathpdf = mfoldername(mfilename('fullpath'),fol_name);
%     filename = fullfile(pathpdf,savename);
%     str = strcat("export_fig '",filename,"' -pdf");
%     eval(str)
% end
% 
% 
% 
% 
% 
% end

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