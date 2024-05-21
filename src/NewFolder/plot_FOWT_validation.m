clc;clear; close all;

% load result
load('WEIS-P2-Q2.mat')

C = materialColors;

blue = 'k';%C.blue(9,:);
yellow = C.yellow(10,:);
red = C.red(7,:);
grey = C.grey(7,:);

color_map = cell(2,1);
color_map{1} = blue;
color_map{2} = red;
color_map{3} = yellow;


fol_name = 'plot_WEIS-P2-Q2';
saveflag = ~false;

%% plot wind speed

nr = 1;

for i = 1:nr
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [131 49 570 413];
commonFigureProperties

results_cell_ = results_cell;


t = results_cell_{1};
U = results_cell_{2};

plot(t,U(:,1),'linewidth',linewidth,'color',blue)

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


for i = 1:nr

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [314 69 570 413];
    commonFigureProperties
    
    results_cell_ = results_cell;
    
    t = results_cell_{1};
    X_cell = results_cell_{3};
    
    X_act = X_cell{1};
    X_dfsm = X_cell{2};
    
    plot(t,X_dfsm(:,1),'linewidth',linewidth,'color',red)
    plot(t,X_act(:,1),'linewidth',linewidth,'color',blue)
    %plot(t,X_dfsm(:,1)-X_act(:,1),'linewidth',linewidth/2,'color',grey)
    

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


for i = 1:nr


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [114 46 570 413];
commonFigureProperties

results_cell_ = results_cell;


t = results_cell_{1};
X_cell = results_cell_{3};

X_act = X_cell{1};
X_dfsm = X_cell{2};


    plot(t,X_dfsm(:,2)+0,'linewidth',linewidth,'color',red)
    plot(t,X_act(:,2)+0,'linewidth',linewidth,'color',blue)
    %plot(t,X_dfsm(:,2)-X_act(:,2),'linewidth',linewidth/2,'color',grey)

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
for i = 1:nr


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [114 69 570 413];
commonFigureProperties

results_cell_ = results_cell;

t = results_cell_{1};
Y_cell = results_cell_{5};

Y_act = Y_cell{1};
Y_dfsm = Y_cell{2};


    plot(t,Y_dfsm(:,1),'linewidth',linewidth,'color',red)
    plot(t,Y_act(:,1),'linewidth',linewidth,'color',blue)
    %plot(t,Y_dfsm(:,1)-Y_act(:,1),'linewidth',linewidth/2,'color',grey)

xlabel('Time [s]');ylabel('TwrBsMyt [kNm]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['TwrBsMyt_',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

%-----------------------------------------------------
fft_act = fft(Y_act(:,1));
fft_dfsm = fft(Y_dfsm(:,1));

nt = length(t);
dt = t(2)-t(1);
Fs = 1/dt;

P2_act = abs(fft_act/nt);
P2_dfsm = abs(fft_dfsm/nt);

P1_act = P2_act(1:(nt+1)/2);
P1_act(2:end-1) = 2*P1_act(2:end-1);
P1_dfsm = P2_dfsm(1:(nt+1)/2);
P1_dfsm(2:end-1) = 2*P1_dfsm(2:end-1);

f = Fs*(0:(nt)/2);

% hf = figure;
% hf.Color = 'w';
% hold on;
% hf.Position = [131 49 570 413];
% commonFigureProperties
% 
% plot(f,P1_dfsm,'linewidth',linewidth,'color',red)
% plot(f,P1_act,'linewidth',linewidth,'color',blue)
% 
% taskflag = 'axes'; commonFigureTasks;
% ha.XScale = 'log';
% ha.YScale = 'log';
%xlim([0,10^3])

% %xlim([0,5])

end
%%
for i = 1:nr


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [131 49 570 413];
commonFigureProperties

results_cell_ = results_cell;

t = results_cell_{1};
Y_cell = results_cell_{5};

Y_act = Y_cell{1};
Y_dfsm = Y_cell{2};


plot(t,Y_dfsm(:,2)+0,'linewidth',linewidth,'color',red)
plot(t,Y_act(:,2)+0,'linewidth',linewidth,'color',blue)
%plot(t,Y_dfsm(:,2)-Y_act(:,2),'linewidth',linewidth/2,'color',grey)

xlabel('Time [s]');ylabel('NcIMURAys [rad/s2]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['NcIMURAys_',num2str(i)];
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

fft_act = fft(Y_act(:,2));
fft_dfsm = fft(Y_dfsm(:,2));

nt = length(t);
dt = t(2)-t(1);
Fs = 1/dt;

P2_act = abs(fft_act/nt);
P2_dfsm = abs(fft_dfsm/nt);

P1_act = P2_act(1:(nt+1)/2);
P1_act(2:end-1) = 2*P1_act(2:end-1);
P1_dfsm = P2_dfsm(1:(nt+1)/2);
P1_dfsm(2:end-1) = 2*P1_dfsm(2:end-1);

f = Fs*(0:(nt)/2);

% hf = figure;
% hf.Color = 'w';
% hold on;
% hf.Position = [131 49 570 413];
% commonFigureProperties
% 
% plot(f,P1_dfsm,'linewidth',linewidth,'color',red)
% plot(f,P1_act,'linewidth',linewidth,'color',blue)
% 
% taskflag = 'axes'; commonFigureTasks;
% ha.XScale = 'log';
% ha.YScale = 'log';
% %xlim([0,10^3])
% 


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
