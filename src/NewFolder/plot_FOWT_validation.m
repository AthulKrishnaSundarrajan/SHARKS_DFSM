clc;clear; close all;

% load result
load('DFSM_FOWT_validation_Final.mat')

C = materialColors;

blue = C.blue(9,:);
yellow = C.yellow(10,:);
red = C.red(7,:);

color_map = cell(2,1);
color_map{1} = blue;
color_map{2} = red;
color_map{3} = yellow;


fol_name = 'plot_oloc_results/validation/EAB';
saveflag = false;

%% plot wind speed



for i = 1
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [131 49 570 413];
commonFigureProperties

results_cell_ = results_cell;


t = results_cell_{1};
U = results_cell_{2};

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
    hf.Position = [314 69 570 413];
    commonFigureProperties
    
    results_cell_ = results_cell;
    
    t = results_cell_{1};
    X_cell = results_cell_{3};
    
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
hf.Position = [114 46 570 413];
commonFigureProperties

results_cell_ = results_cell;


t = results_cell_{1};
X_cell = results_cell_{3};

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
hf.Position = [114 69 570 413];
commonFigureProperties

results_cell_ = results_cell;

t = results_cell_{1};
Y_cell = results_cell_{5};

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

hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [131 49 570 413];
commonFigureProperties

plot(f,P1_dfsm,'linewidth',linewidth,'color',red)
plot(f,P1_act,'linewidth',linewidth,'color',blue)

taskflag = 'axes'; commonFigureTasks;
ha.XScale = 'log';
ha.YScale = 'log';

%xlim([0,5])

end
%%
for i = 1


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [131 49 570 413];
commonFigureProperties

results_cell_ = results_cell{i};

t = results_cell_{1};
Y_cell = results_cell_{5};

Y_act = Y_cell{1};
Y_dfsm = Y_cell{2};


plot(t,Y_dfsm(:,2),'linewidth',linewidth,'color',red)
plot(t,Y_act(:,2),'linewidth',linewidth,'color',blue)
plot(t,Y_dfsm(:,2)-Y_act(:,2),'linewidth',linewidth,'color','k')

xlabel('Time [s]');ylabel('TwrBsMxt [kNm]');
taskflag = 'axes'; commonFigureTasks;

if saveflag
    savename = ['TwrBsMxt',num2str(i)];
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

hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [131 49 570 413];
commonFigureProperties

plot(f,P1_dfsm,'linewidth',linewidth,'color',red)
plot(f,P1_act,'linewidth',linewidth,'color',blue)

taskflag = 'axes'; commonFigureTasks;
ha.XScale = 'log';
ha.YScale = 'log';



end

%%
for i = 1


hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [114 49 570 413];
commonFigureProperties

results_cell_ = results_cell{i};


t = results_cell_{1};
Y_cell = results_cell_{5};

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
