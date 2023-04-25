clc; clear; close all;

% set seed
rng(4357)

load('TLR_simulations_nn.mat');

nsamples = length(dx_cell_lin);

%ind = randsample(1:nsamples,1);

saveflag = false;
fol_name = 'multi-fid-validation';

%--------------------------------------------------------------------------
% plot controls

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

legend('$u_1$','$u_2$')

fontlegend = 22; nCol = 3;lcn = 'best';
taskflag = 'legend';commonFigureTasks;


if saveflag
    savename = 'control_inputs';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)

end

%-----------------------------------------
t0 = time(1); tf = time(end);

nt = length(dx_cell_nl{ind,1});
time_ = linspace(t0,tf,nt);

dx_act = dx_cell_lin{ind,1};
dx_lin = dx_cell_lin{ind,2};
dx_mf = dx_cell_mf{ind,2};
dx_nl = dx_cell_nl{ind,2};

X_act = X_cell_lin{ind,1};
X_mf = X_cell_mf{ind,2};
X_nl = X_cell_nl{ind,2};

time_x = linspace(t0,tf,length(X_mf));

blue = C.blue(9,:);
yellow = C.amber(8,:);
red = C.red(9,:);
green = C.green(9,:);
orange = C.deeporange(6,:);
black = [0,0,0];


% dx_actual vs dx_linfit
for idx = 1:4

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [1314 469 570 413];
    commonFigureProperties
    
    plot(time_,dx_lin(:,idx),'linewidth',1.5,'color',blue)
    plot(time_,dx_act(:,idx),'linewidth',1.5,'color',yellow)
    xlabel('Time [s]');ylabel(['$\dot{\xi}_',num2str(idx),'$']);

    taskflag = 'axes';commonFigureTasks

    if saveflag
        savename = ['dx_linact_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

    end



end

% plot legend
legend_entries = {'$\hat{f}_{\textrm{linfit}}$','$f_{\textrm{actual}}$'};

%----------------------
hf2 = copyfig(hf);
hf2.Position = [263.4000 417.8000 500.2000 70.4000];
h = findall(hf2,'type','line');
close(hf)

axis off; 
for i = 1:length(h)

    h(i).XData = NaN; %ignore warnings

end

% legend
legend(legend_entries);
fontlegend = 22; nCol = 3;lcn = 'best';
taskflag = 'legend';commonFigureTasks;

% export
if saveflag
    savename = 'legend_dx_linact';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

%------------------------------------------------------------------------
% plot multifidelity vs actual
for idx = 1:4

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [1314 469 570 413];
    commonFigureProperties
    
    plot(time_,dx_lin(:,idx),'linewidth',1.5,'color',blue)
    plot(time_,dx_mf(:,idx),'linewidth',1.5,'color',red)
    plot(time_,dx_act(:,idx),'linewidth',1.5,'color',yellow)
    
    xlabel('Time [s]');ylabel(['$\dot{\xi}_',num2str(idx),'$']);

    taskflag = 'axes';commonFigureTasks

    if saveflag
        savename = ['dx_MFact_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

    end



end

% plot legend
legend_entries = {'$\hat{f}_{\textrm{linfit}}$','$\hat{f}_{\textrm{multi-fid}}$','$f_{\textrm{actual}}$'};

%----------------------
hf2 = copyfig(hf);
hf2.Position = [263.4000 417.8000 500.2000 70.4000];
h = findall(hf2,'type','line');
close(hf)

axis off; 
for i = 1:length(h)

    h(i).XData = NaN; %ignore warnings

end

% legend
legend(legend_entries);
fontlegend = 22; nCol = 3;lcn = 'best';
taskflag = 'legend';commonFigureTasks;

% export
if saveflag
    savename = 'legend_dx_MFact';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

%close all;
for idx = 1:4

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [1314 469 570 413];
    commonFigureProperties

    
    plot(time_x,X_mf(:,idx),'linewidth',1.5,'color',red)
    plot(time_x,X_nl(:,idx),'linewidth',1.5,'color',black)
    plot(time_x,X_act(:,idx),'linewidth',1.5,'color',yellow)

    xlabel('Time [s]');ylabel(['$\xi_',num2str(idx),'$']);

    taskflag = 'axes';commonFigureTasks

    if saveflag
        savename = ['x_MFact_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

    end



end

% plot legend
legend_entries = {'$\hat{f}_{\textrm{multi-fid}}$','$\hat{f}_{\textrm{trad}}$','$f_{\textrm{actual}}$'};

%----------------------
hf2 = copyfig(hf);
hf2.Position = [263.4000 417.8000 500.2000 70.4000];
h = findall(hf2,'type','line');
close(hf)

axis off; 
for i = 1:length(h)

    h(i).XData = NaN; %ignore warnings

end

% legend
legend(legend_entries);
fontlegend = 22; nCol = 3;lcn = 'best';
taskflag = 'legend';commonFigureTasks;

% export
if saveflag
    savename = 'legend_dx_MF-trad-act';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end
%-------------------------------------------------------------
%plot mse error

error_trad = cell(nsamples,1);
error_MF = cell(nsamples,1);

for i = 1:nsamples

    error_trad{i} = (X_cell_nl{i,1} - X_cell_nl{i,2});
    error_MF{i} = (X_cell_mf{i,1}- X_cell_mf{i,2});


end

error_trad = vertcat(error_trad{:});
error_MF = vertcat(error_MF{:});

%----------------------------------------------------------------
nbins = 80;
%------------------------------------------------------------------
x_lim = {[-0.2,0.2],[-0.2,0.2],[-0.2,0.2],[-0.2,0.2]};
for idx = 1:4
    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [1314 469 570 413];
    commonFigureProperties

    x_label = ['$\xi_',num2str(idx),'$'];

    h = histogram(error_trad(:,idx),'NumBins',nbins,'FaceAlpha',0.7,'FaceColor',black,'Normalization','count');
    histogram(error_MF(:,idx),'NumBins',nbins,'FaceAlpha',0.7,'FaceColor',red,'Normalization','count')
    xlabel(x_label); ylabel('Count')
    xlim(x_lim{idx})

    taskflag = 'axes'; commonFigureTasks;

    if saveflag
        savename = ['x_error_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

    end
end

%--------------------------------------------------------------------------
hf2 = copyfig(hf);
hf2.Position = [263.4000 417.8000 759.2000 70.4000];
h = findall(hf2,'type','histogram');
close(hf)
axis off; 
for i = 1:length(h)

    h(i).Data = NaN; %ignore warnings

end

% legend
legend({'$e_{\textrm{trad}}$','$e_{\textrm{multi-fid}}$'})
fontlegend = 22; nCol = 5;lcn = 'best';
taskflag = 'legend';commonFigureTasks;


% export

if saveflag
    savename = 'legend_common';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end



return