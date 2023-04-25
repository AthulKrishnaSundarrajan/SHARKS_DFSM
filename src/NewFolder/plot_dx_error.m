clc; clear; close all;

% load results
load('TLR_dx_error_nn.mat');

% number of simulations
nsim = length(dx_cell_mf);

% initialize storage
dx_error_nl = cell(nsim,1);
dx_error_mf = cell(nsim,1);
dx_error_lin = cell(nsim,1);

X_error_nl = cell(nsim,1);
X_error_mf = cell(nsim,1);
X_error_lin = cell(nsim,1);

for i = 1:nsim
    
    % evaluate the error between derivatives
    dx_error_nl{i} = dx_cell_nl{i,1} - dx_cell_nl{i,2};
    dx_error_mf{i} = dx_cell_mf{i,1} - dx_cell_mf{i,2};
    %dx_error_lin{i} = dx_cell_lin{i,1} - dx_cell_lin{i,2};

    X_error_nl{i} = X_cell_nl{i,1} - X_cell_nl{i,2};
    X_error_mf{i} = X_cell_mf{i,1} - X_cell_mf{i,2};
    %X_error_lin{i} = X_cell_lin{i,1} - X_cell_lin{i,2};

end

% concatenate
dx_error_nl = vertcat(dx_error_nl{:});
dx_error_mf = vertcat(dx_error_mf{:});
dx_error_lin = vertcat(dx_error_lin{:});

X_error_nl = vertcat(X_error_nl{:});
X_error_mf = vertcat(X_error_mf{:});
X_error_lin = vertcat(X_error_lin{:});

% number of bins
nbins = 50;

x_lim ={[-1.0,1.0],[-1.5,1.5],[-1e-5,1e-5],[-3e-6,3e-6]};
y_lim = {[0,0.36],[0,0.3],[0,1],[0,1]};

C = materialColors;

blue = C.blue(9,:);
yellow = C.yellow(10,:);
red = C.red(7,:);

saveflag = false;


% loop through and plot
for idx = 1:4
    hf = figure;
    hf.Color = 'w';
    hf.Position = [1314 469 570 413];

    hold on;
    commonFigureProperties

    x_label = ['$\dot{\xi}_',num2str(idx),'$'];
    
    
    histogram(dx_error_nl(:,idx),'NumBins',nbins,'FaceAlpha',0.7,'FaceColor',blue,'Normalization','probability')
    
    if idx >2
        nbins_ = 1;
    else
        nbins_ = nbins;
    end

    histogram(dx_error_mf(:,idx),'NumBins',nbins,'FaceAlpha',0.7,'FaceColor',yellow,'Normalization','probability');
    xlabel(x_label)
    xlim(x_lim{idx})
    ylim(y_lim{idx})

    taskflag = 'axes'; commonFigureTasks;

    if saveflag
        savename = ['dx_error_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),'plots_final');
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

    end
    

end
close all

nbins = 70;
%------------------------------------------------------------------
x_lim = {[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]};
for idx = 1:4
    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [1314 469 570 413];
    commonFigureProperties

    x_label = ['$\xi_',num2str(idx),'$'];

    h = histogram(X_error_nl(:,idx),'NumBins',nbins,'FaceAlpha',0.7,'FaceColor',blue,'Normalization','probability')
    histogram(X_error_mf(:,idx),'NumBins',nbins,'FaceAlpha',0.7,'FaceColor',yellow,'Normalization','probability')
    xlabel(x_label)
    xlim(x_lim{idx})

    taskflag = 'axes'; commonFigureTasks;

    if saveflag
        savename = ['x_error_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),'plots_final');
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
legend({'$e_{\textrm{linfit}}$','$e_{\textrm{taylor}}$'})
fontlegend = 22; nCol = 5;lcn = 'best';
taskflag = 'legend';commonFigureTasks;


% export
saveflag = ~true;
if saveflag
    savename = 'legend_common';
    pathpdf = mfoldername(mfilename('fullpath'),'plots_final');
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

return
