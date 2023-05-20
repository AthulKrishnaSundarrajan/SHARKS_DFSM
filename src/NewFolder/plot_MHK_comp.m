clc; clear; close all;

load('MHK_simulations_comp.mat');

n = 2;

% control indices
control_ind = 1:2;
output_names{2} = 'CurrentVelX';

hf = figure;
hf.Color = 'w';

for ind = 1:2
    for i = 1:2
        
        % extract
        chan_ = chan{ind};
    
        subplot(2,1,i)
        hold on;
        plot(time,chan_(:,i),'linewidth',1.5)
    
        xlabel('Time [s]');ylabel(output_names{i})
    
    end
end


hf = figure;
hf.Color = 'w';

for ind = 1:2

    ind_ = 1;
    for i = 3:10
        
        % extract
        chan_ = chan{ind};
    
        subplot(4,2,ind_)
        hold on;
        plot(time,chan_(:,i),'linewidth',1.5)
    
        xlabel('Time [s]');ylabel(output_names{i})

        ind_ = ind_+1;
    
    end
end
legend('Constant','Step','Orientation','horizontal','Fontsize',10)

return