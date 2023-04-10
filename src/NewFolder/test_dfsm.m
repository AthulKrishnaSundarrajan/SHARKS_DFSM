function [dfsm,X_cell,dx_cell] = test_dfsm(dfsm,sim_details,ind,plot_flag)

    % function to test the constructed dfsm
    ntest = length(ind);

    time_simulation = zeros(ntest,1);
    time_eval = zeros(ntest,1);

    X_cell = cell(ntest,2);
    dx_cell = cell(ntest,2);
    


    for itest = 1:ntest
        
        % extract
        controls = sim_details(itest).controls;
        states = sim_details(itest).states;
        time = sim_details(itest).time;
        state_derivatives = sim_details(itest).state_derivatives;
        outputs = sim_details(itest).outputs;
        noutputs = sim_details(itest).noutputs;

        % construct interpolating function for inputs
        u_pp = spline(time,controls');
        u_fun = @(t) ppval(u_pp,t);
        
        x0 = states(1,:)';

        % define solution options
        options = odeset('RelTol',1e-9,'AbsTol',1e-9);
        
        % run a simulation using dfsm
        tic
        [T,X] = ode45(@(t,x) ode_dfsm(t,x,u_fun,dfsm),time,x0,options);
        
        X_cell{itest,1} = states;
        X_cell{itest,2} = X;

        time_simulation(itest) = toc;
        
        % stack states and controls
        inputs = [controls,states];
        
        % evaluate dfsm
        tic
        dx_dfsm = evaluate_dfsm(inputs,dfsm,'deriv');
        time_eval(itest) = toc;
        
        % transpose
        dx_dfsm = dx_dfsm';
        dx_cell{itest,1} = state_derivatives;
        dx_cell{itest,2} = dx_dfsm;
        
        % evaluate outputs if any
        if ~isempty(outputs)
            outputs_dfsm = evaluate_dfsm(inputs,dfsm,'output');
            outputs_dfsm = outputs_dfsm';
        end
        
        %% plot

        if plot_flag
            %------------------------------------------------------------------
            % plot inputs
            hf = figure;
            hf.Color = 'w';
            sgtitle('Controls')
    
            nc = size(controls,2);
            
            for idx = 1:nc
                subplot(nc,1,idx)
                plot(time,controls(:,idx),"LineWidth",1)
                xlabel('Time [s]')
                ylabel(sim_details(itest).control_names{idx})
            end
            
            %------------------------------------------------------------------
            % plot state derivatives
            hf = figure;
            hf.Color = 'w';
            hold on;
            sgtitle('State derivatives')
    
            nx = size(states,2);
    
            for idx = 1:nx
                subplot(nx,1,idx)
                hold on;
                plot(time,dx_dfsm (:,idx),"LineWidth",1)       
                plot(time,state_derivatives(:,idx),"LineWidth",1)
    
    
                xlabel('Time [s]')
                ylabel(['d',sim_details(itest).state_names{idx}])
                legend('DFSM','OG')
            end
    
            %------------------------------------------------------------------
            % plot states
            hf = figure;
            hf.Color = 'w';
            hold on;
            sgtitle('State')
            
            for idx = 1:nx
                subplot(nx,1,idx)
                hold on;
                plot(T,X(:,idx),"LineWidth",1)
                plot(time,states(:,idx),"LineWidth",1)
                xlabel('Time [s]')
                ylabel(sim_details(itest).state_names{idx})
                legend('DFSM','OG')
            end
            
            %------------------------------------------------------------------
            % plot outputs if any
            if ~isempty(outputs)
    
                % plot
                hf = figure;
                hf.Color = 'w';
                hold on;
                sgtitle('Outputs')
                
                
                for idx = 1:noutputs
                    subplot(noutputs,1,idx)
                    hold on;
                    plot(time,outputs_dfsm(:,idx),"LineWidth",1)
                    plot(time,outputs(:,idx),"LineWidth",1)
                    xlabel('Time [s]')
                    ylabel(sim_details(itest).output_names{idx})
                    legend('DFSM','OG')
                end
    
            end
    
            %------------------------------------------------------------------
        end
    end
    dfsm.test_simulation = mean(time_simulation);
    dfsm.test_eval = mean(time_eval);

end