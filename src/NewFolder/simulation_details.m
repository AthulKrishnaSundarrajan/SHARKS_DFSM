% Function used to load openfast simulation results, extract required
% states and controls,and filter signals as required

function sim_details =  simulation_details(simulation_files,reqd_states,reqd_controls,reqd_outputs,scale_outputs,filter_flag,filter_args,add_dx2,tmin,tmax)

    % get the number of simulations
    nsim = length(simulation_files);
    
    % length of states and controls
    nstates = length(reqd_states);
    ncontrols = length(reqd_controls);
    
    % number of inputs and outputs
    if add_dx2
        nstates = 2*nstates;       
    end

    ninputs = nstates+ncontrols;
    noutputs = length(reqd_outputs);
    nderiv = nstates;

    % initialize struct to store details
    sim_details = cell(nsim,1);

    for isim = 1:nsim
        
        % name of the file
        iname = simulation_files{isim};
        
        % extract the channels from OpenFAST simulation
        sim_detail = load_openfast_sim(iname,reqd_states,reqd_controls,reqd_outputs,scale_outputs,tmin,tmax);

        % filter if needed
        if filter_flag
            
            % filter
            sim_detail = filter_signals(sim_detail,filter_args);

        end

        % construct polynomial approximation for the states and evaluate
        % state derivatives

        sim_detail = approx_derivatives(sim_detail,add_dx2);

        % add names and length
        sim_detail.nstates = nstates;
        sim_detail.ncontrols = ncontrols;
        sim_detail.nderiv = nderiv;
        
        sim_detail.ninputs = ninputs;
        sim_detail.noutputs = noutputs;
        sim_detail.fun_name = 'FOWT';
        sim_detail.dx_act = [];

        sim_details{isim} = sim_detail;

    end
    
    % convert a cell to struct
    sim_details = cell2struct_dfsm(sim_details,nsim);

end
