% Function used to load openfast simulation results, extract required
% states and controls,and filter signals as required

function sim_details =  simulation_details(simulation_files,reqd_states,reqd_controls,filter_flag,filter_args,tmin,tmax)

    % get the number of simulations
    nsim = length(simulation_files);
    
    % length of states and controls
    nstates = length(reqd_states);
    ncontrols = length(reqd_controls);
    
    % number of inputs and outputs
    ninputs = nstates+ncontrols;
    noutputs = nstates;

    % initialize struct to store details
    sim_details = cell(nsim,1);

    for isim = 1:nsim
        
        % name of the file
        iname = simulation_files{isim};
        
        % extract the channels from OpenFAST simulation
        sim_detail = load_openfast_sim(iname,reqd_states,reqd_controls,tmin,tmax);

        % filter if needed
        if filter_flag
            
            % filter
            sim_detail = filter_signals(sim_detail,filter_args);

        end

        % construct polynomial approximation for the states and evaluate
        % state derivatives

        sim_detail = approx_derivatives(sim_detail);

        % add names and length
        sim_detail.nstates = nstates;
        sim_detail.ncontrols = ncontrols;
        sim_detail.ninputs = ninputs;
        sim_detail.noutputs = noutputs;

        sim_details{isim} = sim_detail;

    end
    
    % convert a cell to struct
    sim_details = cell2struct(sim_details,nsim);

end

function sim_details = cell2struct(sim_details_cell,nsim)

% get length
sim_details(nsim) = struct();

% loop through
for i = 1:nsim
    
    % initialize
    sim_details(i).time = [];sim_details(i).states = [];sim_details(i).controls = [];
    sim_details(i).state_names = [];sim_details(i).control_names = [];sim_details(i).state_derivatives = [];
    sim_details(i).nstates = [];sim_details(i).ncontrols = [];sim_details(i).ninputs = [];sim_details(i).noutputs = [];
    
    % assign
    sim_details(i) = sim_details_cell{i};

end

end