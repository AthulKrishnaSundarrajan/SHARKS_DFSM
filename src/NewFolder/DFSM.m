% function to construct the DFSM given inputs and options
function dfsm =  DFSM(sim_details,dfsm_options)

    % extract dfsm options
    ltype = dfsm_options.ltype;
    ntype = dfsm_options.ntype;
    lsamples = dfsm_options.lsamples;
    nsamples = dfsm_options.nsamples;
    sampling_type = dfsm_options.sampling_type;
    train_test_split = dfsm_options.train_test_split;

    % get the train/test splits
    train_sp = train_test_split(1);
    test_sp = train_test_split(2);

    % number of simulations
    nsim = length(sim_details);
    ninputs = sim_details(1).ninputs;
    nderiv = sim_details(1).nderiv;
    noutputs = sim_details(1).noutputs;
    
    % extract the train test simulations
    train_split = 1:floor(train_sp*nsim);
    test_split = floor(train_sp*nsim)+1:nsim;
    
    train_samples = sim_details(train_split);
    test_samples = sim_details(test_split);

    % add number of inputs and outputs to dfsm struct
    dfsm.ninputs = ninputs;
    dfsm.noutputs = noutputs;
    dfsm.nderiv = nderiv;

    scale_flag = dfsm_options.scale_flag;

    % check the type of function
    if isempty(ltype)    % no linear function, construct a nonlinear DFSM
        
        % no linear function
        dfsm.lin = [];
        dfsm.ltype = ltype;
        dfsm.lin_sample = 0;
        dfsm.lin_construct = 0;

        % sample the inputs and outputs
        tic
        [input_sampled,state_dx_sampled,output_sampled,inputs,state_dx,outputs] = sample_data(train_samples,sampling_type,nsamples);
        dfsm.nonlin_sample = toc;

        % scale
        if scale_flag
            input_max = max(inputs,[],1);
            output_max = max(state_dx,[],1);
        else
            input_max = ones(1,ninputs);
            output_max = ones(1,noutputs);
        end

        input_sampled = input_sampled./input_max;

        state_dx_sampled = state_dx_sampled./output_max;

       
        dfsm.scaler_input = input_max;
        dfsm.scaler_output = output_max;
        
        % error
        error_ind = true(noutputs,1);
        dfsm.error_ind = error_ind;

        % check the options
        tic
        nonlin = construct_nonlinear_SM(input_sampled,state_dx_sampled,ntype,error_ind,noutputs);
        dfsm.nonlin_construct = toc;

        % store
        dfsm.nonlin = nonlin;
        dfsm.ntype = ntype;
    
    else 

        % sample the inputs and outputs
        tic
        [input_sampled,state_dx_sampled,output_sampled,inputs,state_dx,outputs] = sample_data(train_samples,sampling_type,lsamples);
        dfsm.lin_sample = toc;

        % scale
        if scale_flag
            input_max = max(inputs,[],1);
            dx_max = max(state_dx,[],1);
            %output_max = max(outputs,[],1);
        else
            input_max = ones(1,ninputs);
            dx_max = ones(1,nderiv);
        end

        inputs = inputs./input_max;
        input_sampled = input_sampled./input_max;

        state_dx = state_dx./dx_max;
        state_dx_sampled = state_dx_sampled./dx_max;

        % construct a linear fit
        tic
        switch ltype
            
           
            case 'LTI'
            
                AB = linsolve(inputs,state_dx);
                
                if ~isempty(outputs)
                    CD = linsolve(inputs,outputs);
                else
                    CD = [];
                end

                dfsm.deriv.AB = AB;
                op.CD = CD;

                dx_error = state_dx_sampled - (input_sampled)*AB;

                if ~isempty(outputs)
                    output_error = output_sampled - input_sampled*CD;

                end
           

            case 'LPV'

                nstates = train_samples(1).nstates;
                
                % extract wind speed
                wind = inputs(:,nstates+1);
                wmin = min(wind)+1;
                wmax = max(wind)-1;

                AB = construct_LPV(inputs,state_dx,wind,wmin,wmax);
                

                dx_error = evaluate_LPV(AB,wmax,wmin,nstates,input_sampled,nderiv);

                dfsm.deriv.lin = AB;
                dfsm.lpv.wmin = wmin;
                dfsm.lpv.wmax = wmax;
                dfsm.lpv.nstates = nstates;

        end
        dfsm.lin_construct = toc;

        % store linear model
        
        dfsm.ltype = ltype;
        dfsm.scaler_input = input_max;
        dfsm.scaler_output = dx_max;

        % evaluate the error between linear model and nonlinear model

        error_mean = mean(dx_error,1);

        % check mean value of error. if error <1e-5, then do not construct
        % corrective function for that derivative
        error_ind = abs(error_mean) > 1e-5; 
        dx_error = dx_error(:,error_ind);
        dfsm.deriv.error_ind = error_ind;

        % corrective function
        tic;
        nonlin = construct_nonlinear_SM(input_sampled,dx_error,ntype,error_ind,nderiv);
        dfsm.nonlin_construct = toc;

        if ~isempty(outputs)
            op.error_ind = true(noutputs,1);
            op.nonlin = construct_nonlinear_SM(input_sampled,output_error,ntype,true(noutputs,1),noutputs);
        else
            op.nonlin = [];
        end

        dfsm.deriv.nonlin = nonlin;
        dfsm.op = op;
        dfsm.ntype = ntype;
    
    end


    ind = randsample(length(train_samples),2);
    dfsm = test_dfsm(dfsm,train_samples,ind);

    % test and plot the predictions
    dfsm = test_dfsm(dfsm,test_samples,1:length(test_samples));



end

function nonlin = construct_nonlinear_SM(input,output,ntype,error_ind,noutputs)

% function to construct nonlinear surrogate model

switch ntype

        case 'RBF'
            % mean square error goal
            mse_goal = 0.001; 

            % train a neural network
            nonlin = newrb(input',output',mse_goal);

        case 'GPR'

            nonlin = cell(noutputs,1);
            ind = 1;
            for i = 1:noutputs
                if error_ind(i)
                     nonlin{i} = fitrgp(input,output(:,ind),'KernelFunction','squaredexponential','OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        struct('AcquisitionFunctionName','expected-improvement-plus','UseParallel',1)); 
                     ind = ind+1;
                end
            
            end

        case 'NN'

            trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation.

            % create a Fitting Network
            hiddenLayerSize = 60;
            net = fitnet(hiddenLayerSize,trainFcn);
            
            % setup Division of Data for Training, Validation, Testing
            net.divideParam.trainRatio = 70/100;
            net.divideParam.valRatio = 15/100;
            net.divideParam.testRatio = 15/100;
            net.trainParam.max_fail = 5000;
            
            % train the Network
            net = train(net,input',output');
            nonlin = net;


end


end

function [input_sampled,state_dx_sampled,output_sampled,inputs,state_dx,outputs] = sample_data(train_samples,sampling_type,n_samples)

    % function to sample the inputs and outputs
    % two sampling stratergies available:
    %   1. equi-distant sampling - 'ED'
    %   2. k-means clustering sampling - 'KM'
    
    % extract data from structure
    if isstruct(train_samples)
        [inputs,state_dx,outputs] = struct2cell_dfsm(train_samples);

        % concatenate
        inputs = vertcat(inputs{:});
        state_dx = vertcat(state_dx{:});
        outputs = vertcat(outputs{:});

    elseif iscell(train_samples)

        inputs = train_samples{1};
        state_dx = train_samples{2};
        outputs = train_samples{3};


    end

    if isnan(n_samples)
        input_sampled = inputs;
        state_dx_sampled = state_dx;
        output_sampled = outputs;

        return
    end


    % sample
    switch sampling_type
 
    
        case 'KM'
            % perform k means clustering
            [input_sampled,state_dx_sampled,output_sampled] = perform_KM(inputs,state_dx,outputs,n_samples);
    
        case 'random'
    
            % get number of inputs
            nt = size(inputs,1);
            
            % generate random
            index = randsample(nt,n_samples);
    
            % extract
            input_sampled = inputs(index,:);
            state_dx_sampled = state_dx(index,:);
            output_sampled = outputs(index,:);

    end


end

function [inputs_cell,state_dx_cell,outputs_cell] = struct2cell_dfsm(sim_details)

% get the number of samples
nsamples = length(sim_details);

inputs_cell = cell(nsamples,1);
state_dx_cell = cell(nsamples,1);
outputs_cell = cell(nsamples,1);

% loop through and extract cell
for isample = 1:nsamples

    controls = sim_details(isample).controls;
    states = sim_details(isample).states;
    state_derivatives = sim_details(isample).state_derivatives;
    outputs = sim_details(isample).outputs;

    % combine
    inputs = [states,controls];
    
    % store
    inputs_cell{isample} = inputs;
    state_dx_cell{isample} = state_derivatives;
    outputs_cell{isample} = outputs;

end
end

function dfsm = test_dfsm(dfsm,sim_details,ind)

    % function to test the constructed dfsm
    ntest = length(ind);

    time_simulation = zeros(ntest,1);
    time_eval = zeros(ntest,1);
    


    for itest = 1:ntest
        
        % extract
        controls = sim_details(itest).controls;
        states = sim_details(itest).states;
        time = sim_details(itest).time;
        state_derivatives = sim_details(itest).state_derivatives;
        outputs = sim_details(itest).outputs;
        noutputs = sim_details(itest).noutputs;
        nderiv = sim_details(itest).nderiv;

        u_pp =spline(time,controls');
        u_fun = @(t) ppval(u_pp,t);
        
        x0 = states(1,:)';

        % define solution options
        options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        
        tic
        [T,X] = ode45(@(t,x) ode_dfsm(t,x,u_fun,dfsm),[time(1),time(end)],x0,options);
        time_simulation(itest) = toc;

        inputs = [states,controls];
        
        tic
        state_dx = evaluate_dfsm(inputs,dfsm.ltype,dfsm.ntype,dfsm.deriv.AB,dfsm.deriv.nonlin,dfsm.nderiv,dfsm.deriv.error_ind,dfsm.scaler_input,dfsm.scaler_output);
        time_eval(itest) = toc;

        state_dx = state_dx';
        
        if ~isempty(outputs)
            outputs_dfsm = evaluate_dfsm(inputs,dfsm.ltype,dfsm.ntype,dfsm.op.CD,dfsm.op.nonlin,dfsm.noutputs,dfsm.op.error_ind,dfsm.scaler_input,dfsm.scaler_output);
            outputs_dfsm = outputs_dfsm';
        end

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

        % plot
        hf = figure;
        hf.Color = 'w';
        hold on;
        sgtitle('State derivatives')

        nx = size(states,2);

        for idx = 1:nx
            subplot(nx,1,idx)
            hold on;
            plot(time,state_dx (:,idx),"LineWidth",1)       
            plot(time,state_derivatives(:,idx),"LineWidth",1)


            xlabel('Time [s]')
            ylabel(['d',sim_details(itest).state_names{idx}])
            legend('DFSM','OG')
        end

        
        
        % plot
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


    end

    dfsm.test_simulation = mean(time_simulation);
    dfsm.test_eval = mean(time_eval);

end

function dx = ode_dfsm(t,x,u_fun,dfsm)
    
    % controls
    u = u_fun(t);
    u = u(:);

    % combine
    inputs = [x',u'];
    
    dx = evaluate_dfsm(inputs,dfsm.ltype,dfsm.ntype,dfsm.deriv.AB,dfsm.deriv.nonlin,dfsm.nderiv,dfsm.deriv.error_ind,dfsm.scaler_input,dfsm.scaler_output);

end

function[input_sampled,state_dx_sampled,output_sampled] = perform_KM(inputs,state_dx,outputs,n_samples)

    % get clusters
    [idx,input_sampled] = kmeans(inputs,n_samples,'MaxIter',1000);

    state_dx_sampled = zeros(n_samples,size(state_dx,2));
    output_sampled = zeros(n_samples,size(outputs,2));

    for icluster = 1:n_samples
        
        % find the index corresponding to given cluster
        ind_cluster = (idx == icluster);
        
        % find the state derivative values for those clusters
        state_dx_cluster = state_dx(ind_cluster,:);

        if ~isempty(outputs)
         output_cluster = outputs(ind_cluster,:);

         output_sampled(icluster,:) = mean(output_cluster,1);
        end
        
        % the state derivative value at the cluster centroid is the
        % mean
        state_dx_sampled(icluster,:) = mean(state_dx_cluster,1);

    end

end

function dx = evaluate_dfsm(inputs,ltype,ntype,lin,nonlin,noutputs,error_ind,scaler_input,scaler_output)

    % extract dfsm 
%     ltype = dfsm.ltype;ntype = dfsm.ntype;
%     nonlin = dfsm.deriv.nonlin;

    % number of points
    nt = size(inputs,1);

%     noutputs = dfsm.noutputs;
%     nderiv = dfsm.nderiv;
%     error_ind = dfsm.error_ind;
%     scaler_input = dfsm.scaler_input;
%     scaler_output = dfsm.scaler_output;

    inputs = inputs./scaler_input;

    % initialize 
    if isempty(ltype)

        dx_lin = zeros(nt,nderiv);
        
    else

        switch ltype

            case 'LTI'

                % evaluate linear part of the dfsm
                dx_lin = inputs*lin;

            case 'LPV'

                dx_lin = evaluate_LPV(dfsm.lpv.lin,dfsm.lpv.wmax,dfsm.lpv.wmin,dfsm.lpv.nstates,inputs,noutputs);

        end
    end

    % initialize
    dx_nonlin = zeros(noutputs,nt);

    switch ntype


        case 'RBF'
 
            dx_nonlin(error_ind,:) = nonlin(inputs');

        case 'GPR'

            dx_nonlin = zeros(nt,noutputs);

            for i = 1:noutputs

                if error_ind(i)
                    gpri = nonlin{i};
    
                    dx_nonlin(:,i) = predict(gpri,inputs);
                end

            end

            dx_nonlin = dx_nonlin';

        case 'NN'

            dx_nonlin(error_ind,:) = nonlin(inputs');

    end

    dx = dx_lin' + dx_nonlin;
    %dx = dx.*scaler_output';


end

function lin = construct_LPV(inputs,outputs,wind,wmin,wmax)

    % function to construct LPV model

    % get the maximum and minimum values of the wind speed
    maxW = max(wind);
    minW = min(wind);
    
    % get the size of inputs and outputs
    ninputs = size(inputs,2);
    noutputs = size(outputs,2);
    
    % create a grid of wind speed values
    W = linspace(wmin,wmax,50);
    
    % initialize offset
    offset = 3;
    
    % initialize storage array
    IW_ = zeros(length(W),1);

    AB_interp = zeros(50,ninputs,noutputs);

for k = 1:length(W)

    m = W(k);

    if m+offset > maxW
        offset = maxW - m;
    end

    if m-offset < minW
        offset = m-minW;
    end

    % find all values within small range
    Iw = (m+offset >= wind) & (m-offset <= wind);
    IW_(k) = sum(Iw);

    % extract
    input_ = inputs(Iw,:);
    output_ = outputs(Iw,:);

    AB_interp(k,:,:) = linsolve(input_,output_);

end

lin_pp = interp1(W,AB_interp,'nearest','pp');
lin = @(w) ppval(lin_pp,w);


end

function dx = evaluate_LPV(lin,wmax,wmin,nstates,inputs,noutputs)

% function to evaluate the LPV model for a given set of inputs
% works for both direct evaluation and through ode

% get the size of inputs
nt = size(inputs,1);

% initialize derivative
dx = zeros(nt,noutputs);

% go through the inputs and evaluate lpv model
for i = 1:nt

    % extract wind
    w = inputs(i,nstates+1);
    
    % check value to see if the wind speed is between limits
    if w > wmax
        w = wmax;
    elseif w<wmin
        w = wmin;
    end
    
    % evaluate derivative function
    dx(i,:) = inputs(i,:)*lin(w);
end

end
