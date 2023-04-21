% function to construct the DFSM given inputs and options
function dfsm =  DFSM(sim_details,dfsm_options)

    % extract dfsm options
    ltype = dfsm_options.ltype;
    ntype = dfsm_options.ntype;
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

    % sample the inputs and outputs
    if isfield(dfsm_options,'use_samples')
        input_sampled = dfsm_options.use_samples.input_sampled;
        dx_sampled = dfsm_options.use_samples.dx_sampled;
        output_sampled = dfsm_options.use_samples.output_sampled;

        [~,~,~,inputs,state_dx,outputs] = sample_data(train_samples,sampling_type,nan);


    else
        tic
        [input_sampled,dx_sampled,output_sampled,inputs,state_dx,outputs] = sample_data(train_samples,sampling_type,nsamples);
        dfsm.nonlin_sample = toc;
    end


    % check the type of function
    if isempty(ltype)    % no linear function, construct a nonlinear DFSM
        
        % no linear function
        dfsm.lin = [];
        dfsm.ltype = ltype;
        dfsm.lin_sample = 0;
        dfsm.lin_construct = 0;

        
        % scale
        if scale_flag
            input_max = max(inputs,[],1);
            dx_max = max(state_dx,[],1);
        else
            input_max = ones(1,ninputs);
            dx_max = ones(1,nderiv);
        end

        input_sampled = input_sampled./input_max;

        dx_sampled = dx_sampled./dx_max;

       
        dfsm.scaler_input = input_max;
        dfsm.scaler_output = dx_max;
        
        % error
        error_ind = true(nderiv,1);
        dfsm.deriv.error_ind = error_ind;

        % check the options
        tic
        [nonlin,nonlin_construct] = construct_nonlinear_SM(input_sampled,dx_sampled,ntype,error_ind,nderiv);
        dfsm.nonlin_construct = nonlin_construct;

        % store
        dfsm.deriv.AB = [];
        dfsm.deriv.nonlin = nonlin;
        dfsm.ntype = ntype;
    
    else 

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
        dx_sampled = dx_sampled./dx_max;

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

                dx_error = dx_sampled - (input_sampled)*AB;

                if ~isempty(outputs)
                    output_error = output_sampled - input_sampled*CD;
                end
           

            case 'LPV'
                
                % number of states
                nstates = train_samples(1).nstates;
                
                % extract wind speed
                wind = inputs(:,1);
                wmin = min(wind)+1;
                wmax = max(wind)-1;
                
                % construct the LPV model
                AB = construct_LPV(inputs,state_dx,wind,wmin,wmax);
                
                % evaluate the error between original model and DFSM model
                dx_error = evaluate_LPV(AB,wmax,wmin,nstates,input_sampled,nderiv);
                
                % store values
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

        % construct corrective function
        %tic;
        if ~isempty(ntype)
            [nonlin,nonlin_construct] = construct_nonlinear_SM(input_sampled,dx_error,ntype,error_ind,nderiv);

        else
            nonlin = [];nonlin_construct = 0;
        end
        
        dfsm.nonlin_construct = nonlin_construct;
        
        % construct nonlinear corrective function for the outputs (if any)
        if ~isempty(outputs)
            op.error_ind = true(noutputs,1);
            op.nonlin = construct_nonlinear_SM(input_sampled,output_error,ntype,true(noutputs,1),noutputs);
        else
            op.nonlin = [];
        end
        
        % store
        dfsm.deriv.nonlin = nonlin;
        dfsm.op = op;
        dfsm.ntype = ntype;
    
    end
    
    % save input samples
    dfsm.input_sampled = input_sampled;
    dfsm.dx_sampled = dx_sampled;
    dfsm.output_sampled = output_sampled;

%     ind = randsample(length(train_samples),1);
%     dfsm = test_dfsm(dfsm,train_samples,ind);

    % test and plot the predictions
%     if ~isempty(test_samples)
%         ind = randsample(length(test_samples),1);
%         dfsm = test_dfsm(dfsm,test_samples,ind);
%     end


end

function [nonlin,nonlin_construct] = construct_nonlinear_SM(input,output,ntype,error_ind,noutputs)

% function to construct nonlinear surrogate model

% initialize
nonlin = cell(noutputs,1);
nonlin_construct = zeros(noutputs,1);

switch ntype

       case 'RBF'
            % mean square error goal
            mse_goal = 1e-3; 

            ind = 1;
            for i = 1:noutputs
                tic
                if error_ind(i)
                     nonlin{i} = newrb(input',(output(:,ind))',mse_goal);
                     ind = ind+1;
                end
                nonlin_construct(i) = toc;
            end



        case 'GPR'

            ind = 1;
            for i = 1:noutputs
                tic
                if error_ind(i)
                     nonlin{i} = fitrgp(input,output(:,ind),'KernelFunction','squaredexponential','OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        struct('AcquisitionFunctionName','expected-improvement-plus','UseParallel',1)); 
                     ind = ind+1;
                end
                nonlin_construct(i) = toc;
            
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


