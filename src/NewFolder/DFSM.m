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
    test = train_test_split(2);

    % number of simulations
    nsim = length(sim_details);
    ninputs = sim_details(1).ninputs;
    noutputs = sim_details(1).noutputs;
    
    % extract the train test simulations
    train_split = 1:floor(train_sp*nsim);
    test_split = floor(train_sp*nsim)+1:nsim;
    
    train_samples = sim_details(train_split);
    test_samples = sim_details(train_split);

    % add number of inputs and outputs to dfsm struct
    dfsm.ninputs = ninputs;
    dfsm.noutputs = noutputs;

    % check the type of function
    if isempty(ltype)    % no linear function, construct a nonlinear DFSM
        
        % no linear function
        dfsm.lin = [];
        dfsm.ltype = ltype;
    
        % sample the inputs and outputs
        [input_sampled,output_sampled,~,~] = sample_data(train_samples,sampling_type,nsamples);

        % check the options
        switch ntype

            case 'RBF'
                % mean square error goal
                mse_goal = 0.001; 

                % train a neural network
                nonlin = newrb(input_sampled',output_sampled',mse_goal); 

            case 'GPR'

                nonlin = cell(noutputs,1);
            
                for i = 1:noutputs
                     nonlin{i} = fitrgp(input_sampled,output_sampled(:,i)); % ,'KernelFunction','squaredexponential', 'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        % struct('AcquisitionFunctionName','expected-improvement-plus'), 'Verbose',1,
                end

            case 'NN'

                trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation.

                % create a Fitting Network
                hiddenLayerSize = 15;
                net = fitnet(hiddenLayerSize,trainFcn);
                
                % setup Division of Data for Training, Validation, Testing
                net.divideParam.trainRatio = 70/100;
                net.divideParam.valRatio = 15/100;
                net.divideParam.testRatio = 15/100;
                
                % train the Network
                net = train(net,input_sampled',output_sampled');
                nonlin = net;


        end

        % store
        dfsm.nonlin = nonlin;
        dfsm.ntype = ntype;
    
    else 

         % sample the inputs and outputs
        [input_sampled,output_sampled,inputs,outputs] = sample_data(train_samples,sampling_type,lsamples);

        % construct a linear fit
        lm = linsolve(input_sampled,output_sampled);

        % store linear model
        dfsm.lin = lm;
        dfsm.ltype = ltype;

        % evaluate the error between linear model and nonlinear model
        dx_error = outputs - inputs*lm;

        % sample the input and output points
        [input_error,dx_error_sampled] = perform_KM(inputs,dx_error,nsamples);

        % corrective function
        switch ntype

            case 'RBF'
                % mean square error goal
                mse_goal = 0.001; 

                % train a neural network
                nonlin = newrb(input_error',dx_error_sampled',mse_goal);

            case 'GPR'

                nonlin = cell(noutputs,1);
            
                for i = 1:noutputs
                     nonlin{i} = fitrgp(input_error,dx_error_sampled(:,i)); % ,'KernelFunction','squaredexponential', 'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
                        % struct('AcquisitionFunctionName','expected-improvement-plus'), 'Verbose',1,
                end

            case 'NN'


        end

        dfsm.nonlin = nonlin;
        dfsm.ntype = ntype;
    
    end

    % test and plot the predictions
    test_dfsm(dfsm,test_samples,1)



end

function [input_sampled,output_sampled,inputs,outputs] = sample_data(train_samples,sampling_type,n_samples)

    % function to sample the inputs and outputs
    % two sampling stratergies available:
    %   1. equi-distant sampling - 'ED'
    %   2. k-means clustering sampling - 'KM'
    
    % extract data from structure
    [inputs,outputs] = struct2cell_dfsm(train_samples);
    
    % concatenate
    inputs = vertcat(inputs{:});
    outputs = vertcat(outputs{:});

    % sample
    switch sampling_type
 
    %     case 'ED'
    %         
    %         % get the length
    %         nt = size(inputs,1);
    % 
    %         index = 
    % 
    %         sampling_index = 
    
        case 'KM'
            % perform k means clustering
            [input_sampled,output_sampled] = perform_KM(inputs,outputs,n_samples);
    
        case 'random'
    
            % get number of inputs
            nt = size(inputs,1);
            
            % generate random
            index = randsample(nt,n_samples);
    
            % extract
            input_sampled = inputs(index,:);
            output_sampled = outputs(index,:);

    end


end

function [inputs_cell,outputs_cell] = struct2cell_dfsm(sim_details)

% get the number of samples
nsamples = length(sim_details);

inputs_cell = cell(nsamples,1);
outputs_cell = cell(nsamples,1);

% loop through and extract cell
for isample = 1:nsamples

    controls = sim_details(isample).controls;
    states = sim_details(isample).states;
    state_derivatives = sim_details(isample).state_derivatives;

    % combine
    inputs = [states,controls];
    
    % store
    inputs_cell{isample} = inputs;
    outputs_cell{isample} = state_derivatives;

end
end

function test_dfsm(dfsm,sim_details,ind)

    % function to test the constructed dfsm
    ntest = length(ind);


    for itest = 1:ntest
        
        % extract
        controls = sim_details(itest).controls;
        states = sim_details(itest).states;
        time = sim_details(itest).time;
        state_derivatives = sim_details(itest).state_derivatives;
        noutputs = sim_details(itest).noutputs;

        u_fun = @(t) interp1(time,controls,t,'spline');
        x0 = states(1,:)';

        % define solution options
        options = odeset('RelTol',1e-5,'AbsTol',1e-5);

        [T,X] = ode45(@(t,x) ode_dfsm(t,x,u_fun,dfsm),[time(1),time(end)],x0,options);

        inputs = [states,controls];

        outputs = evaluate_dfsm(dfsm,inputs);
        outputs = outputs';


        % plot
        hf = figure;
        hf.Color = 'w';
        hold on;
        sgtitle('State derivatives')
        
        idx = 1;
        subplot(2,1,idx)
        hold on;
        plot(time,outputs(:,idx),"LineWidth",1)
        plot(time,state_derivatives(:,idx),"LineWidth",1)
        legend('DFSM','OF')

        idx = idx +1;
        subplot(2,1,idx)
        hold on;
        plot(time,outputs(:,idx),"LineWidth",1)
        plot(time,state_derivatives(:,idx),"LineWidth",1)
        legend('DFSM','OF')
        
        % plot
        hf = figure;
        hf.Color = 'w';
        hold on;
        sgtitle('State')
        
        idx = 1;
        subplot(2,1,idx)
        hold on;
        plot(T,X(:,idx),"LineWidth",1)
        plot(time,states(:,idx),"LineWidth",1)
        legend('DFSM','OF')

        idx = idx +1;
        subplot(2,1,idx)
        hold on;
        plot(T,X(:,idx),"LineWidth",1)
        plot(time,states(:,idx),"LineWidth",1)
        legend('DFSM','OF')



    end
end

function dx = ode_dfsm(t,x,u_fun,dfsm)
    
    % controls
    u = u_fun(t);
    u = u(:);

    % combine
    inputs = [x',u'];
    
    dx = evaluate_dfsm(dfsm,inputs);

end

function[input_sampled,output_sampled] = perform_KM(inputs,outputs,n_samples)
    % get clusters
    [idx,input_sampled] = kmeans(inputs,n_samples,'MaxIter',1000);

    output_sampled = zeros(n_samples,size(outputs,2));

    for icluster = 1:n_samples
        
        % find the index corresponding to given cluster
        ind_cluster = (idx == icluster);
        
        % find the state derivative values for those clusters
        output_cluster = outputs(ind_cluster,:);
        
        % the state derivative value at the cluster centroid is the
        % mean
        output_sampled(icluster,:) = mean(output_cluster,1);


    end

end

function dx = evaluate_dfsm(dfsm,inputs)

    % extract dfsm 
    ltype = dfsm.ltype;ntype = dfsm.ntype;
    lin = dfsm.lin; nonlin = dfsm.nonlin;

    % number of points
    nt = size(inputs,1);
    ninputs = dfsm.ninputs;
    noutputs = dfsm.noutputs;

    if isempty(ltype)

        dx_lin = zeros(nt,noutputs);
        
    else
        % evaluate linear part of the dfsm
        dx_lin = inputs*lin;
    end


    switch ntype

        case 'RBF'

            dx_nonlin = nonlin(inputs');

        case 'GPR'

            dx_nonlin = zeros(nt,noutputs);

            for i = 1:noutputs

                gpri = nonlin{i};

                dx_nonlin(:,i) = predict(gpri,inputs);

            end

            dx_nonlin = dx_nonlin';

        case 'NN'

            dx_nonlin = nonlin(inputs');

    end

    dx = dx_lin' + dx_nonlin;


end