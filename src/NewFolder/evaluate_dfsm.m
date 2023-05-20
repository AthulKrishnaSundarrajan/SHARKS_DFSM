function dx = evaluate_dfsm(inputs,dfsm,fun_type)

    % number of points
    nt = size(inputs,1);

    % extract dfsm options
    ltype = dfsm.ltype;
    ntype = dfsm.ntype;
    error_ind = dfsm.deriv.error_ind;
    scaler_input = dfsm.scaler_input;


    % based on the function type, extract the dfsm
    switch fun_type

        case 'deriv'

            lin = dfsm.deriv.AB;
            nonlin = dfsm.deriv.nonlin;
            noutputs = dfsm.nderiv;

        case 'output'

            lin = dfsm.op.CD;
            nonlin = dfsm.op.nonlin;
            noutputs = dfsm.noutputs;

    end
   
    % scale inputs
    inputs = inputs./scaler_input;

    % initialize 
    if isempty(ltype)

        % if no linear function then set dx_lin as zero
        dx_lin = zeros(nt,noutputs);
      
    else
        % else evaluate dx_lin
        switch ltype

            case 'LTI'

                % evaluate linear part of the dfsm
                dx_lin = inputs*lin;

            case 'LPV'
                % evaluate lpv 
                dx_lin = evaluate_LPV(lin,dfsm.lpv.wmax,dfsm.lpv.wmin,inputs,noutputs);

        end
    end

    % initialize
    dx_nonlin = zeros(noutputs,nt);

    switch ntype

        case 'RBF'
 
           for i = 1:noutputs

               if error_ind(i)
                   rbfi = nonlin{i};

                   dx_nonlin(i,:) = rbfi(inputs');
               end

           end

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

            for i = 1:noutputs

               if error_ind(i)
                   nni = nonlin{i};

                   dx_nonlin(i,:) = nni(inputs');
               end

           end



    end
    
    % add 
    dx = dx_lin' + dx_nonlin;

end