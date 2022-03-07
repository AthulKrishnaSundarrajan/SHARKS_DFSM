function data = approximateStateDerivatives(data)

% TODO: save original data

plotflag = false;
filterflag = true;
d_PtfmPitchflag = trur;

% number of input DLCs
nDLCs = length(data);

% go through each DLC
for iCase = 1:nDLCs

    % extract
    t = data(iCase).time;
    x = data(iCase).states;
    u = data(iCase).inputs;


    % filter data
    if filterflag
        u_ = smoothData(t,u,[1]);
        x_ = smoothData2(t,x);
    else
        u_ = u;
        x_ = x;
    end

    if plotflag

        u_org = u;
        x_org = x;

%         figure; hold on
%         umax =  max(abs([u_;u]),[],1);
%         plot(t,u_./umax); plot(t,u./umax); 
%     
%         figure; hold on
%         xmax =  max(abs([x_;x]),[],1);
%         plot(t,x_./xmax); plot(t,x./xmax); 


        C = materialColors;
                ind = 1;

        hf = figure; hold on
        hf.Color  = 'w';
        hf.Position = [1000 918 720 420];

        t_more = linspace(min(t),max(t),1e6);
        pp_x = spline(t,x_');
        xx = ppval(pp_x,t_more);

        xmax =  max(abs([x_;x]),[],1);
        plot(t_more,xx(ind,:),'linewidth',1.5,'Color',C.blue(10,:))
        plot(t,x_(:,ind),'.','Color',C.blue(5,:),'markersize',4);


        xlim([t(1) t(1)+5])

        ha = gca;
        ha.FontSize = 16;
        ha.LineWidth = 1;
        xlabel('Time (s)')
        ylabel('Platform Pitch [deg]')
        legend('Polynomial Approx.','OpenFast Data')

    end

    % update
    u = u_;
    x = x_;

    % construct polynomial interpolates of the states
    pp_x = spline(t,x');

    % take derivatives of polynomial interpolates
    pp_Dx = fnder(pp_x,1);
    pp_D2x = fnder(pp_x,2);

    % calculate approximate derivatives
    Dx = ppval(pp_Dx,t);
    D2x = ppval(pp_D2x,t);

    if plotflag

        C = materialColors;

        pp_x_org = spline(t,x_org');
        pp_Dx_org = fnder(pp_x_org,1);
        Dx_org = ppval(pp_Dx_org,t);

        hf = figure; hold on
        hf.Color  = 'w';
        hf.Position = [1000 918 720 420];
       

        Dxmax =  max(abs(Dx),[],2);

        ind = 1;

%         t_more = linspace(min(t),max(t),1e6);
%         Dx_ = ppval(pp_Dx,t_more);
%         Dx_org_ = ppval(pp_Dx_org,t_more);
%         plot(t_more,Dx_(ind,:)./Dxmax(ind),'linewidth',1,'Color',C.blue(10,:))
% %         plot(t_more,Dx_org_(ind,:)./Dxmax(ind),'linewidth',0.5,'Color',C.red(10,:))
% 
% 
%         plot(t,Dx(ind,:)./Dxmax(ind),'.','Color',C.blue(5,:)); 
% %         plot(t,Dx_org(ind,:)./Dxmax(ind),'Color',C.red(5,:)); 

        t_more = linspace(min(t),max(t),1e6);
        Dx_ = ppval(pp_Dx,t_more);
        plot(t_more,Dx_(ind,:),'linewidth',1.5,'Color',C.blue(10,:))
%         plot(t,Dx(ind,:),'.','Color',C.blue(5,:)); 

        plot([100 105],[0 0],'--k')

        xlim([t(1) t(1)+5])

        ha = gca;
        ha.FontSize = 16;
        ha.LineWidth = 1;
        xlabel('Time (s)')
        ylabel('Derivative of Platform Pitch [deg/s]')
%         legend()

    end

    if d_PtfmPitchflag

        % update states
        data(iCase).states = [x,Dx(1,:)'];
        data(iCase).state_names{end+1} = char(strcat("d_",data(iCase).state_names{1}));
    
        % assign derivatives
        data(iCase).state_derivatives = [D2x(1,:)',Dx(2,:)'];
        data(iCase).state_derivative_names = {'d2_PtfmPitch','d_GenSpeed'};

    else

        % assign derivatives
        data(iCase).state_derivatives = Dx';
        data(iCase).state_derivative_names = {'d1_PtfmPitch','d_GenSpeed'};

    end

end

end

% integrate
%     pp_Ix = fnint(pp_x,ppval(pp_x,t(1)));