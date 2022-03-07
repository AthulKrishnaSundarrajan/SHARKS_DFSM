% 
function X_ = smoothData(t,X,varargin)

if isempty(varargin)

    I = 1:size(X,2);

else

    I = varargin{1};

    if isempty(I)
        X_ = X;
    end

end


%
t_f     = 1;

% 
dt      = t(2) - t(1);

% 
nb      = floor(t_f/dt);

% 
b = ones(nb,1)/nb;

% 
Xf = filtfilt(b,1,X(:,I));

X_ = X;
X_(:,I) = Xf;


plotflag = false;

if plotflag 

figure
plot(X_,X,'.')

figure; hold on
plot(t, X,'-','linewidth',2)
plot(t, X_,'-','linewidth',2)

end

end