function alpha = isoAvrami(Param,x)
% Kim et al Modified Avrami for Isothermal Crystallisation
% Param[1 2 3 4] = [c psi nc th];
%% Parameters
c =   Param(1);
psi = Param(2);
nc = Param(3);
th = Param(4);

%% state variables
B = x(:,1);
T = x(:,2)+273.15;
t = x(:,3);

%% Line for disabling extra model parameters
% th = 0;
%  B = 1;

%% constants
Tm0 = 220 +273.15;

%% equations
tc = c * exp(psi*Tm0./(T.*(Tm0-T)));

if (t - th) < 0
    alpha = 0;
else
    alpha = B .*(1 - exp(-((t-th)./tc).^nc)); 
end