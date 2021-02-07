function alpha = kimAvrami(Param,x)
% Kim et al Modified Avrami for Isothermal Crystallisation
% Param[1 2 3 4] = [aeq K theta nc];
aeq =   Param(1);
K =     Param(2);
theta = Param(3);
nc =    Param(4);

B = x(:,1);
T = x(:,2)+273.15;
t = x(:,3);

alpha = aeq.*B.*(1 - exp(-K.*(t-theta).^nc) );    

alpha(t <= theta) = 0;
