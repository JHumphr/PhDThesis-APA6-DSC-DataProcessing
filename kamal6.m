function yhat = kamal6(Param,x,BI)
% Kamal Sourour version of below isothermal (6 parameter version, solving
% for k1, k2, m, n)

R = 8.314;

if isempty(BI) == 1; BI = [1,1,1,1,1,1]; end

A1 = Param(1)./BI(1);
Ea1 = Param(2)./BI(2);
A2 = Param(3)./BI(3);
Ea2 = Param(4)./BI(4);
m = Param(5)./BI(5);
n = Param(6)./BI(6);
% disp(x)
B = x(:,1);
% t = x(:,2);
T = x(:,2)+273.15;

k1 = A1.*exp(-Ea1./R./T);
k2 = A2.*exp(-Ea2./R./T);

yNth = k1.*(1-B).^n;
yAuto = (k2.*B.^m).*(1-B).^n;
yhat = yNth + yAuto;
