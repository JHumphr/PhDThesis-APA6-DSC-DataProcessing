function [ym,xm] = kamal6ODEfunc(Param,Tspan)
% given Parameters, Tdata span and tspan solves for true ymodel value
% (based off ODE)
dBdtm = []; Bm = []; tm = []; Tm = [];

Tuni = unique(Tspan);
tspan = 0:0.1:40;

for T = Tuni'
    [tms,Bms]=ode23s(@(t,B) kamal6ODE(t,B,Param,T),tspan,[0]);

    dBdts = kamal6ODE(tms,Bms,Param,T);
    Tms =    T.*ones(length(dBdts),1);
    dBdtm = [dBdtm; dBdts];
    Bm =    [Bm; Bms];
    tm =    [tm; tms];
    Tm =    [Tm; Tms];
    clear dBdts Bs ts
end

ym = real(dBdtm);
xm = [Bm, tm, Tm]; %this is currently not passed out..
end


function dBdt = kamal6ODE(t,B,Param,T)
% Kamal Sourour version of below isothermal (6 parameter version, solving
% for k1, k2, m, n)

R = 8.314;

A1 = Param(1);
Ea1 = Param(2);
A2 = Param(3);
Ea2 = Param(4);
m = Param(5);
n = Param(6);

T = T+273.15;

k1 = A1.*exp(-Ea1./R./T);
k2 = A2.*exp(-Ea2./R./T);

yNth = k1.*(1-B).^n;
yAuto = (k2.*B.^m).*(1-B).^n;
dBdt = yNth + yAuto;
dBdt = dBdt;
end %ode