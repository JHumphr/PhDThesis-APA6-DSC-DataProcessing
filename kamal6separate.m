function [yNth, yAuto] = kamal6(Param,x)
% Kamal Sourour version of below isothermal (6 parameter version, solving
% for k1, k2, m, n)

%HOUGEN Hougen-Watson model for reaction kinetics.
%   YHAT = HOUGEN(BETA,X) gives the predicted values of the
%   reaction rate, YHAT, as a function of the vector of 
%   parameters, BETA, and the matrix of data, X.
%   BETA must have 5 elements and X must have three
%   columns.
%
%   The model form is:
%   y = (b1*x2 - x3/b5)./(1+b2*x1+b3*x2+b4*x3)
%
%   Reference:
%      [1]  Bates, Douglas, and Watts, Donald, "Nonlinear
%      Regression Analysis and Its Applications", Wiley
%      1988 p. 271-272.


%   Copyright 1993-2004 The MathWorks, Inc. 

%   B.A. Jones 1-06-95.

R = 8.314;

A1 = Param(1);
Ea1 = Param(2);
A2 = Param(3);
Ea2 = Param(4);
m = Param(5);
n = Param(6);
% disp(x)
B = x(:,1);
% t = x(:,2);
T = x(:,2)+273.15;

k1 = A1.*exp(-Ea1./R./T);
k2 = A2.*exp(-Ea2./R./T);


% disp(k1)
% disp(k2)
yNth = k1.*(1-B).^n;
yAuto = (k2).*B.^m.*(1-B).^n;
