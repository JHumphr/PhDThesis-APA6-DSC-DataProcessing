    function [KS0,SE0,KS,SE] = kamalFitODE(folderIndex,fitIndex,lb,ub)
% Function takes folderIndex and fitIndex (as defined by KamalFitGUI), loads FitParamFull.mat
% and outputs an updated set of KS parameters. Additionally outputs actually stdError of the
% initial parameters. May need to alter ode solver manually in this code. (as well as kamalODEfunc) 
% for troublesome datasets.

%% Load Fit Param Full
FitParamFull = load('FitParamFull.mat');
FitParamFull = FitParamFull.FitParamFull;

%% Preset Variables for Testing % Disable once operating
% folderIndex = 5;
% fitIndex = 1;
% lb = [];
% ub = [];

%% If Inputs Empty

if ~exist('lb', 'var') || isempty(lb)  
 lb = [1E4, 1E4, 1E4, 1E4, 0.5,  0.5]; end
if ~exist('ub', 'var') || isempty(ub)
 ub = [1E9, 1E5, 1E9, 1E5, 1.5,  1.5]; end
if ~exist('fitIndex', 'var') || isempty(fitIndex)
    fitIndex = 1; end

%% Load Gauss Params & Kamal Params
   gaussArray = FitParamFull{folderIndex}.gParam; %5x11 (Temp x [Param ...];
    kParam = FitParamFull{folderIndex}.kParam;
    KS0 = kParam(fitIndex,1:6); %initial parameters
    
       if(folderIndex == 4 || folderIndex == 6)
           KS0(1) = 5.5E5;
           KS0(3) = 5.5E5;
        end
        
    if fitIndex == 5 %for constant activation energy set, set bounds to the average activation energy
            
        
       lb(2) = KS0(2); 
       lb(4) = KS0(4);
       ub(2) = KS0(2); 
       ub(4) = KS0(4);
    

    end
    
    if fitIndex == 7 %change back to 7 to restore normal fucntionalitilty
        KS0(5) = 1.13;
        KS0(6) = 1.01;
        
        KS0(2) = 66500;
        KS0(4) = 52100; 
        
        
       lb(2) = KS0(2) - 0; 
       lb(4) = KS0(4) - 0;
       lb(5) = KS0(5) - 0;
       lb(6) = KS0(6) - 0;
       
       ub(2) = KS0(2) + 0; 
       ub(4) = KS0(4) + 0;
       ub(5) = KS0(5) + 0;
       ub(6) = KS0(6) + 0; 
    end
    
    stdErrorApparent = kParam(fitIndex,7);
    
    TArray = [1,130;2,140;3,150;4,160;5,170];
    tdata = (0:0.1:40)'; %minutes <--------------------!
    ydata = [];
    deleted = 0;

    % Create Gaussian Curve Again
    for TIndex = TArray(1:end,1)'
        if isnan(gaussArray(TIndex,1))==1
        TArray(TIndex-deleted,:) = [];    
        deleted = deleted+1;
        else
        [f, f1, f2] = SAGauss(tdata,gaussArray(TIndex,1:8));
        f1 = f1.*6/144; % scale to dBdt; f1 is in W/g // 144 J/g == Conv/s
        f1 = f1./0.1;
        Bfdata = cumtrapz(f1)/trapz(f1);
        Bfdata = [ydata; Bfdata];
        ydata = [ydata; f1];
        ydataPlot(:,TIndex) = f1;
        end
    end
    Tspan = TArray(:,2)';
    
KS0 = KS0./[1E5 1E4 1E5 1E4 1 1];
lb = lb./[1E5 1E4 1E5 1E4 1 1];
ub = ub./[1E5 1E4 1E5 1E4 1 1];

    [KS,SE,ymodel,Xmodel]=kamalFitODE2(tdata,ydata,Tspan,KS0,lb,ub,Bfdata);
KS = KS.*[1E5 1E4 1E5 1E4 1 1];  
    
figure
    hold on
    % plot Gaussian to fit
    plot(tdata,ydataPlot,'LineWidth',2.5)
    gcbo
%     plot(tdata,kamal6(KS0,(Bmodel,Tmodel))

% True Solution at KS0 (ODE Solver)

KS0 = KS0.*[1E5 1E4 1E5 1E4 1 1];
[ym0,xm0] = kamal6ODEfunc(KS0(1:end),Tspan');
       hold on
       for ii = 0:1:(length(xm0))/401-1
%        plot(xm0(1+401*ii:401*(ii+1),2),ym0(1+401*ii:401*(ii+1)),'k-','LineWidth',2.5)
       end
       clear ii
       resnorm0 = sum((ym0 - ydata).^2);
       SE0 = resnorm0./sqrt(length(ydata)-6);
       
%     plot(Xmodel(:,2),kamal6(KS0(1:6),[Xmodel(:,1),Xmodel(:,3)]),'o')
length(Xmodel)

        for ii = 0:1:(length(xm0))/401-1
        plot(Xmodel(1+401*ii:401*(ii+1),2),ymodel(1+401*ii:401*(ii+1)),'k:','LineWidth',2.5)
        end
%     plot(Xmodel(:,2),kamal6(KS(1:6),[Xmodel(:,1),Xmodel(:,3)]),'^')
    disp('KS0:')
    disp(KS0)
    disp('KS:')
    disp(KS)    
    fprintf('Apparent Standard Error %f \n',stdErrorApparent)
    fprintf('True Standard Error at KS0 %f \n',SE0)
    fprintf('Minimised Standard Error at KS %f \n',SE)
    
    xlabel('time (minutes)')
    ylabel('dB/dt (1/s)')
    grid minor
    grid
    legend('130 C - Data','140 C - Data','150 C - Data','160 C - Data','170 C - Data','Refit Model Parameters')
    
    %,'Original Model Parameters'
FitParamFull{folderIndex}.kParamODE(fitIndex,:) = [KS SE SE0 stdErrorApparent];    
save('FitParamFull.mat','FitParamFull');     

% load handel
% sound(y,Fs)
end

function [KS,SE,ymodel,Xmodel] = kamalFitODE2(tdata,ydata,Tdata,KS0,lb,ub,Bfdata)
% Takes output KamalParam from KamalFitGUI, and dB/dt data with tdata and Tdata to
% minimiser the error between (y'm(ym) - y'd), that is dB/dt as a function
% of conversion both calculated by the model. Less the data value.
% [t y T] should be nx3 of data values; KS0 is 1x6 of KS parameters
%% 


%%
% tspan = (0:40/60:40);
Tspan = unique(Tdata);
% length(tdata)
tspan = tdata;
yd= ydata;

func = @(Param)objFunc1(Param,tspan,Tspan,yd,Bfdata);
global loopCounter
loopCounter = 1;
% [KS,resnorm] = lsqcurvefit(func,KS0,tdata,ydata,lb,ub);
% options = optimoptions(@fmincon,'OptimalityTolerance',1e-5);


[KS,~,~,~,~,grad,hess] = fmincon(func,KS0,[],[],[],[],lb,ub,[]);%,options);
[objective,ym,xm,resnorm]=objFunc1(KS,tdata,Tspan,yd,Bfdata);
% disp(grad)
% disp(hess)
ymodel = ym;
Xmodel= xm;
KS = KS;
SE = 0;
SE = resnorm./sqrt(length(ydata)-6);
end

function [objective,ym,xm,resnorm] = objFunc1(Param,tspan,Tspan,yd,Bfdata)
% given Parameters, Tdata span and tspan solves for true ymodel value
% (based off ODE)
dBdtm = []; Bm = []; tm = []; Tm = []; ydnorm = []; Bfdnorm = [];
global odeCounter
odeCounter = 1;
global loopCounter
if mod(loopCounter,100) == 0; fprintf('Loop Counter: %d \n',loopCounter); end
loopCounter = loopCounter + 1;

k = 1;
for T = Tspan
    % Run ODE fitting
    opts = odeset('RelTol',1e-2);
    [tms,Bms]=ode45(@(t,B) kamal6ODE(t,B,Param,T),tspan,[0],opts);

    % Updated Outputs
    dBdts = kamal6ODE(tms,Bms,Param,T);
    Tms =    T.*ones(length(dBdts),1);
    dBdtm = [dBdtm; dBdts];
    Bm =    [Bm; Bms];
    tm =    [tm; tms];
    Tm =    [Tm; Tms];
    
    % Create ydata normalised weighted vector
    ydweight = yd(((k-1)*length(dBdts)+1):k*(length(dBdts)));
    Bfdweight = Bfdata(((k-1)*length(dBdts)+1):k*(length(dBdts)));
    ydweight = ydweight./max(ydweight);
    
    Bfdnorm = [Bfdnorm; Bfdweight];
    ydnorm = [ydnorm; ydweight];
    k = k+1;
    clear ydweight dBdts Bs ts    
end

ym = real(dBdtm);

xm = [Bm, tm, Tm]; 

assignin('base','ydvm',[yd ym]);
% assignin('base','ym',ym);
Bm = real(Bm);
% weightFunc = (1-Bfdnorm)+ydnorm;%(1-(Bfdnorm);%ydnorm+(1-Bfdnorm);%+ydnorm;%ydnorm;
% disp(1*(Bfdnorm < 0.95))
weightFunc = ((1-Bfdnorm+0.5).^2).*(Bfdnorm < 0.99); 

% weightFunc = 1; % turnoff weighting
objective = sum(weightFunc.*(ym - yd).^2);
resnorm = sum((ym - yd).^2);
end

function dBdt = kamal6ODE(t,B,Param,T)
% Kamal Sourour version of below isothermal (6 parameter version, solving
% for k1, k2, m, n)
global odeCounter
global loopCounter

if mod(odeCounter,1000000) == 0; fprintf('ode Counter: %d loop %d \n',odeCounter,loopCounter); end
odeCounter = odeCounter + 1;

R = 8.314;

A1 = Param(1)*1E5;
Ea1 = Param(2)*1E4;
A2 = Param(3)*1E5;
Ea2 = Param(4)*1E4;
m = Param(5);
n = Param(6);

T = T+273.15;

k1 = A1.*exp(-Ea1./R./T);
k2 = A2.*exp(-Ea2./R./T);

yNth = k1.*(1-B).^n;
yAuto = (k2.*B.^m).*(1-B).^n;
dBdt = yNth + yAuto;
dBdt = dBdt;  % rate per minute is 10* rate per increment (set to 0.1)
end %ode

function [func, SE] = objFunc2(ymodel,ydata)

func = sum((ymodel-ydata).^2);
SE = func/sqrt(length(ydata)-6);
end