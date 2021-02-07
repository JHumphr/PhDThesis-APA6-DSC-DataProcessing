%% Takes in the Results "Vector" from LoadSeparation.m to Compute Kamal-Sourour Model
clear all
close all
clc
%% Manual Exclusion of Data:
Exc = true(10,5,3);

% Temperature index:
% 130 - 1; 140 -2; 150 -3; 160-4; 170-5; 180 -6;   Index = (T-120)/10;
TRNG = [130 140 150 160 170 180];
%%%%%%%[ 1   2   3   4   5   6 ];
%
% Exclusion List is done manually, this is for clearly irrelvant data
% (negative, too small, etc.) This was done based on the polymerisation
% curves outputted by Load Separation Exc(Folder,Temp,RunNumber)
% Folder 1:
% T = 140: Exclude 1&2
Exc(1,2,1) = false;
Exc(1,2,2) = false;
% T = 150: Exclude 1&3:
Exc(1,3,1) = false;
Exc(1,3,3) = false; %(this one is marginal)
% T = 160: Exclude 1:
Exc(1,4,1) = false;
% T = 170: Exclude 3:
Exc(1,5,1) = false;
Exc(1,5,3) = false;

% Folder 2:
% T=140: Exclude All
Exc(2,2,1) = false;
Exc(2,2,2) = false;
% Exc(2,2,3) = false;
%T=150: Three very different results... auto correlation?
Exc(2,3,3) = false;
%T=160: Similar as 150.. 
Exc(2,4,2) = false; 
%T=170:
Exc(2,5,2) = false;

%Folder 3:
%T=140:
Exc(3,2,2) = false;
%T=150:
Exc(3,3,1) = false;
% Exc(3,3,3) = false; %this one is marginal
%T=160:
Exc(3,4,1) = false; %does this work?
%T=170:

%Folder 4:
%T=140:
% Exc(4,2,1) =false;
Exc(4,2,2) =false;
%T=150:
Exc(4,3,3) =false; %does this work?
%T=160:
Exc(4,4,1) =false;
%T=170:
Exc(4,5,2) = false;

%Folder 5:
%T=130:
Exc(5,1,2) =false;
%T=140:
Exc(5,2,2) =false;
%T=150:
Exc(5,3,2:3) =false; %these are okay but wobbly
%T=160:
% perfect dataset..
%T=170:
% good dataset..

%Folder 6:
%T=130:
Exc(6,1,1:2) =false; %maybe these lower T samples have  aproblem with the polymerisation separationg step
%T=140:
Exc(6,2,1:2) =false;
%T=150:
Exc(6,3,3) =false;
%T=160:
Exc(6,4,1:2) =false;
%T=170:
Exc(6,5,2) =false;

%Folder 7:
%T=130:
Exc(7,1,1:2) =false; %might be possible ot keep 2
%T=140:
Exc(7,2,2) = false;
%T=150:
% Exc(7,3,3) = false; %probably okay
%T=160:
Exc(7,4,2:3) = false;
%T=170:
%okay

%Folder 8:
%T=130:
Exc(8,1,2) = false;
%T=140:
Exc(8,2,1) = false;
%T=150:
%T=160:
Exc(7,4,1) = false;
%T=170:


%Folder 9:
%T=130:
Exc(9,1,2) = false;
%T=140:
%T=150:
Exc(9,3,3) = false;
%T=160:

%T=170:
Ex(9,5,2) = false;

%Folder 10: %these became completely excluded after some change in the
%previous code
%T=130:
Exc(10,1,1:2) = false; %its difficult to say which is better (1 or 2)
%T=140:
Exc(10,2,1:2) = false;
%T=150:
Exc(10,3,1) = false;
Exc(10,3,2) = false;
%T=160:
Exc(10,4,2) = false;
%T=170:
%ok

% Exc = true(10,5,3);

%% 
load('Vector.mat')

for f = 1:1:length(Vector) %for each folder
    V = Vector{f}; % [Ti Rn k Hmelt] where k is a vector of the Gaussian parameters 1-8
    Tspan = unique(Vector{f}(:,1));
    DataVar = []; %[t T B dBdt]
    
    for i = 1:1:length(Tspan)
        T = Tspan(i);
        Tind = (T-120)/10; % find Tind to access Exc data
        
        ExcRn = [Exc(f,Tind,1); Exc(f,Tind,2); Exc(f,Tind,3)]; %Exclude for given folder/Temp
        Tlogic = (V(:,1)==T); % logical when T = loop temperature
        Vredc = V(Tlogic,:); %Vredc is now only the loop temperature rows
        
        Vredc = Vredc.*ExcRn(1:1:size(Vredc,1)); %Vredc is now free from the 'skipped' data
        Vredc(Vredc==0) = []; %remove zeros for size calculation
        
        %% Averaging
        skip = false;
        
        if size(Vredc,1) > 1 % Only average if there's more than 1 parameter set...
            Aavg = mean(Vredc); % Is averaging Gaussian parameters the same as averaging data?
            Aavg(1) = T;
            Aavg(2) = NaN; %Run number no longer makes sense NaN so it can't accidentally be used
        elseif size(Vredc,1) > 0 %only pass on data if there are any entries at all
            Aavg = Vredc;
        end
        
        k = Aavg(3:end-1); % last entry is Hmelt
        
        if isempty(Vredc) == 1 % skip empties
           skip = true; 
        end
        
        
        if skip == false
        %% Find f, f1, f2 again  and plot
        tspan = 1:1/60:40; %timespan in minutes
        [fun, fun1, fun2] = SAGauss(tspan,k); % take the averaged values and find fun1 (poly curve)
        
        figure(f)
        hold on
        plot(tspan,fun1,'DisplayName',num2str(T))
        xlabel('Time in minutes')
        ylabel('f1 - Heat Flow')
        
        %% Convert from heat flows to conversion/rate of conversion
        %This can be achived by considering the relative cumulative
        %integral (cumINT/INT):
        
        B = cumtrapz(fun1)/trapz(fun1); %(This is where we could scale by DOCfrac
        dBdt = diff(B);%/(1/60); % in per minute %can remove 1/60 to have seconds as sampling rate is /1s
        dBdt = [0 dBdt];
        Tones = T.*ones(length(B),1);
        temp = [tspan' Tones B' dBdt'];

        DataVar = [DataVar; ...
                    temp]; %concatenate with temperature
        end     
    end
                  
             xdata = [DataVar(:,3) DataVar(:,2)]; % B, T
              
             RSFact = 1E8; %Rescale Factor
             ydata = DataVar(:,4); %dBdt
              
%              lb = [1E6, 1E4, 1E6, 1E4, 1.0, 1.0]; % why are these intial guess so large
%              p0 = [1E9, 1E5, 1E9, 1E5, 1.5, 1.5];
%              ub = [1E11, 1E6, 1E11, 1E6, 2.5, 2.5];
             
             lb = [1E6, 1E2, 1E6, 1E2, 1.0, 1.0]; % why are these intial guess so large
             p0 = [1E13, 1E5, 1E13, 1E5, 1.5, 1.5];
             ub = [1E18, 1E8, 1E18, 1E8, 2.5, 2.5];
%               
             opts = optimoptions(@lsqcurvefit);
             opts.MaxFunctionEvaluations = 1200; %double the evaluations
%               global T
%               T = CONC{i}(:,2)+273.15;
              [Pp, resnorm, residual] = lsqcurvefit(@kamal6,p0,xdata,ydata,lb,ub,opts); %,resnorm,residual,exitflag,output,lambda,jacobian
              resnorm = resnorm./length(ydata);
%               
              Pp = [Pp resnorm];
              LsqRes(f,:) = Pp.*[1./RSFact 1 1./RSFact 1 1 1 1]; %C1 C20P A1 Ea1 A2 Ea2 m n
              ymodel = kamal6(LsqRes(f,1:end),xdata);
%                            
              B = DataVar(:,3);
              T = DataVar(:,2);
%               
              figure(f+20)
              Tu = unique(T);
%               
              for Tind = 1:length(Tu)
              Tin = Tu(Tind).*ones(length(B),1);
              ymodel = kamal6(LsqRes(f,1:end),[B(2:end-1) Tin(2:end-1)]); 
                 hold on
%               plot(B(1:end-2),ymodel.*RSFact)
%               plot(B(1:5:end-2),ydata(1:5:end),'o')
%               
              end
              
        legend
end