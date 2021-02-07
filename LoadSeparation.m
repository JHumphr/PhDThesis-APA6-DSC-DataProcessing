close all
clear all
Vector = cell(1,12);
% cd('C:\Users\jarra\Documents\00_Thesis Work\03_DSC Data & Processing\DSC')
addpath('C:\Users\jarra\Documents\00_Thesis Work\03_DSC Data & Processing\DSC')
load('intPoints.mat')
foldConc = ones(10,2);

%% Conditionally Preinitialising/Loading Pan Mass for efficiency
    PMCheck = 0;
    if isfile('PanMass.mat')
    load('PanMass.mat');
    else
    PanMass = [];  
    end

%%
fold = dir;
DOCfrac = 0.999; % The TGA data for this also exists for some of the runs...
%% Open Folder
foldIDX =0;
for f = 1:1:12%length(fold)-11                                                 %loop through data folders
   if strcmp(fold(f).name,'.') == 1 || strcmp(fold(f).name,'..') == 1       %Exclude parent and current directory
   elseif fold(f).isdir == 0                                                %Exclude non-directories
   else
       foldIDX = foldIDX + 1;
       cd(fold(f).name)
       temp = char(fold(f).name);
       foldConc(foldIDX,1) = str2double(temp(6:8));
       foldConc(foldIDX,2) = str2double(temp(19:21));
       cd('.\matlab')                                                       %subfolder matlab which has the .txt version
       %% Load File
       datafile = dir;

       fileIDX = 0; %figure index
       for d = 1:1:length(datafile) %loop through data files
            if strcmp(datafile(d).name,'.') == 1 || strcmp(datafile(d).name,'..') == 1 || strcmp(datafile(d).name,'Vector.mat') == 1 || strcmp(datafile(d).name,'PanMass.mat') == 1  
%             || strcmp(datafile(d).name,'PanMass.mat') == 1 %elseif datafile(d).isdir == 0 %Exclude nondirectories
            else
                
            %% Load Raw Data using textscan    
            fileIDX = fileIDX + 1;                                          % increment figure index
            dfLength = length(datafile);
            
            fileID = fopen(datafile(d).name);
            hdrRows = 68;
            hdrData = textscan(fileID, '%s %s %s', hdrRows);
            S = textscan(fileID,'%f %f %f','HeaderLines',hdrRows);          %skip the 68 headerlines
            
            %% Load PanMass using textscan or load old dataset
            if isempty(PanMass) || PMCheck == 1
                idx = find(strcmp(hdrData{1,1},'Size')); %this slows everything down A LOT; consider extracting this info and storing in an array
                m = hdrData{1,2}(idx);
                m = str2double(m);
                PanMass(foldIDX,fileIDX) = m;
                PMCheck = 1;
            else
                m = PanMass(foldIDX,fileIDX);
            end
            fclose(fileID);
            
            clear idx hdrRows hdrData dfLength
            
            %% Extract metadata from Filename
            str = datafile(d).name(1:3); %take the first 3 characters which is the temperature
            str2 = datafile(d).name(5); %the fifth character is the run number
            Ti = str2double(str); % Isotherm temperature
            Rn = str2double(str2);
            PlotData{foldIDX,fileIDX}.Ti = Ti;
            PlotData{foldIDX,fileIDX}.Rn = Rn;
            PlotData{foldIDX,fileIDX}.filename = datafile(d).name;
            clear str str2
            
            t = S{1}(:);
            T = S{2}(:);
            hF = S{3}(:);
            hF = hF./m; %normalise all heat flow
            
            if isempty(hF) == 1; continue; end
            
            T(t <= 0) = NaN; 
            hF(t <= 0) = NaN;
            t(t <= 0) = NaN;
                
            clear S m 
            
            PlotData{foldIDX,fileIDX}.Raw = [t, hF];
            
            %% Extract melting peak
            % Only consider Tiso+10 and onwards
            logicTm = (T >= 175); % higher than 190 degree data samples might have 
            % trouble with this line
            
            Tm = T(logicTm);        % Remove below 175
            I = find(Tm >= 275,1);  % I = the FIRST index where Tm >= 270
            Tm = Tm(1:I);           % To first entry > 275 C
            
            tm = t(logicTm);        % Remove below 175
            tm = tm(1:I);           % To first entry > 275 C
            
            hFm = hF(logicTm);      % Remove below 175
            hFm = hFm(1:I);         % To first entry > 275 C
            
            clear I logicTm
            %% Integrate melt peak area
            x = Tm; %linear section first
            y = -hFm; %make endotherm up


            %% Two Region Method 
            %this is probably a step which would appreciate a better algorithm
            %indeed 1/6/19, this Tlow (180) temperatures gives trouble for some
            %170 degree datasets...
            Tlow = 180; Thigh = 260; %integration range
            
            %Exceptions
            % For 170 degree datasets add a Tlow of 185
            if abs(Ti - 170) <= 1; Tlow = 185; end
            
            ylinP1 = y(abs(x-Tlow) <= 1.5);  % Region 1: ~180 C
            ylinP2 = y(abs(x-Thigh) <= 1.5); % Region 2: ~260 C

            ylin = [ylinP1(1); ylinP2(1)];
            xlin = [Tlow; Thigh];
            
            Xlin = [ones(length(xlin),1) xlin];
            
            B = Xlin\ylin; % Linear Regression
            yline = B(2).*x+B(1);
            
            clear B xlin ylin  ylinP1 ylinP2 Xlin Tlow Thigh
            
            %% Calculate Hmelt
            IybVAL = trapz(x,y-yline); % Reduce integration range
            %units !! integral of T (in C) and Q (in mW) = mJC/s
            Hmelt = IybVAL;    %mJK/s * 60 s / 20 K  (heating rate)
            Hmelt = Hmelt*3;%./m; % mJ/mg == J/g % consider./m earlier (entire heat flow)
            % check Hmelt with Thermal Analysis Software
%             subplot(1,4,2)
            PlotData{foldIDX,fileIDX}.hFlowTrace = [x, y];
            PlotData{foldIDX,fileIDX}.hFlowBaseline = [x, yline];
            PlotData{foldIDX,fileIDX}.Hmelt = Hmelt;

            clear IybVAL yline
            
            %% Extract Reaction Cycle      
            % Integration points
            intP1 = intPoints{f}(d,1);
            intP2 = intPoints{f}(d,2); %time in minutes at integration point
%             disp(f)
%             disp(d)
            
%             disp(intP1)
%             disp(intP2)
%             %also remove data points at high times (>intP2 minutes)
         
            logicti = (t <= intP2);                                      %
            tRxn = t(logicti);                                              % I recently swapped this with the below
            TRxn = T(logicti);                                              % separation step, may need to undo
            hFRxn = hF(logicti);                                            %
            
            logicTi = (abs(TRxn-Ti) <= 3); %logical expression for data points when T ~= isothermal T
            tRxn = tRxn(logicTi);
            TRxn = TRxn(logicTi);
            hFRxn = hFRxn(logicTi);
            
            clear logicTi logicti
            %% Create Polymerisation Curve Baseline
%             % Algorithmic approach
%             tRxnBL =  tRxn(tRxn >= 41 & tRxn <= 41.5); % Use last 30 seconds of cycle
%             hFRxnBL = hFRxn(tRxn >= 41 & tRxn <= 41.5);
                        
            tRxnBLP1 = tRxn(tRxn >= intP1-0.1 & tRxn <= intP1+0.1);
            tRxnBLP2 = tRxn(tRxn >= intP2-0.1 & tRxn <= intP2+0.1);
            
            hFRxnBLP1 = hFRxn(tRxn >= intP1-0.01 & tRxn <= intP1+0.01);
            hFRxnBLP2 = hFRxn(tRxn >= intP2-0.01 & tRxn <= intP2+0.01);
            
            %starting point exception
            if mean(hFRxnBLP1) > 0 % this exception salavages some data that
                                   % may be otherwise cut-off [questionable whether we want this]
               hFRxnBLP1 = zeros(length(hFRxnBLP1),1);
            end
            
            tRxnBL =  [mean(tRxnBLP1)  mean(tRxnBLP2)]';
            hFRxnBL = [mean(hFRxnBLP1) mean(hFRxnBLP2)]';
            
            % Regress baseline
            tRxnBL2 = [ones(length(tRxnBL),1) tRxnBL];
            B = tRxnBL2\hFRxnBL; %baseline gradient/intercept
            
            heatFlowRxnBLine = B(2).*tRxn + B(1);
            tRxnNR = tRxn;
            
            %% Subtract Baseline
            hFRxnR = hFRxn - heatFlowRxnBLine;
            
            
            %Remove Negative Reaction
            hFRxnR = hFRxnR(hFRxnR >= 0); % negative reaction doesn't make sense...
            tRxnR = tRxn(hFRxnR >= 0);
            
%             %Remove Reaction Before and after integration points
%             hFRxnR = hFRxnR(tRxnR >= intP1);
%             tRxnR  = tRxnR(tRxnR  >= intP1);
            
            hFRxnR = hFRxnR(tRxnR <= intP2);
            tRxnR  = tRxnR(tRxnR  <= intP2); 
            
            Htotal = trapz(tRxnR,hFRxnR).*60;
            
%             if isempty(tRxnR) == false %this is the code which resets onset to t = 0
%             tRxnR = tRxnR-tRxnR(1);
%             end
            
            PlotData{foldIDX,fileIDX}.hFlowRxn = [tRxn, hFRxn];
            PlotData{foldIDX,fileIDX}.hFlowRxnBL = [tRxn,heatFlowRxnBLine]; 
            PlotData{foldIDX,fileIDX}.Htotal = Htotal;

            PlotData{foldIDX,fileIDX}.hFRxnRed = [tRxnR,hFRxnR];

%             clear B tRxnBL2 tRxnBL
            %% Data Process
            x = tRxnR(1:end-15)-tRxnR(1); % note tspan is irregular
%             disp(x)
            y = hFRxnR(1:end-15); %drop last 20 entries (mW)           

            weightFunc = 1;%((y+0.1)./max(y+0.1)); %scaling by y/ymax emphasises the phenomena; 0.1 is the floor weighting
%             weightFunc = 1;
            
            format short g  
            func = @(k)sum((weightFunc).*(((k(1)*(exp(-(x-k(2)).^2/(2*k(3)^2)).*(x<=k(2))+exp(-(x-k(2)).^2/(2*k(4)^2)).*(x>k(2)))  + ...  
                             k(5)*(exp(-(x-k(6)).^2/(2*k(7)^2)).*(x<=k(6))+exp(-(x-k(6)).^2/(2*k(8)^2)).*(x>k(6)))) - y).^2));            
            
                         
            %     ||-Gauss1-| |-Gauss2-||
            k0 = [0.6 7 1 1 0.5 8 1 1];
            
            %% Define k0 intelligently
            diffy = diff(y);
            i_cp1 = find(diffy <=0,1);
            k0(1) = y(i_cp1); %the first critical point
            k0(2) = x(i_cp1);
            i1_hwhm = find((k0(1)/2-y<=0),1);
            t1_hwhm = t(i1_hwhm);
            k0(3) = (k0(2)-t1_hwhm)/(sqrt(2*log(2))); %from gaussian func. Full width half max calc.
            k0(4) = k0(3);                            %approximate gaussian as symmetric at first
            
%             diffy = diffy(y>=0.01);
%             disp(diffy)
            i_cp2 = find(diffy >= 1E-5,1,'last');
%             disp(i_cp2)
            k0(5) = y(i_cp2); %the last critical point
            k0(6) = x(i_cp2);
            if k0(6) <= k0(2)
               k0(6) = k0(2)+k0(4).*1.5; %move it right by 1.5*the RMS  
               k0(5) = k0(1).*0.65; 
            end
            
            i2_hwhm = find((k0(5)/2-y<=0),1,'last');
            t2_hwhm = t(i2_hwhm);
            k0(8) = abs(t2_hwhm-k0(6))/(sqrt(2*log(2))); %from gaussian func. Full width half max calc.
            k0(7) = k0(4); %k0(8); %approximate gaussian as symmetric at first
            
            PlotData{foldIDX,fileIDX}.k0 = k0;
                      
            %% 

            Kmin= [ 0  0  1E-3 1E-3   0   0   1E-3   1E-3 ];
            Kmax= [inf 40 inf inf   inf   40  10  10];
            
            options = optimoptions('fmincon','Display','off');
            [k,fval,exitflag] = fmincon(func,k0,[],[],[],[],Kmin,Kmax,@(k)CONSTRAIN(x,k,Hmelt),options); % I can get goodness of fit indicators from here...
            %k = fmincon(f,k0,[],[],[],[],[],[],@ys); 
            
            %% Added to calculate R2
            TSS = sum(y - (mean(y)).^2);

            %%
            RSS = sum(((k(1)*(exp(-(x-k(2)).^2/(2*k(3)^2)).*(x<=k(2))+exp(-(x-k(2)).^2/(2*k(4)^2)).*(x>k(2)))  + ...  
                             k(5)*(exp(-(x-k(6)).^2/(2*k(7)^2)).*(x<=k(6))+exp(-(x-k(6)).^2/(2*k(8)^2)).*(x>k(6)))) - y).^2); %unweighted error function
            standardError = sqrt(RSS./(length(x) - 8)); %careful that this is now a weighted standard error
            R2 = 1 - RSS/TSS;

            f1c = k(1)*(exp(-(x-k(2)).^2/(2*k(3)^2)).*(x<=k(2))+exp(-(x-k(2)).^2/(2*k(4)^2)).*(x>k(2)));  
            f2c = k(5)*(exp(-(x-k(6)).^2/(2*k(7)^2)).*(x<=k(6))+exp(-(x-k(6)).^2/(2*k(8)^2)).*(x>k(6)));  
            
            %% Plot & Parameters
            PlotData{foldIDX,fileIDX}.gaussCurve = [x, f1c, f2c, f1c+f2c]; %actual values
            PlotData{foldIDX,fileIDX}.meta = [foldIDX, fileIDX, f, d]; %general meta data dump, f,d are the loop indicies as defining int points
            PlotData{foldIDX,fileIDX}.gaussParameters = [Ti Rn k standardError R2];
            PlotData{foldIDX,fileIDX}.exitFlag = exitflag;
%             Vector{end + 1, :} = [Ti Rn k Hmelt fval]; %Results Matrix
    
            end
     %%   Code for adding intPoints
     
     for i=1:1 %Supressed intPoints Code
%         disp(f)
%         disp(d)
%         figure(1)
%         grid
%         grid minor
%         figure(3)
%         grid
%         grid minor
%         axis([4 12 -0.5 0.5])
%         intPoints{f}(d,1) = input('t0');
%         axis([35 55 -0.5 0.5])
%         intPoints{f}(d,2) = input('tf');
     end
       end
       
       %% Once the above has been finished for the entire folder
% plot to check

% tspan = 0:1/60:40;

% Tprev = 0;
%       for i = 1:size(Vector{f-2},1)
%          k = Vector{f-2}(i,3:10); %parameters
%          Ti = Vector{f-2}(i,1); %iso temp
%          Rn = Vector{f-2}(i,2); %run number
%          [func, func1, func2] = SAGauss(tspan,k);
%       if abs(Ti-Tprev) >= 1
% %             figure
%       end
      %uncomment for plots
%          hold on
% %          plot(tspan,func1,'DisplayName',strcat(num2str(Rn),'-P'))
% %          axis([0 40 -0.8 0.8])
% %          plot(tspan,func2,'DisplayName',strcat(num2str(Rn),'-C'))
% %          plot(tspan,func,'-','DisplayName',strcat(num2str(Rn),'-Total'))
% %          plot(x,y,'o')
%          %also plot raw data
%          str = num2str(Ti);
%          str2 = num2str(f-2); 
%          str = [str ' ' str2];
%          title(str)
% %          legend
%          Tprev = Ti;

%       end

       cd ..
       %%
      cd ..
   end
end
%%
clear B
% save('Vector.mat','Vector')
save('PanMass.mat','PanMass')
save('PlotData.mat','PlotData')
save('foldConc.mat','foldConc')
% save('intPoints.mat','intPoints')
disp('PlotData.mat now contains a cell array of the reduced dataset parameters')

function [c,ceq] = CONSTRAIN(x,k,Hmelt)
% in nonlinfit ceq is minimised. c can be used for inequalities
c1 =  k(2) - k(6); % poly peak is before crys peak (bp - bc) < 0; bp < bc

% itk2 = find(x - k(2) <= 0.1,1);  % I = the FIRST index where Tm >= 270
% xc2 = x(1:itk2);  %(this is too ineffienct)
% clear itk2
xc2 = (0:k(2)/3:k(2))';
c2 = - (k(1)*(exp(-(xc2-k(2)).^2/(2*k(3)^2)))) + (k(5)*(exp(-(xc2-k(6)).^2/(2*k(7)^2)))); % until t = k(2) poly > crys

xc3 = (k(6)+2*k(8))'; 
c3 =  (k(1)*(exp(-(xc3-k(2)).^2/(2*k(4)^2)))) - (k(5)*(exp(-(xc3-k(6)).^2/(2*k(8)^2)))); % after t = k(6) crys > poly

c = [c1; c2; c3];
%% Crystallinity from Hmelt Method
Ifc = (quadgk(@(x)(k(5)*(exp(-(x-k(6)).^2/(2*k(7)^2)).*(x<=k(6))+exp(-(x-k(6)).^2/(2*k(8)^2)).*(x>k(6)))),0,40));%integrate from 0 or from first time in the isotherm 
ceq = Hmelt - Ifc*60;  % Hmelt is in J/g; If is the integral of time x heat; W*min/g   |60s/min| 

%% Polymerisation Enthalpy Method
%Alternate condition of the area under poly curve = lit value (-144
%Karger-Kocsis) would look something like:
% Ifp = (quadgk(@(x)k(1)*(exp(-(x-k(2)).^2/(2*k(3)^2)).*(x<=k(2))+exp(-(x-k(2)).^2/(2*k(4)^2)).*(x>k(2))),0,40)); 
% Hpoly = 144; % (Karger-Kocsis)
% ceq = Hpoly - Ifp*60; % + Hmelt - Ifc*60;

end
