function varargout = KamalFitGUI(varargin)
% KAMALFITGUI MATLAB code for KamalFitGUI.fig
%      KAMALFITGUI, by itself, creates a new KAMALFITGUI or raises the existing
%      singleton*.
%
%      H = KAMALFITGUI returns the handle to a new KAMALFITGUI or the handle to
%      the existing singleton*.
%
%      KAMALFITGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KAMALFITGUI.M with the given input arguments.
%
%      KAMALFITGUI('Property','Value',...) creates a new KAMALFITGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before KamalFitGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to KamalFitGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help KamalFitGUI

% Last Modified by GUIDE v2.5 07-Jun-2019 01:06:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KamalFitGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @KamalFitGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before KamalFitGUI is made visible.
function KamalFitGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to KamalFitGUI (see VARARGIN)
PlotData = load('PlotData.mat');
handles.PlotData = PlotData.PlotData;
ExcLoad = load('dataExclusion.mat');
handles.ExcLoad = ExcLoad.Exc;
handles.sortedGaussParamters = cell(3,5);

% Predefine figure symbols and colours here
handles.colorTriplet = [[0 0 1]; [0 0.5 0]; [1 0 0]; [0 0.75 0.75]; [0.75 0 0.75]];%[[146 0 0]; [146 73 0]; [219 209 1]; [36 255 36]; [219 255 109]]/255;
% colorTriplet = 
handles.shapeSet = {'o'; 'v'; 's'; '^'; 'd'};
handles.lineSet = {'-'; '--';'-.';':';'-'};

% Initialise all Axes Here
grid(handles.axes2,'on'); grid(handles.axes9,'on'); grid(handles.axes10,'on'); grid(handles.axes11,'on'); grid(handles.axes12,'on');
handles.axes2.YMinorGrid ='on'; handles.axes9.YMinorGrid ='on'; handles.axes10.YMinorGrid ='on'; handles.axes11.YMinorGrid ='on'; handles.axes12.YMinorGrid ='on';
handles.axes2.XMinorGrid ='on'; handles.axes9.XMinorGrid ='on'; handles.axes10.XMinorGrid ='on'; handles.axes11.XMinorGrid ='on'; handles.axes12.XMinorGrid ='on';
xlabel(handles.axes2,'time (minutes)');xlabel(handles.axes9,'time (minutes)');xlabel(handles.axes10,'time (minutes)');xlabel(handles.axes11,'time (minutes)');xlabel(handles.axes12,'time (minutes)');
ylabel(handles.axes2,'Heat Flow (J/g)');ylabel(handles.axes9,'Heat Flow (J/g)');ylabel(handles.axes10,'Monomer Conversion (B)');ylabel(handles.axes11,'Rate of Reaction (dB/dt)');ylabel(handles.axes12,'Fractional Crystallinity (a/aeq)');

% Main Driver...
popupmenu1_Callback(hObject,eventdata,handles)

% Choose default command line output for KamalFitGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes KamalFitGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = KamalFitGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Determine the selected data set
str = get(handles.popupmenu1,'String');
val = get(handles.popupmenu1,'Value');
foldIDX = val;

%Set current data to the selected data set
cla(handles.axes2); cla(handles.axes9); cla(handles.axes10); cla(handles.axes11); cla(handles.axes12);

Exc = squeeze(handles.ExcLoad(foldIDX,:,:))';

uit1 = uitable(handles.uitable1);
uit1.Data = double(Exc)';
uit1.RowName = {'130 C', '140 C', '150 C', '160 C', '170 C'};

handles.foldIDX = val;
handles.Exc = Exc; 

handles.sortedGaussParameters = gaussPlotHeatMap(hObject, eventdata, handles);
handles.DataVar = averageGauss(hObject, eventdata, handles);

KamalFit(hObject, eventdata, handles);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.checkbox2,'Value',0)
set(handles.checkbox3,'Value',0)
% Hint: get(hObject,'Value') returns toggle state of checkbox1

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.checkbox1,'Value',0)
set(handles.checkbox3,'Value',0)
% Hint: get(hObject,'Value') returns toggle state of checkbox2

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.checkbox2,'Value',0)
set(handles.checkbox1,'Value',0)
% Hint: get(hObject,'Value') returns toggle state of checkbox3

% --- Plots raw Polymerisation Gauss Curves & Presents Gauss Parameters in
% tables, create kV a 3x5 matrix (Rn x Temp) for current folder (set by
% popout1
%% --- Extracts gaussParameters/Sorts Data
function gaussP_CellArr = gaussPlotHeatMap(hObject, eventdata, handles) %this is no longer a heat map... dont want to change the name for fear of breaking things...
Exc = handles.Exc;
PlotData = handles.PlotData;
foldIDX = handles.foldIDX;
x = 0:0.1:40; %tspan
%%
 PlotData{4,12} = [];                                                                                                % <------------------ Hard coded data correction
%%


axes(handles.axes2)
for i = 1:12
%     disp(PlotData{foldIDX,i})
    if isempty(PlotData{foldIDX,i}) == 0 % If the data set exists then:
    gaussP_Vector = PlotData{foldIDX,i}.gaussParameters;
    
    % add these to end of gaussParameters
    Hmelt = PlotData{foldIDX,i}.Hmelt;
    Htotal = PlotData{foldIDX,i}.Htotal;
    
    gaussP_Vector(end+1) = Hmelt;
    gaussP_Vector(end+1) = Htotal;
    
    Ti = gaussP_Vector(1);
    TIDX = (Ti - 120)/10;
    Rn = gaussP_Vector(2);
    gaussP = gaussP_Vector(3:10);
    stdError = gaussP_Vector(11);
    
    if Exc(Rn,TIDX) == 1 % If we're not excluding the data set then:
        [~,f1,~] = SAGauss(x,gaussP);
        hold on
        plot(x,f1); hold off
        gaussP_CellArr{Rn,TIDX} = gaussP_Vector;
    else % If we're excluding the data set then:
        gaussP_CellArr{Rn,TIDX} = NaN(1,length(gaussP_Vector)); %[NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
    end
    end
end

%% This loop groups data into Isothermal Temperatures
for i = 1:size(gaussP_CellArr,2) %for each temperature IDX
    tablehandle = sprintf('uitable%d',i+1);
    uiti = uitable(handles.(tablehandle));
        for j = 1:size(gaussP_CellArr,1) %for each run number (IDX)
            if isempty(gaussP_CellArr{j,i}) == 0
                gaussP_TempArr(j,:) = double(gaussP_CellArr{j,i});
            else
                gaussP_TempArr(j,:) = NaN(14,1); % this might cause errors when adding parameters/etc to PlotData.gaussianParameters
            end  
        end
    uiti.Data = gaussP_TempArr(:,3:end);
end

gaussP_CellArr=gaussP_CellArr;

%% --- Averages Gaussian Paramters
    function DataVar = averageGauss(hObject, eventdata, handles)
gaussP_Sorted = handles.sortedGaussParameters;
%                 130     140      150      160       170
%                 130     140      150      160       170
colorTriplet = [[0 0 1]; [0 0.5 0]; [1 0 0]; [0 0.75 0.75]; [0.75 0 0.75]];

%colorTriplet = [[146 0 0]; [146 73 0]; [219 209 1]; [36 255 36]; [219 255 109]]/255;
shapeSet = {'o'; 'v'; 's'; '^'; 'd'};
lineSet = {'-'; '--';'-.';':';'-'};
%                 Red    Green     Cyan     Blue   Magenta

% disp(Kv)
x = 0:0.1:40;
axes(handles.axes9)
DataVar = [];
for i = 1:size(gaussP_Sorted,2) %for each temperature
    gaussP_TempArr2 = cell2mat(gaussP_Sorted(:,i));
%     disp(gaussP_TempArr2)
    if isempty(gaussP_TempArr2) == 1; gaussP_TempArr2 = [NaN NaN NaN]; end %weird exclusion for empties...

    if mean(isnan(gaussP_TempArr2)) == 1 %if all are NaN so are the output paramters
        gaussP_AvgTempArr2 = NaN(1,14); % [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
    else
        for k = 1:length(gaussP_TempArr2) % for each parameter
%             disp(k)
            gaussP_kVec = gaussP_TempArr2(:,k);
            gaussP_kVec(isnan(gaussP_kVec)) = [];
            gaussP_AvgTempArr2(k) = mean(gaussP_kVec);
        end

    end
    gaussP_AvgTempArr2 = gaussP_AvgTempArr2(3:end);

    disp(gaussP_AvgTempArr2)
    Ktable(i,:) = gaussP_AvgTempArr2(1:end);
%     disp(Ktable)
    [f,f1,f2] = SAGauss(x,gaussP_AvgTempArr2(1:8));
    
    axes(handles.axes9)
    hold on
    plot(x(1:2:end),f1(1:2:end),'LineStyle',lineSet{i},'LineWidth',1.5,'Color','k'); hold off % swap colorTriplet(i,:) for 'k' for black plots
 
%% Convert Gaussians to Phenomena
B = cumtrapz(f1)/trapz(f1);%*(6/144); %This is where we could scale by DOCfrac   %remove /trapz(f1) %60 is from sampling rate, 144 is from Karger-Kocsis %* /trapz(f1);%
    if sum(B>1) > 0
        disp('Warning: Conversion Greater than 1 set to 1')
    end
B = B./max(B); %stop gap to prevent weird results for Hpoly > 144 [...]

a = cumtrapz(f2)/trapz(f2); %Relative crystallinity a/aeq
dBdt = diff(B)./diff(x);%/(1/60); % in per minute %can remove 1/60 to have seconds as sampling rate is /1s
disp(trapz(f1))
dBdt = dBdt.*(trapz(f1).*6./144);
dBdt = [dBdt(1) dBdt]; % here 6/144 scales by are underthe curve and we've removed the /trapz(f1) normalisation from the B calculation.
T = i.*10+120;
onesT = T.*ones(length(B),1);
temp = [x' onesT B' dBdt' a'];

axes(handles.axes10)
hold on; plot(x(1:2:end),B(1:2:end),'LineStyle',lineSet{i},'LineWidth',1.5,'Color',colorTriplet(i,:)); hold off % colorTriplet(i,:) or 'k'

axes(handles.axes11)
hold on; plot(x(1:2:end),dBdt(1:2:end),'LineStyle','none','Marker',shapeSet{i},'Color',colorTriplet(i,:)); hold off % colorTriplet(i,:) or 'k'

DataVar = [DataVar; temp];
    
end

assignin('base','Ktable',Ktable)
axes(handles.axes9); legend('130 C', '140 C', '150 C', '160 C', '170 C'); title('Averaged Heat Flow Gaussians')
axes(handles.axes10); 
legend('130 C', '140 C', '150 C', '160 C', '170 C'); 
title('Conversion')
axes(handles.axes11); legend('130 C', '140 C', '150 C', '160 C', '170 C'); title('Polymerisation Model Fitting')
uitT = uitable(handles.uitable7);
uitT.Data = Ktable(:,1:4); 

set(handles.uitable7,'UserData',Ktable)

%% Crystallisation Gaussian Parameters are converted here just not pushed to the GUI


%% Plot Export
              axes(handles.axes10)
              savefigure(handles,'Conv')
              
              axes(handles.axes9)
              savefigure(handles,'AvgGauss')

%% --- Executes Kamal & Kim Model Fittings
function KamalFit(~, ~, handles)
%% Define handles

            colorTriplet = handles.colorTriplet;
            shapeSet = handles.shapeSet;
            lineSet = handles.lineSet;

    
%% Initialise Variables [t T B dBdt a/a0] 
            DataVar = handles.DataVar; % from DataVar Created by averageGauss
            B = DataVar(:,3);
            T = DataVar(:,2);
            dBdt = DataVar(:,4);
            t = DataVar(:,1);
            a = DataVar(:,5);

            B(isnan(dBdt)) = [];
            T(isnan(dBdt)) = [];
            t(isnan(dBdt)) =[];
            a(isnan(dBdt)) =[];
            dBdt(isnan(dBdt)) = [];

            
            T(isnan(B)) =[];
            dBdt(isnan(B)) =[];
            t(isnan(B)) = [];
            a(isnan(B)) =[];
            B(isnan(B)) =[];

            RSFact = 1;%E6; % rescales to avoid underflow
            ydata = dBdt.*RSFact;

            xdata = [B T t]; % B, T            

%% Get checkbox values to contrains Kamal Fitting Parameters      
            freeCheckbox            = get(handles.checkbox5,'Value');

            fixedeaCheckbox         = get(handles.checkbox4,'Value');
            fixedidxCheckbox        = get(handles.checkbox6,'Value');
            fixedidxeaCheckbox      = get(handles.checkbox7,'Value');
            
            fixidxeaINPUTCheckbox   = get(handles.checkbox8,'Value');
            
            fixedAvgeaCheckbox      = get(handles.checkbox9,'Value');
            fixedAvgidxCheckbox     = get(handles.checkbox10,'Value');
            fixedAvgbothTCheckbox   = get(handles.checkbox11,'Value');
            
            % need to add checkboxes for other fit conditions

            fitConditions = {'Free'; 'Fixed Ea Indicies'; 'Fixed Ea'; 'Fixed index'};
            clear lb lbloop p0 p0loop ub ubloop AI kamalParam resnorm residual kamalParamRN
            fitIDX = [];

             % double the number of evaluatoins              
            opts = optimoptions(@lsqcurvefit);
            opts.MaxFunctionEvaluations = 2400; %double the evaluations

%% Parameter Constraints for Kamal Model Fitting
for ParameterConstraints = 1
%Free Parameters
            axes(handles.axes11)
            % values before March2020
if freeCheckbox == 1
            % all paramters are variables
            %          [ A1  Ea1  A2   Ea2   m    n ]
             lb(1,:) = [1E4, 2E4, 1E4, 2E4, 0.5,  0.5]; %lower bound
             p0(1,:) = [1E6, 6E4, 1E6, 5E4, 1.55,  1.21]; %initial value
             ub(1,:) = [1E7, 7E4, 1E7, 7E4, 2.0,  2.0]; %upper bound
             fitIDX = [fitIDX 1];
end


% if freeCheckbox == 1
%             all paramters are variables
%                      [ A1  Ea1  A2   Ea2   m    n ]
%              lb(1,:) = [1E3, 2E4, 1E3, 2E4, 0.5,  0.5]; %lower bound
%              p0(1,:) = [1E6, 6E4, 1E6, 5E4, 1.55,  1.21]; %initial value
%              ub(1,:) = [1E7, 1E5, 1E7, 1E5, 2.0,  2.0]; %upper bound
%              fitIDX = [fitIDX 1];
% end


%Fixed Ea to Teuwen
if fixedeaCheckbox == 1
             lb(2,:) = [1E4, 6.86E4, 1E4, 5.94E4, 0.5, 0.5]; % why are these intial guess so large
             p0(2,:) = [1E5, 6.86E4, 1E5, 5.94E4, 1, 1];
             ub(2,:) = [1E7, 6.86E4, 1E7, 5.94E4, 2.0, 2.0];
             fitIDX = [fitIDX 2];
end

%Fixed m,n to Teuwen
if fixedidxCheckbox == 1
            % all paramters are variables
            %     [ A1   Ea1  A2  Ea2   m   n ]
             lb(3,:) = [1E4, 2E4, 1E4, 2E4, 1.55, 1.21]; % why are these intial guess so large
             p0(3,:) = [1E5, 6E4, 1E6, 5E4, 1.55, 1.21];
             ub(3,:) = [1E7, 7E4, 1E7, 7E4, 1.55, 1.21];
             fitIDX = [fitIDX 3];
end

%
if fixedidxeaCheckbox == 1
            % all paramters are variables
            %           [ A1   Ea1  A2  Ea2   m   n ]
             lb(4,:) = [1E4, 6.86E4, 1E4, 5.94E4, 1.55, 1.21]; % why are these intial guess so large
             p0(4,:) = [1E6, 6.86E4, 1E6, 5.94E4, 1.55, 1.21];
             ub(4,:) = [1E7, 6.86E4, 1E7, 5.94E4, 1.55, 1.21];
             fitIDX = [fitIDX 4];
end

if fixidxeaINPUTCheckbox == 1
            % all paramters are variables
            %     [ A1   Ea1  A2  Ea2   m   n ]
            DataArray = get(handles.uitable9,'Data');
            A1 = DataArray(1); A1l = A1; A1u = A1;
            A2 = DataArray(2); A2l = A2; A2u = A2;
            Ea1 = DataArray(3); Ea1l = Ea1; Ea1u = Ea1;
            Ea2 = DataArray(4); Ea2l = Ea2; Ea2u = Ea2;
            m =   DataArray(5); ml = m; mu = m;
            n =   DataArray(6); nl = n; nu = n;
            % if 0 set to normal bounds
            if A1 == 0;  A1l = 1E4;  A1 = 1E6;  A1u = 1E8; end
            if A2 == 0;  A2l = 1E4;  A2 = 1E6;  A2u = 1E8; end
            if Ea1 == 0; Ea1l = 1E4; Ea1 = 4E4; Ea1u = 7.5E4; end
            if Ea2 == 0; Ea2l = 1E4; Ea2 = 4E4; Ea2u = 7.5E4; end
            if m   == 0; ml = 0.5; m = 1.0; mu = 2.0; end
            if n   == 0; nl = 0.5; n = 1.0; nu = 2.0; end
            
            lb(8,:) = [A1l,  Ea1l, A2l, Ea2l, ml, nl]; 
            p0(8,:) = [A1,  Ea1,  A2, Ea2,  m,  n];
            ub(8,:) = [A1u, Ea1u, A2u, Ea2u, mu, nu];

            fitIDX = [fitIDX 8];
end
 %% these must be manually defined!
%Fixed Ea to Avg Value
if fixedAvgeaCheckbox == 1
    tempEa1Avg = 6.41E4;
    tempEa2Avg = 5.06E4; 
             lb(5,:) = [1E4, tempEa1Avg, 1E4, tempEa2Avg, 0.5, 0.5]; 
             p0(5,:) = [1E5, tempEa1Avg, 1E5, tempEa2Avg, 1, 1];
             ub(5,:) = [1E7, tempEa1Avg, 1E7, tempEa2Avg, 2.0, 2.0];
             fitIDX = [fitIDX 5];
end

%Fixed m,n to Avg Value
if fixedAvgidxCheckbox == 1
    tempmAvg = 1.148;
    tempnAvg =  1.077;
            % all paramters are variables
            %           [ A1   Ea1  A2  Ea2   m   n ]
             lb(6,:) = [1E4, 2E4, 1E4, 2E4, 0.83, 0.71]; 
             p0(6,:) = [1E5, 5E4, 1E5, 5E4, 0.83, 0.71];
             ub(6,:) = [1E7, 7E4, 1E7, 7E4, 0.83, 0.71];
             fitIDX = [fitIDX 6];
end

% both to AvgValues
if fixedAvgbothTCheckbox == 1
    tempmAvg = 1.148;
    tempnAvg =  1.077;
    tempEa1Avg = 6.84E4;
    tempEa2Avg = 5.02E4; 
            % all paramters are variables
            %           [ A1   a1  A2  Ea2   m   n ]
             lb(7,:) = [1E4, tempEa1Avg, 1E4, tempEa2Avg, tempmAvg, tempnAvg]; 
             p0(7,:) = [1E5, tempEa1Avg, 1E5, tempEa2Avg, tempmAvg, tempnAvg];
             ub(7,:) = [1E7, tempEa1Avg, 1E7, tempEa2Avg, tempmAvg, tempnAvg];
             fitIDX = [fitIDX 7];
    clear tempmAvg tempnAvg
end
end

%% Kamal Fitting for Each Set of Parameter Constrains
for fitIdx = fitIDX % for each set of fit conditions run model fitting
% grab bounds and p0 from arrays using fitIdx; Rescale
     AI = [RSFact, 1, RSFact, 1, 1, 1];
     BI = [1E-4, 1E-4, 1E-4, 1E-4, 1, 1];
     lbloop = BI.*AI.*lb(fitIdx,:);
     p0loop = BI.*AI.*p0(fitIdx,:);
     ubloop = BI.*AI.*ub(fitIdx,:);
    
% Regression Weighting
      weightFunc = 1;%(1-xdata(:,1)).^4;%ydata.*(1-xdata(:,1)).^4;%;./max(ydata);

% % % Unweighted Version: (lsqcurvefit, unweighted)
%     
%     [kamalParam, resnorm, ~] = lsqcurvefit(@kamal6,p0loop,xdata,ydata,lbloop,ubloop,opts); %,resnorm,residual,exitflag,output,lambda,jacobian
%      resnormUnweighted = resnorm;     
%                 
% % Weighted Version: (lsqcurvefit, unweighted)
    ydataweighted = ydata.*sqrt(weightFunc);
    
    kamal6weight2 = @(Param,xdata)(kamal6(Param,xdata,BI).*sqrt(weightFunc));
    [kamalParam, resnorm, ~] = lsqcurvefit(kamal6weight2,p0loop,xdata,ydataweighted,lbloop,ubloop,opts); %,resnorm,residual,exitflag,output,lambda,jacobian
     
% New Version (fmincon (or fminsearch), weighted)
% 
%       weightFunc = 1;
%       objFunc = @(Param)sum(weightFunc.*(ydata - kamal6(Param,xdata)).^2); 
%       
%       [kamalParam] = fmincon(@(Param)objFunc(Param), p0loop, [],[],[],[],lbloop,ubloop);
%       resnorm = sum((ydata-kamal6(kamalParam,xdata)).^2);
%       disp(resnorm) 
%       disp(resnormold)


      stdErr = resnorm/sqrt(length(ydata)-6);
      row = [kamalParam./BI stdErr];
      rowScaled = row.*[1./RSFact 1 1./RSFact 1 1 1 1./RSFact.^2]; %note the error must be scaled by the scalefactor squared due to the ^2 in the fitting algorithm

      % True Solution (ODE Solver)
      tspan = xdata(:,3);
      Tspan = xdata(:,2);
      [ym,xm] = kamal6ODEfunc(rowScaled(1:end-1),Tspan);
    
      hold on
      plot(xm(:,2),ym,'.','Color',[0.4 0.4 0.4])
       
      truestdErr = sum((ym-ydata).^2)/sqrt(length(ydata)-6);
      kamalParamRN(fitIdx,1:8) = [rowScaled, truestdErr];

%       kamalParamRN(fitIdx,:) = rowScaled;
      clear row rowScaled
      

      % Store parameters in uitable8
      uitKP = uitable(handles.uitable8);
      uitKP.Data = [kamalParamRN(fitIdx,1:2) kamalParamRN(fitIdx,5); kamalParamRN(fitIdx,3:4) kamalParamRN(fitIdx,6)];% kamalParamRN(fitIdx,8:11)];
end
% plot last one calculated     
              ymodel = kamal6(kamalParamRN(fitIDX,1:6),[B(2:end-1) T(2:end-1)],[]);
              
              hold on
              
              TU = unique(T(2:end-1));
              for Tu = (TU)' % for each unique temperature
                  TuIdx = (Tu-120)/10; % define index by formula
                  ymodelT = ymodel(T(2:end-1) == Tu);
                  plot(t(1:length(ymodelT)),ymodelT,'LineWidth',1.5,'LineStyle','-','Color','k');
                  %'LineStyle',lineSet{TuIdx},
                  legend('130 C - Data','140 C - Data', '150 C - Data', '160 C - Data', '170 C - Data', 'Model - Actual', 'Model - Apparent')
              end
              
             savefigure(handles,'KamalFit')

%% Kim-Avrami (Modified Avrami) Model Fitting
          p0c = [120 10 3 0]; % Initial Value
%           TempRange = unique(xdata(:,2));
          isoParamSet = [];

for ncset = 3:3 % 1:4 <-                      Can limit the ncd analysed here.
%     for TIso = TempRange'
%     Ti = (TIso-120)/10;
%                   p0c(4)=ncset;
%                   xdata_reduced = xdata(xdata(:,2) == TIso,:);
%                   a_reduced = a(xdata(:,2) == TIso,:);
% 
%                   [kimParam, resnorm] = lsqcurvefit(@kimAvrami,p0c,xdata_reduced,a_reduced,[1 1E-4 0 p0c(4)],[1 1 40 p0c(4)],opts);
% %                   assignin('base','XData',[xdata a])
%                   stdErr = resnorm/sqrt(length(a_reduced)-3);
% 
%                   uitAP = uitable(handles.uitable10);
%                   uitAP.Data = real(kimParam);
% 
%                   tempRow = [TIso kimParam stdErr];
%                   kimParamSet = [kimParamSet; tempRow];
% 
%                   ycmodel = kimAvrami(kimParam,xdata_reduced);
%                   tcmodel = xdata_reduced(:,3);
% 
%                   axes(handles.axes12)
%                   hold on
%                   title('Crystallisation Model Fitting')
%                   plot(tcmodel(1:5:end),ycmodel(1:5:end),'LineWidth',1.5,'LineStyle',lineSet{Ti},'Color','k')
% 
%                   plot(xdata_reduced(1:3:end,3),a_reduced(1:3:end),'LineStyle','none','Marker',shapeSet{Ti},'Color',colorTriplet(Ti,:)) %where does the perfectly linear dataset come from
%                   legend('Model','Gaussian','Location','SouthEast')
% 
%                   tempArray(Ti,:) = tempRow;
%     end
%                   kAParam(ncset,:,:) = tempArray; %Tiso [4 param] stdErr

%% New Section
                  p0c(3)=ncset;

                  [avramiParam, resnorm] = lsqcurvefit(@isoAvrami,p0c,xdata,a,[-inf -inf p0c(3) 0],[inf inf p0c(3) inf],opts);
                    disp(avramiParam)
                  stdErr = resnorm/sqrt(length(xdata)-3);

                  uitAP = uitable(handles.uitable10);
                  uitAP.Data = real(avramiParam);

                  tempRow = [avramiParam stdErr];
                  isoParamSet = [isoParamSet; tempRow];

                  ycmodel = isoAvrami(avramiParam,xdata);
                  tcmodel = xdata(:,3);

                  axes(handles.axes12)
                  hold on
                  title('Crystallisation Model Fitting')
%                   plot(tcmodel(1:5:end),ycmodel(1:5:end),'ko');%,'LineWidth',1.5,'LineStyle',lineSet{Ti},'Color','k')

%                   plot(xdata(1:3:end,3),a(1:3:end))%,'LineStyle','none','Marker',shapeSet{Ti},'Color',colorTriplet(Ti,:)) %where does the perfectly linear dataset come from
                  legend('Model','Gaussian','Location','SouthEast')
                  shapeSet = {'o'; 'v'; 's'; '^'; 'd'};

              TU = unique(T(2:end-1));
              for Tu = (TU)' % for each unique temperature
                  TuIdx = (Tu-120)/10; % define index by formula
                  
                  xdataT = xdata(xdata(2:end-1,2) == Tu,3);
                  aT = a(xdata(2:end-1,2) == Tu);
                  plot(xdataT(1:3:end),aT(1:3:end),'LineStyle','none','Marker',shapeSet{TuIdx},'Color',colorTriplet(TuIdx,:))%,'LineStyle','none','Marker',shapeSet{Ti},'Color',colorTriplet(Ti,:)) %where does the perfectly linear dataset come from
                  ycmodelT = ycmodel(T(2:end-1) == Tu);
                  plot(t(2:length(ycmodelT)),ycmodelT(2:end),'LineWidth',1.5,'LineStyle','-','Color',colorTriplet(TuIdx,:));
                  %'LineStyle',lineSet{TuIdx},
                   legend('130 C','','140 C','','150 C','', '160 C','','170 C','','Location','SouthEast')
              end

%                   tempArray(Ti,:) = tempRow;

                  isoAvramiParam(ncset,:,:) = tempRow; %Tiso [4 param] stdErr

%%
end
          axes(handles.axes12)
          savefigure(handles,'KimFit')
%% Set output handles
          set(handles.uitable10,'Data',isoParamSet)   
          handles.kamalParam = kamalParamRN;
          set(handles.uitable8,'UserData',kamalParamRN);
          set(handles.uitable10,'UserData',isoAvramiParam);  %what's up with this redundancy?
          
          
%% Button % Checkbox Callbacks              
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles) % Save
str = get(handles.popupmenu1,'String');
val = get(handles.popupmenu1,'Value');

gParam = get(handles.uitable7,'UserData');
kParam = get(handles.uitable8,'UserData');
kAParam = get(handles.uitable10,'UserData');
% disp(kParam)

foldID = val;

load('FitParamFull.mat');
% FitParamFull = FitParamFull.FitParamFull;

FitParamFull{foldID}.gParam = gParam;
FitParamFull{foldID}.kParam = kParam;
FitParamFull{foldID}.kAParam = kAParam;

save('FitParamFull.mat','FitParamFull');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles) % Open Parameter PLot
ParameterPlotGUI
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles) % Run and Save for all folders
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Determine the selected data set
str = get(handles.popupmenu1,'String');
val = get(handles.popupmenu1,'Value');

for foldIDX = 1:10
    %Set current data to the selected data set
    cla(handles.axes2); cla(handles.axes9); cla(handles.axes10); cla(handles.axes11); cla(handles.axes12);
%     set(handles.popupmenu1,'String',foldIDX);
    set(handles.popupmenu1,'Value',foldIDX);
    Exc = squeeze(handles.ExcLoad(foldIDX,:,:))';
    if sum(sum(Exc)) == 0; continue; end % skip run if exclusion matrix is empty
    
    uit1 = uitable(handles.uitable1);
    uit1.Data = double(Exc)';
    uit1.RowName = {'130 C', '140 C', '150 C', '160 C', '170 C'};

    handles.foldIDX = foldIDX;
    handles.Exc = Exc; 

    handles.sortedGaussParameters = gaussPlotHeatMap(hObject, eventdata, handles);
    handles.DataVar = averageGauss(hObject, eventdata, handles);

    KamalFit(hObject, eventdata, handles);
    
    gParam = get(handles.uitable7,'UserData');
    kParam = get(handles.uitable8,'UserData');
    kAParam = get(handles.uitable10,'UserData');

    load('FitParamFull.mat');
    % FitParamFull = FitParamFull.FitParamFull;

    FitParamFull{foldIDX}.gParam = gParam;
    FitParamFull{foldIDX}.kParam = kParam;
    FitParamFull{foldIDX}.kAParam = kAParam;

    save('FitParamFull.mat','FitParamFull');

end
guidata(hObject,handles);


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenu1_Callback(hObject, eventdata, handles)

function uitable10_KeyPressFcn(hObject, eventdata, handles)
% disp(eventdata.Key)
if strcmp(eventdata.Key,'return') == 1
popupmenu1_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11
%% Other Functions
function savefigure(handles,figName)
    set(0,'showhiddenhandles','on'); % Make the GUI figure handle visible
    axesH = gca;% findobj(gcf,'type','axes'); % Find the axes object in the GUI

    saveFigHandle = figure('visible','off');
    %               ax = axes;
    copyobj(axesH,saveFigHandle);             
    set(gca,'ActivePositionProperty','outerposition')
    set(gca,'Units','normalized')
    set(gca,'OuterPosition',[0 0 1 1])
    set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
    set(gcf,'PaperPosition',[.25 .25 16 12]);
    title('') % remove figure title
    legend('130 C','','140 C','','150 C','', '160 C','','170 C','','Location','SouthEast')
%     legend('130 C','140 C', '150 C', '160 C', '170 C', 'ODE - Actual', 'Model - Apparent')
    %                   legend('130 C','140 C', '150 C', '160 C', '170 C', 'Model')
    %               str = get(handles.popupmenu1,'String');
    foldIDX = get(handles.popupmenu1,'Value');
    foldIDXstr = num2str(foldIDX);
    fileSaveStr = strcat('figures/',figName,'_Folder',foldIDXstr,'.png'); % can change preferred image format here (i.e. epsc, fig, png)
    %               figure('PaperPosition',[.25 .25 8 6]);
    saveas(gcf,fileSaveStr)  
    cla(saveFigHandle)
    axes(axesH)
