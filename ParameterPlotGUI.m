function varargout = ParameterPlotGUI(varargin)
% PARAMETERPLOTGUI MATLAB code for ParameterPlotGUI.fig
%      PARAMETERPLOTGUI, by itself, creates a new PARAMETERPLOTGUI or raises the existing
%      singleton*.
%
%      H = PARAMETERPLOTGUI returns the handle to a new PARAMETERPLOTGUI or the handle to
%      the existing singleton*.
%
%      PARAMETERPLOTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMETERPLOTGUI.M with the given input arguments.
%
%      PARAMETERPLOTGUI('Property','Value',...) creates a new PARAMETERPLOTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ParameterPlotGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ParameterPlotGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ParameterPlotGUI

% Last Modified by GUIDE v2.5 30-Apr-2019 12:02:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ParameterPlotGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ParameterPlotGUI_OutputFcn, ...
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


% --- Executes just before ParameterPlotGUI is made visible.
function ParameterPlotGUI_OpeningFcn(hObject, eventdata, handles, varargin)

UpdatePlots(hObject, eventdata, handles)

% Update stuff
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes ParameterPlotGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ParameterPlotGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function UpdatePlots(hObject, eventdata, handles, Parameters, foldConc)
% Load data
load('FitParamFull.mat') % 10 Cell array with kamal param with b
load('foldConc.mat') % foldConc is a 12 x 2
R = 8.314;

Parameters = FitParamFull;

val = get(handles.popupmenu1,'Value');
fitConstrain = get(handles.popupmenu2,'Value');
table = [];

%% For First Set
if val == 1
% C1 Loop (C20P constant)
for foldIDX = 1:5
    %% Load Data for Specific Concentration Set
Pgauss = Parameters{foldIDX}.gParam;
Pkamal = Parameters{foldIDX}.kParam(fitConstrain,:);
xConc = foldConc(foldIDX,1);
xTicks = foldConc(1:5,1);
    %% Save and Export Table
table = [table; Pkamal];
save('table.mat', 'table')
    %% Plot Kamal Parameters    
% A1
axes(handles.axes1)
hold on
title('A1')
plot(xConc,Pkamal(1),'ro')
A1 = Pkamal(1);

% A2
axes(handles.axes2)
hold on
title('A2')
plot(xConc,Pkamal(3),'ro')
A2 = Pkamal(3);

% m
axes(handles.axes3)
hold on
title('m')
plot(xConc,Pkamal(5),'ro')

% Ea1
axes(handles.axes4)
hold on
title('Ea1')
plot(xConc,Pkamal(2),'ro')
Ea1 = Pkamal(2);

% Ea2
axes(handles.axes5)
hold on
title('Ea2')
plot(xConc,Pkamal(4),'ro')
Ea2 = Pkamal(4);

% n
axes(handles.axes6)
hold on
title('n')
plot(xConc,Pkamal(6),'ro')

% k1
axes(handles.axes7)
hold on
title('k1')
k1 = A1.*exp(-Ea1./R./(150+273));
plot(xConc,k1,'ro')

% k2
axes(handles.axes8)
hold on
title('k2')
k2 = A2.*exp(-Ea2./R./(150+273));
plot(xConc,k2,'ro') 
% safd
    %% Plot Gaussian Parameters  
        %% Initialise (Load Parameters) 
    % Initialise colour and symbols for plotting
    %                 130     140      150      160       170
    colorTriplet = [[1 0 0]; [0 1 0]; [0 0 1]; [1 1 0]; [1 0 1]];
    shapeSet =       {'o';     'v';     's';     '^';     'd'};
    %                 Red    Green     Blue     Yellow   Magenta

    xrange = [0.5 1.5];
    axes(handles.axes9) % old handle 
        %% b plot              
    figure(1)
    xlabelStr = 'C1 Concentration (mol%)';
            %% polymerisation  
    paramVal = 2;
    axisRange = [xrange 0 15];
    ylabelStr = 'b_p (min)';

    % plot
    sub1 = subplot(2,2,1);
    for Ti = 1:5
        hold on
        plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
    end % plot for each temperature; as black symbol defined from symbol list

    % dress figure
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
    axis(axisRange)
    xticks(xTicks)

    clear paramVal axisRange ylabelStr
            %% crystallisation 
    % axes(handles.axes10) %old handle
    paramVal = 6;
    ybound = 30;
    axisRange = [xrange 0 ybound];
    ylabelStr = 'b_c (min)';

    % plot
    sub3 = subplot(2,2,3);
    % title('Gaussian b - crystallisation')
    % plot with except for out of range values
    for Ti = 1:5
        hold on

        bc = Pgauss(Ti,paramVal);
        bcOutOfBounds = bc>ybound;
        bc(bc>=ybound) = ybound;

        plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        if bcOutOfBounds == 1 %plot a + as well
        plot(xConc,bc,'k+')    
        end
    end

    % dress
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub3.XMinorGrid = 'on'; sub3.YMinorGrid = 'on'; grid on; 
    axis(axisRange)
    xticks(xTicks)
    
    clear paramVal axisRange ylabelStr    
            %% polymerisation vs crystallisation
    figure(5)
    paramVal1 = 2;
    paramVal2 = 6;
%     axisRange = [xrange 0 15];
    ylabelStr = 'b_c (min)';

    % plot
%     sub1 = subplot(2,2,1);
    for Ti = 1:5
        hold on
        plot(Pgauss(Ti,paramVal1),Pgauss(Ti,paramVal2),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
    end % plot for each temperature; as black symbol defined from symbol list

    % dress figure
    xlabel('b_p (min)')
    ylabel(ylabelStr)
    sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
%     axis(axisRange)
%     xticks(xTicks)
            %% polymerisation/crystallisation
    figure(6)
    paramVal1 = 2;
    paramVal2 = 6;
    axisRange = [xrange];
    ylabelStr = 'b_c / b_p';

    % plot
    sub1 = subplot(2,1,1);
    for Ti = 1:5
        hold on
        plot(xConc,Pgauss(Ti,paramVal2)./Pgauss(Ti,paramVal1),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
    end % plot for each temperature; as black symbol defined from symbol list

    % dress figure
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
    xlim(axisRange)
    xticks(xTicks)

    clear paramVal axisRange ylabelStr
            %% polymerisation - crystallisation
    figure(7)
    paramVal1 = 2;
    paramVal2 = 6;
    axisRange = [xrange];
    ylabelStr = 'b_c - b_p';

    % plot
    sub1 = subplot(2,1,1);
    for Ti = 1:5
        hold on
        plot(xConc,Pgauss(Ti,paramVal2)-Pgauss(Ti,paramVal1),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
    end % plot for each temperature; as black symbol defined from symbol list

    % dress figure
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
    xlim(axisRange)
    xticks(xTicks)

    clear paramVal axisRange ylabelStr 
        %% a plot
    figure(2)
    xlabelStr = 'C1 Concentration (mol%)';    
            %% polymerisation  
            paramVal = 1;
            axisRange = [xrange];
            ylabelStr = 'a_p (W/g)';

            % plot
            sub1 = subplot(2,2,1);
            for Ti = 1:5
                hold on
                plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
            end % plot for each temperature; as black symbol defined from symbol list

            % dress figure
            xlabel(xlabelStr)
            ylabel(ylabelStr)
            sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
            xlim(axisRange)
            xticks(xTicks)

            clear paramVal axisRange ylabelStr           
            %% crystallisation 
    % axes(handles.axes10) %old handle
    paramVal = 5;
    ybound = 30;
    axisRange = [xrange];
    ylabelStr = 'a_c (W/g)';

    % plot
    sub3 = subplot(2,2,3);
    % title('Gaussian b - crystallisation')
    % plot with except for out of range values
    for Ti = 1:5
        hold on

        bc = Pgauss(Ti,paramVal);
        bcOutOfBounds = bc>ybound;
        bc(bc>=ybound) = ybound;

        plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        if bcOutOfBounds == 1 %plot a + as well
        plot(xConc,bc,'k+')    
        end
    end

    % dress
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub3.XMinorGrid = 'on'; sub3.YMinorGrid = 'on'; grid on; 
    xlim(axisRange)
    xticks(xTicks)            
        %% c1 plot
    figure(3)
    xlabelStr = 'C1 Concentration (mol%)';    
            %% polymerisation  
            paramVal = 3;
            axisRange = [xrange];
            ylabelStr = 'c1_p (min)';

            % plot
            sub1 = subplot(2,2,1);
            for Ti = 1:5
                hold on
                plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
            end % plot for each temperature; as black symbol defined from symbol list

            % dress figure
            xlabel(xlabelStr)
            ylabel(ylabelStr)
            sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
            xlim(axisRange)
            xticks(xTicks)

            clear paramVal axisRange ylabelStr           
            %% crystallisation 
    % axes(handles.axes10) %old handle
    paramVal = 7;
    ybound = 30;
    axisRange = [xrange];
    ylabelStr = 'c1_c (min)';

    % plot
    sub3 = subplot(2,2,3);
    % title('Gaussian b - crystallisation')
    % plot with except for out of range values
    for Ti = 1:5
        hold on

        bc = Pgauss(Ti,paramVal);
        bcOutOfBounds = bc>ybound;
        bc(bc>=ybound) = ybound;

        plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        if bcOutOfBounds == 1 %plot a + as well
        plot(xConc,bc,'k+')    
        end
    end

    % dress
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub3.XMinorGrid = 'on'; sub3.YMinorGrid = 'on'; grid on; 
    xlim(axisRange)
    xticks(xTicks)                  
        %% c2 plot
    figure(4)
    xlabelStr = 'C1 Concentration (mol%)';    
            %% polymerisation  
            paramVal = 4;
            axisRange = [xrange];
            ylabelStr = 'c2_p (min)';

            % plot
            sub1 = subplot(2,2,1);
            for Ti = 1:5
                hold on
                plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
            end % plot for each temperature; as black symbol defined from symbol list

            % dress figure
            xlabel(xlabelStr)
            ylabel(ylabelStr)
            sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
            xlim(axisRange)
            xticks(xTicks)

            clear paramVal axisRange ylabelStr           
            %% crystallisation 
    % axes(handles.axes10) %old handle
    paramVal = 8;
    ybound = 30;
    axisRange = [xrange];
    ylabelStr = 'c2_c (min)';

    % plot
    sub3 = subplot(2,2,3);
    % title('Gaussian b - crystallisation')
    % plot with except for out of range values
    for Ti = 1:5
        hold on

        bc = Pgauss(Ti,paramVal);
        bcOutOfBounds = bc>ybound;
        bc(bc>=ybound) = ybound;

        plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        if bcOutOfBounds == 1 %plot a + as well
        plot(xConc,bc,'k+')    
        end
    end

    % dress
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub3.XMinorGrid = 'on'; sub3.YMinorGrid = 'on'; grid on; 
    xlim(axisRange)
    xticks(xTicks)                  
        %% Enthalpies 
        figure(8)
        
        paramVal1 = 10;
        paramVal2 = 11;
        axisRange = [xrange];
        ylabelStr1 = 'H_m (J/g)';
        ylabelStr2 = 'H_p (J/g)';

        % plot
        sub1 = subplot(2,2,1);
        for Ti = 1:5
            hold on
            plot(xConc,Pgauss(Ti,paramVal1),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        end % plot for each temperature; as black symbol defined from symbol list

        % dress figure
        xlabel(xlabelStr)
        ylabel(ylabelStr1)
        sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
        xlim(axisRange)
        xticks(xTicks)

            % plot
        sub2 = subplot(2,2,2);
        for Ti = 1:5
            hold on
            plot(xConc,Pgauss(Ti,paramVal2)-Pgauss(Ti,paramVal1),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        end % plot for each temperature; as black symbol defined from symbol list

        % dress figure
        xlabel(xlabelStr)
        ylabel(ylabelStr2)
        sub2.XMinorGrid = 'on'; sub2.YMinorGrid = 'on'; grid on
        xlim(axisRange)
        xticks(xTicks)

        clear paramVal axisRange ylabelStr 
        %% Parameter Summariser
for Ti = 1:5
    for paramIDX = 1:11
    pGtemp = Pgauss(Ti,paramIDX);

    % Creates a summary for bc
    SummaryArray(foldIDX,Ti,paramIDX) = pGtemp ;
    end
end
end
%% For Second Set
elseif val == 2
% C20P Loop (C1 constant)
for foldIDX = 6:10
    %% Load Data for Specific Concentration Set
Pgauss = Parameters{foldIDX}.gParam;
Pkamal = Parameters{foldIDX}.kParam(fitConstrain,:);
xConc = foldConc(foldIDX,2);
xTicks = foldConc(6:10,2);

    %% Save and Export Table
table = [table; Pkamal];
save('table.mat', 'table')
    %% Plot Kamal Parameters    
% A1
axes(handles.axes1)
hold on
title('A1')
plot(xConc,Pkamal(1),'ro')
A1 = Pkamal(1);

% A2
axes(handles.axes2)
hold on
title('A2')
plot(xConc,Pkamal(3),'ro')
A2 = Pkamal(3);

% m
axes(handles.axes3)
hold on
title('m')
plot(xConc,Pkamal(5),'ro')

% Ea1
axes(handles.axes4)
hold on
title('Ea1')
plot(xConc,Pkamal(2),'ro')
Ea1 = Pkamal(2);

% Ea2
axes(handles.axes5)
hold on
title('Ea2')
plot(xConc,Pkamal(4),'ro')
Ea2 = Pkamal(4);

% n
axes(handles.axes6)
hold on
title('n')
plot(xConc,Pkamal(6),'ro')

% k1
axes(handles.axes7)
hold on
title('k1')
k1 = A1.*exp(-Ea1./R./(150+273));
plot(xConc,k1,'ro')

% k2
axes(handles.axes8)
hold on
title('k2')
k2 = A2.*exp(-Ea2./R./(150+273));
plot(xConc,k2,'ro') 
% safd
    %% Plot Gaussian Parameters  
        %% Initialise (Load Parameters) 
    % Initialise colour and symbols for plotting
    %                 130     140      150      160       170
    colorTriplet = [[1 0 0]; [0 1 0]; [0 0 1]; [1 1 0]; [1 0 1]];
    shapeSet =       {'o';     'v';     's';     '^';     'd'};
    %                 Red    Green     Blue     Yellow   Magenta

    xrange = [0.3 1.3];
    axes(handles.axes9) % old handle 
        %% b plot              
    figure(1)
    xlabelStr = 'C20P Concentration (mol%)';
            %% polymerisation  
    paramVal = 2;
    axisRange = [xrange 0 15];
    ylabelStr = 'b_p (min)';

    % plot
    sub1 = subplot(2,2,2);
    for Ti = 1:5
        hold on
        plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
    end % plot for each temperature; as black symbol defined from symbol list

    % dress figure
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
    axis(axisRange)
    xticks(xTicks)

    clear paramVal axisRange ylabelStr
            %% crystallisation 
    % axes(handles.axes10) %old handle
    paramVal = 6;
    ybound = 30;
    axisRange = [xrange 0 ybound];
    ylabelStr = 'b_c (min)';

    % plot
    sub3 = subplot(2,2,4);
    % title('Gaussian b - crystallisation')
    % plot with except for out of range values
    for Ti = 1:5
        hold on

        bc = Pgauss(Ti,paramVal);
        bcOutOfBounds = bc>ybound;
        bc(bc>=ybound) = ybound;

        plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        if bcOutOfBounds == 1 %plot a + as well
        plot(xConc,bc,'k+')    
        end
    end

    % dress
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub3.XMinorGrid = 'on'; sub3.YMinorGrid = 'on'; grid on; 
    axis(axisRange)
    xticks(xTicks)
            %% polymerisation vs Crystallisation
    figure(5)
    paramVal1 = 2;
    paramVal2 = 6;
%     axisRange = [xrange 0 15];
    ylabelStr = 'b_c (min)';

    % plot
%     sub1 = subplot(2,2,1);
    for Ti = 1:5
        hold on
        plot(Pgauss(Ti,paramVal1),Pgauss(Ti,paramVal2),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
    end % plot for each temperature; as black symbol defined from symbol list

    % dress figure
    xlabel('b_p (min)')
    ylabel(ylabelStr)
    sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
%     axis(axisRange)
%     xticks(xTicks)
            %% polymerisation/crystallisation
    figure(6)
    paramVal1 = 2;
    paramVal2 = 6;
    axisRange = [xrange];
    ylabelStr = 'b_c / b_p';

    % plot
    sub1 = subplot(2,1,2);
    for Ti = 1:5
        hold on
        plot(xConc,Pgauss(Ti,paramVal2)./Pgauss(Ti,paramVal1),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
    end % plot for each temperature; as black symbol defined from symbol list

    % dress figure
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
    xlim(axisRange)
    xticks(xTicks)

    clear paramVal axisRange ylabelStr 
            %% polymerisation - crystallisation
    figure(7)
    paramVal1 = 2;
    paramVal2 = 6;
    axisRange = [xrange];
    ylabelStr = 'b_c - b_p';

    % plot
    sub1 = subplot(2,1,2);
    for Ti = 1:5
        hold on
        plot(xConc,Pgauss(Ti,paramVal2)-Pgauss(Ti,paramVal1),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
    end % plot for each temperature; as black symbol defined from symbol list

    % dress figure
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
    xlim(axisRange)
    xticks(xTicks)

    clear paramVal axisRange ylabelStr 
        %% a plot
    figure(2)
    xlabelStr = 'C20P Concentration (mol%)';    
            %% polymerisation  
            paramVal = 1;
            axisRange = [xrange];
            ylabelStr = 'a_p (W/g)';

            % plot
            sub1 = subplot(2,2,2);
            for Ti = 1:5
                hold on
                plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
            end % plot for each temperature; as black symbol defined from symbol list

            % dress figure
            xlabel(xlabelStr)
            ylabel(ylabelStr)
            sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
            xlim(axisRange)
            xticks(xTicks)

            clear paramVal axisRange ylabelStr           
            %% crystallisation 
    % axes(handles.axes10) %old handle
    paramVal = 5;
    ybound = 30;
    axisRange = [xrange];
    ylabelStr = 'a_c (W/g)';

    % plot
    sub3 = subplot(2,2,4);
    % title('Gaussian b - crystallisation')
    % plot with except for out of range values
    for Ti = 1:5
        hold on

        bc = Pgauss(Ti,paramVal);
        bcOutOfBounds = bc>ybound;
        bc(bc>=ybound) = ybound;

        plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        if bcOutOfBounds == 1 %plot a + as well
        plot(xConc,bc,'k+')    
        end
    end

    % dress
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub3.XMinorGrid = 'on'; sub3.YMinorGrid = 'on'; grid on; 
    xlim(axisRange)
    xticks(xTicks)            
        %% c1 plot
    figure(3)
    xlabelStr = 'C20P Concentration (mol%)';    
            %% polymerisation  
            paramVal = 3;
            axisRange = [xrange];
            ylabelStr = 'c1_p (min)';

            % plot
            sub1 = subplot(2,2,2);
            for Ti = 1:5
                hold on
                plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
            end % plot for each temperature; as black symbol defined from symbol list

            % dress figure
            xlabel(xlabelStr)
            ylabel(ylabelStr)
            sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
            xlim(axisRange)
            xticks(xTicks)

            clear paramVal axisRange ylabelStr           
            %% crystallisation 
    % axes(handles.axes10) %old handle
    paramVal = 7;
    ybound = 30;
    axisRange = [xrange];
    ylabelStr = 'c1_c (min)';

    % plot
    sub3 = subplot(2,2,4);
    % title('Gaussian b - crystallisation')
    % plot with except for out of range values
    for Ti = 1:5
        hold on

        bc = Pgauss(Ti,paramVal);
        bcOutOfBounds = bc>ybound;
        bc(bc>=ybound) = ybound;

        plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        if bcOutOfBounds == 1 %plot a + as well
        plot(xConc,bc,'k+')    
        end
    end

    % dress
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub3.XMinorGrid = 'on'; sub3.YMinorGrid = 'on'; grid on; 
    xlim(axisRange)
    xticks(xTicks)                  
        %% c2 plot
    figure(4)
    xlabelStr = 'C20P Concentration (mol%)';    
            %% polymerisation  
            paramVal = 4;
            axisRange = [xrange];
            ylabelStr = 'c2_p (min)';

            % plot
            sub1 = subplot(2,2,2);
            for Ti = 1:5
                hold on
                plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
            end % plot for each temperature; as black symbol defined from symbol list

            % dress figure
            xlabel(xlabelStr)
            ylabel(ylabelStr)
            sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
            xlim(axisRange)
            xticks(xTicks)

            clear paramVal axisRange ylabelStr           
            %% crystallisation 
    % axes(handles.axes10) %old handle
    paramVal = 8;
    ybound = 30;
    axisRange = [xrange];
    ylabelStr = 'c2_c (min)';

    % plot
    sub3 = subplot(2,2,4);
    % title('Gaussian b - crystallisation')
    % plot with except for out of range values
    for Ti = 1:5
        hold on

        bc = Pgauss(Ti,paramVal);
        bcOutOfBounds = bc>ybound;
        bc(bc>=ybound) = ybound;

        plot(xConc,Pgauss(Ti,paramVal),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        if bcOutOfBounds == 1 %plot a + as well
        plot(xConc,bc,'k+')    
        end
    end

    % dress
    xlabel(xlabelStr)
    ylabel(ylabelStr)
    sub3.XMinorGrid = 'on'; sub3.YMinorGrid = 'on'; grid on; 
    xlim(axisRange)
    xticks(xTicks)                  
        %% Enthalpies 
        figure(8)
        
        paramVal1 = 10;
        paramVal2 = 11;
        axisRange = [xrange];
        ylabelStr1 = 'H_m (J/g)';
        ylabelStr2 = 'H_p (J/g)';

        % plot
        sub1 = subplot(2,2,3);
        for Ti = 1:5
            hold on
            plot(xConc,Pgauss(Ti,paramVal1),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        end % plot for each temperature; as black symbol defined from symbol list

        % dress figure
        xlabel(xlabelStr)
        ylabel(ylabelStr1)
        sub1.XMinorGrid = 'on'; sub1.YMinorGrid = 'on'; grid on
        xlim(axisRange)
        xticks(xTicks)

            % plot
        sub2 = subplot(2,2,4);
        for Ti = 1:5
            hold on
            plot(xConc,Pgauss(Ti,paramVal2)-Pgauss(Ti,paramVal1),'Marker',shapeSet{Ti},'MarkerEdgeColor','k')
        end % plot for each temperature; as black symbol defined from symbol list

        % dress figure
        xlabel(xlabelStr)
        ylabel(ylabelStr2)
        sub2.XMinorGrid = 'on'; sub2.YMinorGrid = 'on'; grid on
        xlim(axisRange)
        xticks(xTicks)

        clear paramVal axisRange ylabelStr 
        %% Parameter Summariser
for Ti = 1:5
    for paramIDX = 1:11
    pGtemp = Pgauss(Ti,paramIDX);

    % Creates a summary for bc
    SummaryArray(foldIDX,Ti,paramIDX) = pGtemp ;
    end
end
end
end

    function popupmenu1_Callback(hObject, eventdata, handles)

    cla(handles.axes1); cla(handles.axes2); cla(handles.axes3);
    cla(handles.axes4); cla(handles.axes5); cla(handles.axes6);
    cla(handles.axes7); cla(handles.axes8); cla(handles.axes9);
    cla(handles.axes10);

    UpdatePlots(hObject, eventdata, handles);
    % hObject    handle to popupmenu1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from popupmenu1


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


    % --- Executes on button press in pushbutton1.
    function pushbutton1_Callback(hObject, eventdata, handles)
    popupmenu1_Callback(hObject, eventdata, handles)

