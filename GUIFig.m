function varargout = GUIFig(varargin)
% GUIFig MATLAB code for GUIFig.fig
%      GUIFig, by itself, creates a new GUIFig or raises the existing
%      singleton*.
%
%      H = GUIFig returns the handle to a new GUIFig or the handle to
%      the existing singleton*.
%
%      GUIFig('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIFig.M with the given input arguments.
%
%      GUIFig('Property','Value',...) creates a new GUIFig or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIFig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIFig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIFig

% Last Modified by GUIDE v2.5 20-Mar-2019 20:30:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIFig_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIFig_OutputFcn, ...
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


% --- Executes just before GUIFig is made visible.
function GUIFig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIFig (see VARARGIN)
PlotData = load('PlotData.mat');
PlotData = PlotData.PlotData;

foldIDX = 1;
fileIDX = 1;
    
uit = uitable(handles.uitable1);
RowLabels = {'T';'Rn';'Hmelt';'Htot';'Hpoly'};

uit.Data = [RowLabels];

% plotFunc(foldIDX,fileIDX,handles)
checkboxLabel(foldIDX,handles)
pushbutton4_Callback(hObject, eventdata, handles)

% Choose default command line output for GUIFig
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIFig wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIFig_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Determine the selected data set
str = get(handles.popupmenu4,'String');
val = get(handles.popupmenu4,'Value');
%Set current data to the selected data set

foldIDX = val;

%clear figures and tables
cla(handles.axes1); cla(handles.axes2); cla(handles.axes3); cla(handles.axes4);cla(handles.axes5);

uit = uitable(handles.uitable1);
uit.Data = {'T';'Rn';'Hmelt';'Htot';'Hpoly'};

checkboxLabel(foldIDX,handles)
foldConc = [0.6, 0.6; 0.8, 0.6; 1.0, 0.6; 1.2, 0.6; 1.4, 0.6; 1.2, 0.4; 1.2, 0.6; 1.2, 0.8; 1.2, 1.0; 1.2, 1.2];
concString = sprintf('%2.2f%% C1 & %2.2f%% C20P',foldConc(foldIDX,1),foldConc(foldIDX,2));
set(handles.text4,'String',concString)

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = get(handles.popupmenu4,'String');
val = get(handles.popupmenu4,'Value');
%Set current data to the selected data set

%clear figures and tables
cla(handles.axes1); cla(handles.axes2); cla(handles.axes3); cla(handles.axes4);cla(handles.axes5);

uit = uitable(handles.uitable1);
uit.Data = {'T';'Rn';'Hmelt';'Htot';'Hpoly'};

%set fold IDX
foldIDX = val;
fileIDX = 1:12; %will need to get this from somewhere...
CheckBoxCheck = [];
    for i = fileIDX
       ElementHandle = sprintf('checkbox%d',i);
       CheckBoxCheck = [CheckBoxCheck get(handles.(ElementHandle),'Value')];
    end
    
fileIDX = fileIDX.*CheckBoxCheck;

fileIDX(fileIDX ==0) = [];
Counter = 0;
    
    for fileIDX = fileIDX
    plotFunc(foldIDX,fileIDX,Counter,handles)
    Counter = Counter+1;
    end
figGUI = gcf;

%Save Figure
figure(foldIDX)
figName = 'Gaussian Separation';
xlabel('time (minutes)')
ylabel('Heat Flow W/g')
legend('Heat Flow Data (W/g)','Reduced Polymerisation Gaussian Curve','Reduced Crystallisation Gaussian Curve','Total Modelled Heat Flow','Location','NorthWest')
grid on
grid minor

foldIDXstr = num2str(foldIDX);
fileSaveStr = strcat('figures/',figName,'_Folder',foldIDXstr,'.png'); % can change preferred image format here (i.e. epsc, fig, png)
saveas(gcf,fileSaveStr)
% cla(foldIDX)
   
figure(figGUI)
guidata(hObject,handles)

% --- Updates Plots and Tables
function plotFunc(foldIDX,fileIDX,Counter,handles)
    PlotData = load('PlotData.mat');
    PlotData = PlotData.PlotData;

%Raw
    axes(handles.axes1)
    hold on
    handles.Raw = PlotData{foldIDX,fileIDX}.Raw;
    plot(handles.Raw(:,1),handles.Raw(:,2))
    title('Raw DSC Data')

%Melting & Baseline
    axes(handles.axes2)
    hold on
    handles.hFlowTrace = PlotData{foldIDX,fileIDX}.hFlowTrace;
    handles.hFlowBaseline = PlotData{foldIDX,fileIDX}.hFlowBaseline;

    plot(handles.hFlowTrace(:,1),handles.hFlowTrace(:,2));
    hold on; plot(handles.hFlowBaseline(:,1),handles.hFlowBaseline(:,2),'--');
    title('Melting and Baseline')

%Polymerisation & Baseline
    axes(handles.axes3)
    hold on
    handles.hFlowRxn = PlotData{foldIDX,fileIDX}.hFlowRxn;
    handles.hFlowRxnBL = PlotData{foldIDX,fileIDX}.hFlowRxnBL;
    grid minor
    grid on
    zoom on

    plot(handles.hFlowRxn(:,1),handles.hFlowRxn(:,2));
    plot(handles.hFlowRxnBL(:,1),handles.hFlowRxnBL(:,2),'--');
    title('Polymerisation and Baseline')

%Reduced Polymerisation Curve
    axes(handles.axes4)
    hold on
    handles.hFRxnRed = PlotData{foldIDX,fileIDX}.hFRxnRed;

    plot(handles.hFRxnRed(:,1),handles.hFRxnRed(:,2))
    title('Reduced Polymerisation Curve')

%Gauss Curve
    axes(handles.axes5)
    hold on
    handles.gaussCurve = PlotData{foldIDX,fileIDX}.gaussCurve;

    plot(handles.gaussCurve(:,1),handles.gaussCurve(:,2));
    plot(handles.gaussCurve(:,1),handles.gaussCurve(:,3),':');
    plot(handles.gaussCurve(:,1),handles.gaussCurve(:,4),'--');
    title('Gaussian Model Fitting')
    
    
%Figure
    figure(foldIDX)
    hold on
    plot(handles.hFRxnRed(:,1)+25*Counter,handles.hFRxnRed(:,2)+0.15*Counter,'k')
    idx = length(handles.hFRxnRed(:,1)) - length(handles.gaussCurve(:,1));
    
    plot(handles.hFRxnRed(idx+1:end,1)+25*Counter,handles.gaussCurve(:,2)+0.15*Counter,'b-');
    plot(handles.hFRxnRed(idx+1:end,1)+25*Counter,handles.gaussCurve(:,3)+0.15*Counter,'r-');
    plot(handles.hFRxnRed(idx+1:end,1)+25*Counter,handles.gaussCurve(:,4)+0.15*Counter,'k--');
    
% 
    figure(21)
    hold on
    plot(handles.hFRxnRed(idx+1:end,1),handles.gaussCurve(:,4)-handles.hFRxnRed(idx+1:end,2),'.')
       axes(handles.axes5)
       
    figure(22)
    hold on
    plot(handles.gaussCurve(:,4)./max(handles.gaussCurve(:,4)),handles.gaussCurve(:,4)-handles.hFRxnRed(idx+1:end,2),'.')
       axes(handles.axes5)
       
       
%Table Data
    Ti = PlotData{foldIDX,fileIDX}.Ti;
    Rn = PlotData{foldIDX,fileIDX}.Rn;
    Hmelt = PlotData{foldIDX,fileIDX}.Hmelt;
    Htotal = PlotData{foldIDX,fileIDX}.Htotal; %PlotData... 250 is placeholder
    Hpoly = Htotal - Hmelt;

    uit = uitable(handles.uitable1);

    Data = get(handles.uitable1,'Data');
    uit.Data = [Data {Ti;Rn;Hmelt;Htotal;Hpoly}];
     
% --- Renames checkboxes given folder
function checkboxLabel(foldIDX,handles)
    PlotData = load('PlotData.mat');
    PlotData = PlotData.PlotData;
%     handles
    for i = 1:size(PlotData,2)
        if isempty(PlotData{foldIDX,i})
       ElementHandle = sprintf('checkbox%d',i);
       TempRunString = 'No Data';
       set(handles.(ElementHandle),'String',TempRunString,'Value',0,'Visible','Off')       
        else
       Ti = PlotData{foldIDX,i}.Ti;
       Rn = PlotData{foldIDX,i}.Rn;
       ElementHandle = sprintf('checkbox%d',i);
       TempRunString = sprintf('%d - %d',Ti,Rn);
       set(handles.(ElementHandle),'String',TempRunString,'Value',1,'Visible','On')
       end
    end    

% --- Executes on button press in checkboxCLEAR.
function checkbox1_Callback(hObject, eventdata, handles)
function checkbox2_Callback(hObject, eventdata, handles)
function checkbox3_Callback(hObject, eventdata, handles)
function checkbox4_Callback(hObject, eventdata, handles)
function checkbox5_Callback(hObject, eventdata, handles)
function checkbox6_Callback(hObject, eventdata, handles)
function checkbox7_Callback(hObject, eventdata, handles)
function checkbox8_Callback(hObject, eventdata, handles)
function checkbox9_Callback(hObject, eventdata, handles)
function checkbox10_Callback(hObject, eventdata, handles)
function checkbox11_Callback(hObject, eventdata, handles)
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCLEAR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCLEAR    
    
    
% --- Executes on button press in checkboxCLEAR.
function checkboxCLEAR_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCLEAR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    for i = 1:12
       ElementHandle = sprintf('checkbox%d',i);
       set(handles.(ElementHandle),'Value',0)
    end

% Hint: get(hObject,'Value') returns toggle state of checkboxCLEAR


% --- Executes on button press in checkboxALL.
    function checkboxALL_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxALL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    for i = 1:12
       ElementHandle = sprintf('checkbox%d',i);
       Vis = get(handles.(ElementHandle),'Visible');
       set(handles.(ElementHandle),'Value',1)
    end

% Hint: get(hObject,'Value') returns toggle state of checkboxALL


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run logicalInputGUI


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run KamalFitGUI


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run LoadSeparation

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
    %               legend('130 C','140 C', '150 C', '160 C', '170 C', 'Model')
    %               str = get(handles.popupmenu1,'String');
    foldIDX = get(handles.popupmenu4,'Value');
    foldIDXstr = num2str(foldIDX);
    fileSaveStr = strcat('figures/',figName,'_Folder',foldIDXstr,'.png'); % can change preferred image format here (i.e. epsc, fig, png)
    %               figure('PaperPosition',[.25 .25 8 6]);
    saveas(gcf,fileSaveStr)  
    cla(saveFigHandle)
    axes(axesH)

