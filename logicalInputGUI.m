function varargout = logicalInputGUI(varargin)
% LOGICALINPUTGUI MATLAB code for logicalInputGUI.fig
%      LOGICALINPUTGUI, by itself, creates a new LOGICALINPUTGUI or raises the existing
%      singleton*.
%
%      H = LOGICALINPUTGUI returns the handle to a new LOGICALINPUTGUI or the handle to
%      the existing singleton*.
%
%      LOGICALINPUTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOGICALINPUTGUI.M with the given input arguments.
%
%      LOGICALINPUTGUI('Property','Value',...) creates a new LOGICALINPUTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before logicalInputGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to logicalInputGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help logicalInputGUI

% Last Modified by GUIDE v2.5 13-Mar-2019 00:51:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @logicalInputGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @logicalInputGUI_OutputFcn, ...
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


% --- Executes just before logicalInputGUI is made visible.
function logicalInputGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to logicalInputGUI (see VARARGIN)

load('dataExclusion.mat');
Exc = double(Exc);

for i = 1:10
   tableui = sprintf('uitable%d',i);
   uit = handles.(tableui);
   A = squeeze(double(Exc(i,:,:)));
   uit.Data = logical(A');
   uit.ColumnEditable = true;
   uit.ColumnName = {'130 C', '140 C', '150 C', '160 C', '170 C'};
end


% Choose default command line output for logicalInputGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes logicalInputGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = logicalInputGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i = 1:10
   tableui = sprintf('uitable%d',i);
   ExcFoldi = get(handles.(tableui),'Data');
   Exc(i,:,:) = ExcFoldi';
end

save('dataExclusion.mat','Exc')
