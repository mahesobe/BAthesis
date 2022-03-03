function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 25-Mar-2015 13:05:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in exitbutton.
function exitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to exitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection=0;
setappdata(handles.figure1, 'selection', selection);
uiresume(handles.figure1)



% --- Executes on button press in step7.
function step7_Callback(hObject, eventdata, handles)
% hObject    handle to step7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set(handles.figure1,'Visible','off')
drawnow
setappdata(handles.figure1, 'selection', 7);
set(handles.step7,'string','running','ForegroundColor','red','enable','off');
drawnow
uiresume(handles.figure1)


% --- Executes on button press in step6.
function step6_Callback(hObject, eventdata, handles)
% hObject    handle to step6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1, 'selection', 6);
set(handles.step6,'string','running','ForegroundColor','red','enable','off');
drawnow
uiresume(handles.figure1)


% --- Executes on button press in step5.
function step5_Callback(hObject, eventdata, handles)
% hObject    handle to step5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1, 'selection', 5);
set(handles.step5,'string','running','ForegroundColor','red','enable','off');
drawnow
uiresume(handles.figure1)


% --- Executes on button press in step4.
function step4_Callback(hObject, eventdata, handles)
% hObject    handle to step4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1)
setappdata(handles.figure1, 'selection', 4);
set(handles.step4,'string','running','ForegroundColor','red','enable','off');
drawnow


% --- Executes on button press in step2.
function step3_Callback(hObject, eventdata, handles)
% hObject    handle to step2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1, 'selection', 3);
set(handles.step3,'string','running','ForegroundColor','red','enable','off');
drawnow
uiresume(handles.figure1)


% --- Executes on button press in step2.
function step2_Callback(hObject, eventdata, handles)
% hObject    handle to step2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figure1, 'selection', 2);
set(handles.step2,'string','running','ForegroundColor','red','enable','off');
drawnow
uiresume(handles.figure1)


% --- Executes on button press in step1.
function step1_Callback(hObject, eventdata, handles)
% hObject    handle to step1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1)
setappdata(handles.figure1, 'selection', 1);
set(handles.step1,'string','running','ForegroundColor','red','enable','off');
drawnow


