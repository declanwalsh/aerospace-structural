function varargout = panelGUI(varargin)

T_SKIN_MIN = 0.1;
T_SKIN_DEFAULT = 1;
T_SKIN_MAX = 3;


PITCH_STIFF_MIN = 1;
PITCH_STIFF_DEFAULT = 20;
PITCH_STIFF_MAX = 100;

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @panelGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @panelGUI_OutputFcn, ...
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

%  Create and then hide the UI as it is being constructed.
S.f = figure('Visible','off','Position',[600, 0,400,400]);

slider_t_skin = uicontrol('Style', 'slider', 'min', T_SKIN_MIN, 'max' , T_SKIN_MAX, ...
    'value', T_SKIN_DEFAULT, 'Position', [200, 20, 120, 20], 'Callback', {@panel, S, 't_skin'});
txt_t_skin = uicontrol('Style', 'text', 'Position', [200, 45, 120, 20], ...
    'String', 'Skin Thickness (mm)');

slider_pitch_stiffener = uicontrol('Style', 'slider', 'min', PITCH_STIFF_MIN, 'max' , PITCH_STIFF_MAX, ...
    'value', PITCH_STIFF_DEFAULT, 'Position', [0, 20, 120, 20], 'Callback', {@panel, S, 'pitch_stiffener'});
txt_pitch_stiffener = uicontrol('Style', 'text', 'Position', [0, 45, 120, 20], ...
    'String', 'Stiffener Pitch (mm)');

S.f.Visible = 'on';

function panelGUI_OpeningFcn(hObject, eventdata, handles)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mygui (see VARARGIN)
 
% Choose default command line output for mygui
handles.output = hObject;

handles.t_skin = 1;
handles.pitch_stiffener = 20;
 
% Update handles structure
guidata(hObject, handles);
