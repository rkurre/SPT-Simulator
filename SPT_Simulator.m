function varargout = SPT_Simulator(varargin)
% SPT_SIMULATOR MATLAB code for SPT_Simulator.fig
%      SPT_SIMULATOR, by itself, creates a new SPT_SIMULATOR or raises the existing
%      singleton*.
%
%      H = SPT_SIMULATOR returns the handle to a new SPT_SIMULATOR or the handle to
%      the existing singleton*.
%
%      SPT_SIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPT_SIMULATOR.M with the given input arguments.
%
%      SPT_SIMULATOR('Property','Value',...) creates a new SPT_SIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SPT_Simulator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SPT_Simulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SPT_Simulator

% Last Modified by GUIDE v2.5 16-Apr-2024 10:55:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SPT_Simulator_OpeningFcn, ...
                   'gui_OutputFcn',  @SPT_Simulator_OutputFcn, ...
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


% --- Executes just before SPT_Simulator is made visible.
function SPT_Simulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SPT_Simulator (see VARARGIN)

wArea = str2double(get(handles.wArea, 'String'));
hArea = str2double(get(handles.hArea, 'String'));
im = zeros(hArea, wArea);
imagesc(handles.im, im);
% Choose default command line output for SPT_Simulator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SPT_Simulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SPT_Simulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function dt_Callback(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt as text
%        str2double(get(hObject,'String')) returns contents of dt as a double


% --- Executes during object creation, after setting all properties.
function dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wArea_Callback(hObject, eventdata, handles)
% hObject    handle to wArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wArea as text
%        str2double(get(hObject,'String')) returns contents of wArea as a double


% --- Executes during object creation, after setting all properties.
function wArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hArea_Callback(hObject, eventdata, handles)
% hObject    handle to hArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hArea as text
%        str2double(get(hObject,'String')) returns contents of hArea as a double


% --- Executes during object creation, after setting all properties.
function hArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pxl_sz_Callback(hObject, eventdata, handles)
% hObject    handle to pxl_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pxl_sz as text
%        str2double(get(hObject,'String')) returns contents of pxl_sz as a double


% --- Executes during object creation, after setting all properties.
function pxl_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pxl_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wavelength_Callback(hObject, eventdata, handles)
% hObject    handle to wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wavelength as text
%        str2double(get(hObject,'String')) returns contents of wavelength as a double


% --- Executes during object creation, after setting all properties.
function wavelength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NA_Callback(hObject, eventdata, handles)
% hObject    handle to NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NA as text
%        str2double(get(hObject,'String')) returns contents of NA as a double


% --- Executes during object creation, after setting all properties.
function NA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cycles_Callback(hObject, eventdata, handles)
% hObject    handle to cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cycles as text
%        str2double(get(hObject,'String')) returns contents of cycles as a double


% --- Executes during object creation, after setting all properties.
function cycles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frames_Callback(hObject, eventdata, handles)
% hObject    handle to frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frames as text
%        str2double(get(hObject,'String')) returns contents of frames as a double


% --- Executes during object creation, after setting all properties.
function frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function camNoise_Callback(hObject, eventdata, handles)
% hObject    handle to camNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of camNoise as text
%        str2double(get(hObject,'String')) returns contents of camNoise as a double


% --- Executes during object creation, after setting all properties.
function camNoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ecf_Callback(hObject, eventdata, handles)
% hObject    handle to ecf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ecf as text
%        str2double(get(hObject,'String')) returns contents of ecf as a double


% --- Executes during object creation, after setting all properties.
function ecf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ecf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function QE_Callback(hObject, eventdata, handles)
% hObject    handle to QE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of QE as text
%        str2double(get(hObject,'String')) returns contents of QE as a double


% --- Executes during object creation, after setting all properties.
function QE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to QE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function camOffset_Callback(hObject, eventdata, handles)
% hObject    handle to camOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of camOffset as text
%        str2double(get(hObject,'String')) returns contents of camOffset as a double


% --- Executes during object creation, after setting all properties.
function camOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigE_Callback(hObject, eventdata, handles)
% hObject    handle to sigE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigE as text
%        str2double(get(hObject,'String')) returns contents of sigE as a double


% --- Executes during object creation, after setting all properties.
function sigE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigBG_Callback(hObject, eventdata, handles)
% hObject    handle to sigBG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigBG as text
%        str2double(get(hObject,'String')) returns contents of sigBG as a double


% --- Executes during object creation, after setting all properties.
function sigBG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigBG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pDen_Callback(hObject, eventdata, handles)
% hObject    handle to pDen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pDen as text
%        str2double(get(hObject,'String')) returns contents of pDen as a double


% --- Executes during object creation, after setting all properties.
function pDen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pDen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ltDye_Callback(hObject, eventdata, handles)
% hObject    handle to ltDye (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ltDye as text
%        str2double(get(hObject,'String')) returns contents of ltDye as a double


% --- Executes during object creation, after setting all properties.
function ltDye_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ltDye (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function oligState_Callback(hObject, eventdata, handles)
% hObject    handle to oligState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oligState as text
%        str2double(get(hObject,'String')) returns contents of oligState as a double


% --- Executes during object creation, after setting all properties.
function oligState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oligState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frOlig_Callback(hObject, eventdata, handles)
% hObject    handle to frOlig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frOlig as text
%        str2double(get(hObject,'String')) returns contents of frOlig as a double


% --- Executes during object creation, after setting all properties.
function frOlig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frOlig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dm_Callback(hObject, eventdata, handles)
% hObject    handle to Dm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dm as text
%        str2double(get(hObject,'String')) returns contents of Dm as a double


% --- Executes during object creation, after setting all properties.
function Dm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Do_Callback(hObject, eventdata, handles)
% hObject    handle to Do (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Do as text
%        str2double(get(hObject,'String')) returns contents of Do as a double


% --- Executes during object creation, after setting all properties.
function Do_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Do (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in startSim.
function startSim_Callback(hObject, eventdata, handles)
% hObject    handle to startSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dt = str2double(get(handles.dt, 'String'))/1000;
wAreaPxl = str2double(get(handles.wArea, 'String'));
hAreaPxl = str2double(get(handles.hArea, 'String'));
pxlSzIm = str2double(get(handles.pxl_sz, 'String'))/1000;
nCycles = str2double(get(handles.cycles, 'String')); 
nCycleFrames = str2double(get(handles.frames, 'String'));
lambda = str2double(get(handles.wavelength, 'String'))/1000;
NA = str2double(get(handles.NA, 'String'));
QE = str2double(get(handles.QE, 'String'));
ph2cnts = 1/str2double(get(handles.ecf, 'String'));
camOffset = str2double(get(handles.camOffset, 'String'));
camNoise = str2double(get(handles.camNoise, 'String'));
useGauss = get(handles.useGaussProfile, 'Value');
sigGauss = str2double(get(handles.illGaussSigma, 'String'));
pd = str2double(get(handles.pd, 'String'));
sigE = str2double(get(handles.sigE, 'String'));
sigBG = str2double(get(handles.sigBG, 'String'));
pDen = str2double(get(handles.pDen, 'String'));
ltDye = str2double(get(handles.ltDye, 'String'));
oligState = str2double(get(handles.oligState, 'String'));
frOlig = str2double(get(handles.frOlig, 'String'));
Dm = str2double(get(handles.Dm, 'String'));
Do = str2double(get(handles.Do, 'String'));
zAvg = str2double(get(handles.zAvg, 'String'));
zStd = str2double(get(handles.zStd, 'String'));
nSets = str2double(get(handles.nSets, 'String'));
fileName = get(handles.fileName,'String');
simFrames(dt, wAreaPxl, hAreaPxl, pxlSzIm, nCycles, nCycleFrames,...
    lambda, NA, QE, ph2cnts, camOffset, camNoise, useGauss, sigGauss, ...
    pd, sigE, sigBG, pDen, ltDye, oligState, frOlig, Dm, Do, zAvg, zStd, ...
    nSets, fileName, handles);

function fileName_Callback(hObject, eventdata, handles)
% hObject    handle to fileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileName as text
%        str2double(get(hObject,'String')) returns contents of fileName as a double


% --- Executes during object creation, after setting all properties.
function fileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zAvg_Callback(hObject, eventdata, handles)
% hObject    handle to zAvg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zAvg as text
%        str2double(get(hObject,'String')) returns contents of zAvg as a double


% --- Executes during object creation, after setting all properties.
function zAvg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zAvg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zStd_Callback(hObject, eventdata, handles)
% hObject    handle to zStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zStd as text
%        str2double(get(hObject,'String')) returns contents of zStd as a double


% --- Executes during object creation, after setting all properties.
function zStd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in useGaussProfile.
function useGaussProfile_Callback(hObject, eventdata, handles)
% hObject    handle to useGaussProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useGaussProfile
if get(hObject,'Value')==1
    set(hObject,'Value',1);
    set(handles.useFlattopProfile,'Value',0);
else
    set(hObject,'Value',1);
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in useFlattopProfile.
function useFlattopProfile_Callback(hObject, eventdata, handles)
% hObject    handle to useFlattopProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useFlattopProfile

if get(hObject,'Value')==1
    set(hObject,'Value',1);
    set(handles.useGaussProfile,'Value',0);
else
    set(hObject,'Value',1);
end
% Update handles structure
guidata(hObject, handles);

function pd_Callback(hObject, eventdata, handles)
% hObject    handle to pd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pd as text
%        str2double(get(hObject,'String')) returns contents of pd as a double


% --- Executes during object creation, after setting all properties.
function pd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function illGaussSigma_Callback(hObject, eventdata, handles)
% hObject    handle to illGaussSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of illGaussSigma as text
%        str2double(get(hObject,'String')) returns contents of illGaussSigma as a double


% --- Executes during object creation, after setting all properties.
function illGaussSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to illGaussSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nSets_Callback(hObject, eventdata, handles)
% hObject    handle to nSets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nSets as text
%        str2double(get(hObject,'String')) returns contents of nSets as a double


% --- Executes during object creation, after setting all properties.
function nSets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nSets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
