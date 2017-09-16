function varargout = IMAGE2LOC(varargin)
% IMAGE2LOC MATLAB code for IMAGE2LOC.fig
%      IMAGE2LOC, by itself, creates a new IMAGE2LOC or raises the existing
%      singleton*.
%
%      H = IMAGE2LOC returns the handle to a new IMAGE2LOC or the handle to
%      the existing singleton*.
%
%      IMAGE2LOC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE2LOC.M with the given input arguments.
%
%      IMAGE2LOC('Property','Value',...) creates a new IMAGE2LOC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IMAGE2LOC_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IMAGE2LOC_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IMAGE2LOC

% Last Modified by GUIDE v2.5 27-May-2017 16:07:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IMAGE2LOC_OpeningFcn, ...
                   'gui_OutputFcn',  @IMAGE2LOC_OutputFcn, ...
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


% --- Executes just before IMAGE2LOC is made visible.
function IMAGE2LOC_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IMAGE2LOC (see VARARGIN)

% Choose default command line output for IMAGE2LOC
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IMAGE2LOC wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IMAGE2LOC_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in BUTTON_ADDPROBE.
function BUTTON_ADDPROBE_Callback(hObject, eventdata, handles)
% hObject    handle to BUTTON_ADDPROBE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function SLIDER_THRESHOLD_Callback(hObject, eventdata, handles)
% hObject    handle to SLIDER_THRESHOLD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function SLIDER_THRESHOLD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SLIDER_THRESHOLD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in BUTTON_LOCALISE.
function BUTTON_LOCALISE_Callback(hObject, eventdata, handles)
% hObject    handle to BUTTON_LOCALISE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in MENU_NNP.
function MENU_NNP_Callback(hObject, eventdata, handles)
% hObject    handle to MENU_NNP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MENU_NNP contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MENU_NNP


% --- Executes during object creation, after setting all properties.
function MENU_NNP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MENU_NNP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in BUTTON_SAVELOC.
function BUTTON_SAVELOC_Callback(hObject, eventdata, handles)
% column heading
%image-ID      frame number from total acquisition
%time-point
%cycle         cycle is not relevant here, use 0 throughout
%z-step        for z stacks, maybe not useful here as we use x,y,z coord
%frame         frame number within a cycle, so identical to image-ID here
%accum         not sure
%probe         probe id
%photon-count    total photon count  (use intensity as photon-count)
%photon-count11  total photon countof probe 1 in channel 1
%photon-count12  total photon countof probe 1 in channel 2
%photon-count21  total photon countof probe 2 in channel 1
%photon-count22  total photon countof probe 2 in channel 2
%psfx            psfx error
%psfy            psfy error
%psfz            psfz error
%psf-photon-count  not sure
%x               x localised position
%y               y localised position
%z               z localised position
%stdev          not sure
%amp            not sure
%background11   background photon count for probe 1 channel 1
%background12   background photon count for probe 1 channel 2
%background21   background photon count for probe 2 channel 1
%background22   background photon count for probe 2 channel 2
%maxResidualSlope
%chisq          localisation quality metric
%log-likelihood localisation quality metric
%llr
%accuracy
%fiducial       fiducial marker id, not used in this case
%valid          localisation validity, all valid in this case


function localfit
    [X,Y] = meshgrid(-MdataSize/2:MdataSize/2);
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;
xdata(1,1,:)
lb = [0,-MdataSize/2,0,-MdataSize/2,0];
    ub = [realmax('double'),MdataSize/2,(MdataSize/2)^2,MdataSize/2,(MdataSize/2)^2];
    [x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunction,x0,xdata,Z,lb,ub);
    
   function F = D2GaussFunction(x,xdata)
 F = x(1)*exp(-((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + (xdata(:,:,2)-x(4)).^2/(2*x(5)^2)));