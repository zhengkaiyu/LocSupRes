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

% Last Modified by GUIDE v2.5 03-Nov-2017 22:24:12

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
function IMAGE2LOC_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IMAGE2LOC (see VARARGIN)

% Choose default command line output for IMAGE2LOC
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = IMAGE2LOC_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
handles.VAL_NNPIXEL.Value=1;handles.VAL_NNPIXEL.String=1;
handles.VAL_ZNNPIXEL.Value=1;handles.VAL_ZNNPIXEL.String=1;
handles.VAL_THRESHOLD.Value=30;handles.VAL_THRESHOLD.String=30;
%=========================================================================
% --- Executes on button press in BUTTON_LOADIMG.
function BUTTON_LOADIMG_Callback(~, ~, handles) %#ok<*DEFNU>
global tiffobj rawdata localdata;
try
    % load raw tiff image
    [pathname,~,~]=fileparts(handles.FIGURE_IMG2LOC.Name);
    if isempty(pathname)
        pathname='./';
    end
    % ask for one file to open
    [filename,pathname,~] = uigetfile({'*.tif;*.TIF;*.tiff;*.TIFF','TIFF image file (*.tif,*.TIF,*.tiff,*.TIFF)';...
        '*.*','All Files (*.*)'},...
        'Select Saved Localisation Analysis File',...
        'MultiSelect','off',pathname);
    if ischar(filename)
        filename=cat(2,pathname,filename);
        tiffobj = Tiff(filename,'r');
        xdim=tiffobj.getTag('ImageLength');
        handles.VAL_IMGLENGTH.String=xdim;
        ydim=tiffobj.getTag('ImageWidth');
        handles.VAL_IMGWIDTH.String=ydim;
        bitdim=tiffobj.getTag('BitsPerSample');
        handles.VAL_IMGBITDEPTH.String=bitdim;
        val=tiffobj.getTag('XResolution');
        handles.VAL_IMGXRES.String=val;
        handles.VAL_RESX.Value=1/val;
        rawdata.dx=1/val;
        rawdata.x=linspace(0,(xdim-1)*rawdata.dx,xdim);
        handles.VAL_RESX.String=1/val;
        val=tiffobj.getTag('YResolution');
        handles.VAL_IMGYRES.String=val;
        handles.VAL_RESY.Value=1/val;
        rawdata.dy=1/val;
        rawdata.y=linspace(0,(ydim-1)*rawdata.dy,ydim);
        handles.VAL_RESY.String=1/val;
        temp=tiffobj.getTag('ImageDescription');
        a=regexp(temp,'(images|channels|slices|frames|spacing|min|max)=([\d|.]*)','tokens');
        for p=1:numel(a)
            val=str2double(a{p}{2});
            switch a{p}{1}
                case 'images'
                    handles.VAL_IMGNUM.String=val;
                    handles.VAL_IMGNUM.Value=val;
                    nimg=val;
                case 'channels'
                    handles.VAL_IMGCHANNEL.String=val;
                    handles.VAL_IMGCHANNEL.Value=val;
                    nch=val;
                case 'slices'
                    handles.VAL_IMGSLICE.String=val;
                    handles.VAL_IMGSLICE.Value=val;
                    nslice=val;
                case 'frames'
                    handles.VAL_IMGFRAME.String=val;
                    handles.VAL_IMGFRAME.Value=val;
                    nframe=val;
                case 'spacing'
                    handles.VAL_IMGSPACING.String=val;
                    handles.VAL_IMGSPACING.Value=val;
                    handles.VAL_RESZ.String=val;
                    handles.VAL_RESZ.Value=val;
                    rawdata.dz=val;
                case 'min'
                    handles.VAL_IMGMIN.String=val;
                    handles.VAL_IMGMIN.Value=val;
                case 'max'
                    handles.VAL_IMGMAX.String=val;
                    handles.VAL_IMGMAX.Value=val;
            end
        end
        rawdata.z=linspace(0,(nslice-1)*rawdata.dz,nslice);
        rawdata.val=zeros(nch,xdim,ydim,nslice,nframe,cat(2,'uint',num2str(bitdim)));
        imgidx=1;chidx=1;zidx=1;fidx=1;
        waitbar_handle = waitbar(0,'Please wait...',...
            'Name','Loading TIFF image',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)',...
            'WindowStyle','normal',...
            'Color',[0.2,0.2,0.2]);
        setappdata(waitbar_handle,'canceling',0);
        while ~lastDirectory(tiffobj)
            if getappdata(waitbar_handle,'canceling')
                fprintf('image loading cancelled\n');
                delete(waitbar_handle);       % DELETE the waitbar; don't try to CLOSE it.
                return;
            end
            rawdata.val(chidx,:,:,zidx,fidx)=tiffobj.read();
            nextDirectory(tiffobj);
            if chidx<nch
                chidx=chidx+1;
            else
                chidx=1;
                if zidx<nslice
                    zidx=zidx+1;
                else
                    zidx=1;
                    if fidx<nframe
                        fidx=fidx+1;
                    else
                        disp('extra frames');
                    end
                end
            end
            imgidx=imgidx+1;
            % Report current estimate in the waitbar's message field
            done=imgidx/nimg;
            waitbar(done,waitbar_handle,sprintf('%3.1f%%',100*done));
        end
        delete(waitbar_handle);       % DELETE the waitbar; don't try to CLOSE it.
        handles.SLIDER_C.Visible='on';handles.VAL_C.Visible='on';
        handles.SLIDER_C.Max=nch;
        handles.SLIDER_C.Value=1;
        handles.VAL_C.Value=handles.SLIDER_C.Value;
        handles.VAL_C.String=handles.SLIDER_C.Value;
        handles.SLIDER_C.SliderStep=[1/(nch-1),1];
        handles.SLIDER_Z.Visible='on';handles.VAL_Z.Visible='on';
        handles.SLIDER_Z.Max=nslice;
        handles.SLIDER_Z.Value=1;
        handles.VAL_Z.Value=handles.SLIDER_Z.Value;
        handles.VAL_Z.String=handles.SLIDER_Z.Value;
        handles.SLIDER_Z.SliderStep=[1/(nslice-1),1];
        handles.SLIDER_F.Visible='on';handles.VAL_F.Visible='on';
        handles.SLIDER_F.Max=nframe;
        handles.SLIDER_F.Value=1;
        handles.VAL_F.Value=handles.SLIDER_F.Value;
        handles.VAL_F.String=handles.SLIDER_F.Value;
        handles.SLIDER_F.SliderStep=[1/(nframe-1),1];
        
        imagesc(handles.PANEL_IMAGE,rawdata.x,rawdata.y,squeeze(rawdata.val(1,:,:,1,1)));
        colormap(handles.PANEL_IMAGE,'gray');
        localdata.val=cell(nch,nframe);
        localdata.colname={'image-ID','time-point','cycle','z-step','frame','accum','probe',...
            'photon-count','photon-count11','photon-count11','photon-count12','photon-count21','photon-count22',...
            'psfx','psfy','psfz','psf-photon-count','x','y','z','stdev','amp',...
            'background11','background12','background21','background22','maxResidualSlope',...
            'chisq','log-likelihood','llr','accuracy','fiducial','valid'};
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
        handles.SLIDER_IMGMIN.Max=2^bitdim-1;
        handles.SLIDER_IMGMIN.SliderStep=[1/handles.SLIDER_IMGMIN.Max,10/handles.SLIDER_IMGMIN.Max];
        SLIDER_IMGMIN_Callback(handles.SLIDER_IMGMIN,[],handles);
        handles.SLIDER_IMGMAX.Max=2^bitdim-1;
        handles.SLIDER_IMGMAX.SliderStep=[1/handles.SLIDER_IMGMAX.Max,10/handles.SLIDER_IMGMAX.Max];
        SLIDER_IMGMAX_Callback(handles.SLIDER_IMGMAX,[],handles);
        handles.FIGURE_IMG2LOC.Name=filename;
        msgbox('File Loaded');
    else
        msgbox('File Loading Cancelled');
    end
catch exception
    if exist('waitbar_handle','var')&&ishandle(waitbar_handle)
        delete(waitbar_handle);
    end
    errordlg(exception.message,'Calculation Error','modal');
end

% --- Executes on button press in BUTTON_OPEN.
function BUTTON_OPEN_Callback(~, ~, handles)
global rawdata localdata;
[pathname,~,~]=fileparts(handles.FIGURE_IMG2LOC.Name);
if isempty(pathname)
    pathname='./';
end
% ask for one file to open
[filename,pathname,~] = uigetfile({'*.i2l','image2localisation file  (*.i2l)';...
    '*.*','All Files (*.*)'},...
    'Select Saved Localisation Analysis File',...
    'MultiSelect','off',pathname);
% if files selected
if pathname~=0
    temp = load(cat(2,pathname,filename),'-mat'); % load file
    rawdata=temp.rawdata;
    localdata=temp.localdata;
    info=temp.info;
    handles.VAL_IMGLENGTH.String=info.imglen;
    handles.VAL_IMGWIDTH.String=info.imgwidth;
    handles.VAL_IMGBITDEPTH.String=info.bitdim;
    handles.VAL_IMGNUM.String=info.nimg;
    handles.VAL_IMGCHANNEL.String=info.nch;
    handles.VAL_IMGSLICE.String=info.nslice;
    handles.VAL_IMGFRAME.String=info.nframe;
    handles.VAL_IMGMIN.String=info.minval;
    handles.VAL_IMGMAX.String=info.maxval;
    handles.VAL_IMGXRES.String=info.xres;
    handles.VAL_IMGYRES.String=info.yres;
    handles.VAL_IMGSPACING.String=info.spacing;
    handles.VAL_NNPIXEL.String=info.nn;handles.VAL_NNPIXEL.Value=str2double(handles.VAL_NNPIXEL.String);
    handles.VAL_ZNNPIXEL.String=info.znn;handles.VAL_ZNNPIXEL.Value=str2double(handles.VAL_ZNNPIXEL.String);
    handles.VAL_THRESHOLD.String=info.threshold;handles.VAL_THRESHOLD.Value=str2double(handles.VAL_THRESHOLD.String);
    handles.VAL_RESX.String=info.resolution{1};handles.VAL_RESX.Value=str2double(handles.VAL_RESX.String);
    handles.VAL_RESY.String=info.resolution{2};handles.VAL_RESY.Value=str2double(handles.VAL_RESY.String);
    handles.VAL_RESZ.String=info.resolution{3};handles.VAL_RESZ.Value=str2double(handles.VAL_RESZ.String);
    handles.SLIDER_C.Visible='on';handles.VAL_C.Visible='on';
    handles.SLIDER_C.Max=str2double(info.nch);
    handles.SLIDER_C.Value=1;
    handles.VAL_C.Value=handles.SLIDER_C.Value;
    handles.VAL_C.String=handles.SLIDER_C.Value;
    handles.SLIDER_C.SliderStep=[1/(str2double(info.nch)-1),1];
    handles.SLIDER_Z.Visible='on';handles.VAL_Z.Visible='on';
    handles.SLIDER_Z.Max=str2double(info.nslice);
    handles.SLIDER_Z.Value=1;
    handles.VAL_Z.Value=handles.SLIDER_Z.Value;
    handles.VAL_Z.String=handles.SLIDER_Z.Value;
    handles.SLIDER_Z.SliderStep=[1/(str2double(info.nslice)-1),1];
    handles.SLIDER_F.Visible='on';handles.VAL_F.Visible='on';
    handles.SLIDER_F.Max=str2double(info.nframe);
    handles.SLIDER_F.Value=1;
    handles.VAL_F.Value=handles.SLIDER_F.Value;
    handles.VAL_F.String=handles.SLIDER_F.Value;
    handles.SLIDER_F.SliderStep=[1/(str2double(info.nframe)-1),1];
    imagesc(handles.PANEL_IMAGE,rawdata.x,rawdata.y,squeeze(rawdata.val(1,:,:,1,1)));
    colormap(handles.PANEL_IMAGE,'gray');
    handles.SLIDER_IMGMIN.Max=2^str2double(info.bitdim)-1;
    handles.SLIDER_IMGMIN.SliderStep=[1/handles.SLIDER_IMGMIN.Max,10/handles.SLIDER_IMGMIN.Max];
    SLIDER_IMGMIN_Callback(handles.SLIDER_IMGMIN,[],handles);
    handles.SLIDER_IMGMAX.Max=2^str2double(info.bitdim)-1;
    handles.SLIDER_IMGMAX.SliderStep=[1/handles.SLIDER_IMGMAX.Max,10/handles.SLIDER_IMGMAX.Max];
    SLIDER_IMGMAX_Callback(handles.SLIDER_IMGMAX,[],handles);
    handles.FIGURE_IMG2LOC.Name=cat(2,pathname,filename);
    msgbox(sprintf('%s successfully loaded\n',filename),'Open Image 2 Localisation Analysis File','modal');
end

% --- Executes on button press in BUTTON_SAVE.
function BUTTON_SAVE_Callback(~, ~, handles)
global rawdata localdata;%#ok<NUSED>
[pathname,~,~]=fileparts(handles.FIGURE_IMG2LOC.Name);
if isempty(pathname)
    pathname='./';
end
[filename,pathname,~]=uiputfile({'*.i2l','image2localisation file (*.i2l)';...
    '*.*','All Files (*.*)'},...
    'Select Saved Image 2 Localisation Analysis File',pathname);
if pathname~=0
    filename=cat(2,pathname,filename);
    version='7.3';
    info.imglen=handles.VAL_IMGLENGTH.String;
    info.imgwidth=handles.VAL_IMGWIDTH.String;
    info.bitdim=handles.VAL_IMGBITDEPTH.String;
    info.nimg=handles.VAL_IMGNUM.String;
    info.nch=handles.VAL_IMGCHANNEL.String;
    info.nslice=handles.VAL_IMGSLICE.String;
    info.nframe=handles.VAL_IMGFRAME.String;
    info.minval=handles.VAL_IMGMIN.String;
    info.maxval=handles.VAL_IMGMAX.String;
    info.xres=handles.VAL_IMGXRES.String;
    info.yres=handles.VAL_IMGYRES.String;
    info.spacing=handles.VAL_IMGSPACING.String;
    info.nn=handles.VAL_NNPIXEL.String;
    info.znn=handles.VAL_ZNNPIXEL.String;
    info.threshold=handles.VAL_THRESHOLD.String;
    info.resolution={handles.VAL_RESX.String,handles.VAL_RESY.String,handles.VAL_RESZ.String}; %#ok<STRNU>
    save(filename,'rawdata','localdata','info','-mat',cat(2,'-v',version));
    handles.LOCSUPRES.Name=filename;
    msgbox(sprintf('%s saved in ver %s\n',filename,version),'Save File','modal');
else
    % user interuption
    msgbox(sprintf('File save cancelled\n'),'Save File','modal');
end

%-------------------------------------------------------------------------
function VAL_NNPIXEL_Callback(hObject, ~, ~)
val=str2double(hObject.String);
val=min(max(val,1),5);
hObject.String=val;
hObject.Value=val;

function VAL_THRESHOLD_Callback(hObject, ~, handles)
val=str2double(hObject.String);
val=min(max(val,0),handles.SLIDER_IMGMAX.Max);
hObject.String=val;
hObject.Value=val;
updateimage(handles);

function VAL_RESX_Callback(hObject, eventdata, handles)

function VAL_RESY_Callback(hObject, eventdata, handles)

function VAL_RESZ_Callback(hObject, eventdata, handles)

%-------------------------------------------------------------------------
function SLIDER_C_Callback(~, ~, handles)
handles.SLIDER_C.Value=round(handles.SLIDER_C.Value);
handles.VAL_C.String=handles.SLIDER_C.Value;
updateimage(handles);

function SLIDER_C_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
hObject.Visible='off';
hObject.Min=1;
hObject.Value=1;

function SLIDER_Z_Callback(~, ~, handles)
handles.SLIDER_Z.Value=round(handles.SLIDER_Z.Value);
handles.VAL_Z.String=handles.SLIDER_Z.Value;
updateimage(handles);

function SLIDER_Z_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
hObject.Visible='off';
hObject.Min=1;
hObject.Value=1;

function SLIDER_F_Callback(~, ~, handles)
handles.SLIDER_F.Value=round(handles.SLIDER_F.Value);
handles.VAL_F.String=handles.SLIDER_F.Value;
updateimage(handles);

function SLIDER_F_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
hObject.Visible='off';
hObject.Min=1;
hObject.Value=1;

function VAL_C_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
hObject.Visible='off';

function VAL_Z_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
hObject.Visible='off';

function VAL_F_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
hObject.Visible='off';

function SLIDER_IMGMAX_Callback(hObject, ~, handles)
hObject.Value=max(hObject.Value,handles.SLIDER_IMGMIN.Value);
handles.PANEL_IMAGE.CLim(2)=hObject.Value;
handles.VAL_MAX.String=hObject.Value;

function SLIDER_IMGMAX_CreateFcn(hObject, ~, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
hObject.Min=0;hObject.Max=255;hObject.Value=255;
handles.VAL_MAX.String=hObject.Value;
hObject.SliderStep=[1/hObject.Max,10/hObject.Max];

function SLIDER_IMGMIN_Callback(hObject, ~, handles)
hObject.Value=min(hObject.Value,handles.SLIDER_IMGMAX.Value);
handles.PANEL_IMAGE.CLim(1)=hObject.Value;
handles.VAL_MIN.String=hObject.Value;

function SLIDER_IMGMIN_CreateFcn(hObject, ~, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
hObject.Min=0;hObject.Max=255;hObject.Value=0;
hObject.SliderStep=[1/hObject.Max,10/hObject.Max];
handles.VAL_MIN.String=hObject.Value;

%-------------------------------------------------------------------------
% --- Executes on button press in BUTTON_LOCALISEFRAME.
function BUTTON_LOCALISEFRAME_Callback(hObject, ~, handles)
global rawdata localdata;
try
    nn=handles.VAL_NNPIXEL.Value;
    znn=handles.VAL_ZNNPIXEL.Value;
    cidx=handles.SLIDER_C.Value;
    fidx=handles.SLIDER_F.Value;
    threshold=handles.VAL_THRESHOLD.Value;
    subimgsize=[2*nn+1,2*nn+1,2*znn+1];
    nx=numel(rawdata.x)-nn*2;
    ny=numel(rawdata.y)-nn*2;
    nz=numel(rawdata.z)-znn*2;
    loopsize=nx*ny*nz;
    xstart=1+nn;xstop=numel(rawdata.x)-nn;
    ystart=1+nn;ystop=numel(rawdata.y)-nn;
    zstart=1+znn;zstop=numel(rawdata.z)-znn;
    [i,j,k]=ind2sub(subimgsize,1:prod(subimgsize));
    i=i(:);j=j(:);k=k(:);
    ds=[rawdata.dx,rawdata.dy,rawdata.dz];
    hObject.String='Calculating';pause(0.1);
    fprintf(1,'start...\n');
    local=zeros(loopsize,3);
    intensity=zeros(loopsize,1);
    chisq=zeros(loopsize,1);
    localidx=1;
    tempval=squeeze(rawdata.val(cidx,:,:,:,fidx));
    tempval(tempval<threshold)=0;
    for zpix=zstart:1:zstop
        fprintf(1,'frame %g, channel %g, step %g/%g...\n',fidx,cidx,zpix-nn,nz);
        zorg=rawdata.z(zpix-znn);
        for xpix=xstart:1:xstop
            xorg=rawdata.x(xpix-nn);
            for ypix=ystart:1:ystop
                yorg=rawdata.y(ypix-nn);
                val=tempval(xpix-nn:xpix+nn,ypix-nn:ypix+nn,zpix-znn:zpix+znn);
                Imedian=median(val(:));%need to decide on this criteria
                if Imedian>0
                    [intensity(localidx),local(localidx,1:3)]=roicentroid([i,j,k,val(:)],ds,[xorg,yorg,zorg]);
                    chisq(localidx)=Imedian/mean(val(:));
                else
                    local(localidx,1:3)=[nan,nan,nan];
                    intensity(localidx)=0;
                    chisq(localidx)=1;
                end
                localidx=localidx+1;
            end
        end
    end
    % remove invalid points
    invalid=(intensity==0);
    local(invalid,:)=[];
    intensity(invalid)=[];
    chisq(invalid)=[];
    %x,y,z,photon-count
    localdata.val{cidx,fidx}=[local(:,[1,2,3]),intensity,chisq];
    fprintf(1,'finished\n');
    beep;beep;
    hObject.String='Localise Frame';
    updateimage(handles);
    msgbox(sprintf('Found %g localisation for Channel%g Frame%g.',numel(intensity),cidx,fidx),'Localisation Finished');
catch exception
    errordlg(exception.message,'Calculation Error','modal');
end

function BUTTON_CLEARLOCAL_Callback(~, ~, handles)
global localdata;
cidx=handles.SLIDER_C.Value;
fidx=handles.SLIDER_F.Value;
localdata.val{cidx,fidx}=[];
updateimage(handles);
msgbox(sprintf('Localisation for Channel%g Frame%g cleared.',cidx,fidx),'Localisation Cleard');

function BUTTON_EXPORTFRAME_Callback(~, ~, handles)
global localdata;
[pathname,filename,~]=fileparts(handles.FIGURE_IMG2LOC.Name);
if isempty(pathname)
    pathname=pwd;
end
fidx=handles.SLIDER_F.Value;
filename=cat(2,filename,'_frame',num2str(fidx));
% ask for one file to open
[filename,pathname,~] = uiputfile({'*.csv','Bruker export ascii file (*.csv)';...
    '*.*','All Files (*.*)'},...
    'Select Localisation ASCII File',cat(2,pathname,filesep,filename));
if ischar(filename)
    filename=cat(2,pathname,filename);
    threshold=handles.VAL_THRESHOLD.Value;
    locsize=cellfun(@(x)size(x,1),localdata.val);
    tabval=zeros(sum(locsize(:)),numel(localdata.colname));
    nch=size(localdata.val,1);
    startidx=1;
    for cidx=1:nch
        if ~isempty(localdata.val{cidx,fidx})
            endidx=size(localdata.val{cidx,fidx},1)+startidx-1;
            tabval(startidx:endidx,strcmp(localdata.colname,'probe'))=cidx-1;
            tabval(startidx:endidx,strcmp(localdata.colname,'x'))=localdata.val{cidx,fidx}(:,1)*1000;
            tabval(startidx:endidx,strcmp(localdata.colname,'y'))=localdata.val{cidx,fidx}(:,2)*1000;
            tabval(startidx:endidx,strcmp(localdata.colname,'z'))=localdata.val{cidx,fidx}(:,3)*1000;
            tabval(startidx:endidx,strcmp(localdata.colname,'psf-photon-count'))=localdata.val{cidx,fidx}(:,4);
            tabval(startidx:endidx,strcmp(localdata.colname,'psfx'))=ones(size(localdata.val{cidx,fidx},1),1)*0.025;
            tabval(startidx:endidx,strcmp(localdata.colname,'psfy'))=ones(size(localdata.val{cidx,fidx},1),1)*0.025;
            tabval(startidx:endidx,strcmp(localdata.colname,'psfz'))=ones(size(localdata.val{cidx,fidx},1),1)*0.05;
            tabval(startidx:endidx,strcmp(localdata.colname,'background11'))=ones(size(localdata.val{cidx,fidx},1),1)*threshold;
            tabval(startidx:endidx,strcmp(localdata.colname,'background12'))=ones(size(localdata.val{cidx,fidx},1),1)*threshold;
            tabval(startidx:endidx,strcmp(localdata.colname,'background21'))=ones(size(localdata.val{cidx,fidx},1),1)*threshold;
            tabval(startidx:endidx,strcmp(localdata.colname,'background22'))=ones(size(localdata.val{cidx,fidx},1),1)*threshold;
            tabval(startidx:endidx,strcmp(localdata.colname,'chisq'))=localdata.val{cidx,fidx}(:,5);
            tabval(startidx:endidx,strcmp(localdata.colname,'valid'))=ones(size(localdata.val{cidx,fidx},1),1);
            tabval(startidx:endidx,strcmp(localdata.colname,'accuracy'))=ones(size(localdata.val{cidx,fidx},1),1)*15;
            startidx=endidx+1;
        end
    end
    %csvwrite(filename,tabval);
    fileID = fopen(filename,'w');
    headerfmt=[repmat('%s,',1,numel(localdata.colname)-1),'%s\n'];
    filefmt=[repmat('%0.6g,',1,numel(localdata.colname)-1),'%0.6g\n'];
    fprintf(fileID,headerfmt,localdata.colname{:});
    fprintf(fileID,filefmt,tabval');
    fclose(fileID);
    msgbox('CSV file exported');
else
    msgbox('File Exporting Cancelled');
end

function BUTTON_LOCALISEALL_Callback(hObject, ~, handles)
global rawdata localdata;
try
    nn=handles.VAL_NNPIXEL.Value;
    znn=handles.VAL_ZNNPIXEL.Value;
    %cidx=handles.SLIDER_C.Value;
    nch=handles.SLIDER_C.Max;
    %fidx=handles.SLIDER_F.Value;
    nframe=handles.SLIDER_F.Max;
    threshold=handles.VAL_THRESHOLD.Value;
    subimgsize=[2*nn+1,2*nn+1,2*znn+1];
    nvoxel=prod(subimgsize);
    nx=numel(rawdata.x)-nn*2;
    ny=numel(rawdata.y)-nn*2;
    nz=numel(rawdata.z)-znn*2;
    loopsize=nx*ny*nz;
    xstart=1+nn;xstop=numel(rawdata.x)-nn;
    ystart=1+nn;ystop=numel(rawdata.y)-nn;
    zstart=1+znn;zstop=numel(rawdata.z)-znn;
    [i,j,k]=ind2sub(subimgsize,1:nvoxel);
    i=i(:);j=j(:);k=k(:);
    ds=[rawdata.dx,rawdata.dy,rawdata.dz];
    hObject.String='Calculating';pause(0.1);
    fprintf(1,'start...\n');
    for fidx=1:1:nframe
        for cidx=1:1:nch
            local=zeros(loopsize,3);
            intensity=zeros(loopsize,1);
            chisq=zeros(loopsize,1);
            localidx=1;
            tempval=squeeze(rawdata.val(cidx,:,:,:,fidx));
            tempval(tempval<threshold)=0;
            for zpix=zstart:1:zstop
                fprintf(1,'frame %g, channel %g, step %g/%g...\n',fidx,cidx,zpix-nn,nz);
                zorg=rawdata.z(zpix-znn);
                for xpix=xstart:1:xstop
                    xorg=rawdata.x(xpix-nn);
                    for ypix=ystart:1:ystop
                        yorg=rawdata.y(ypix-nn);
                        val=tempval(xpix-nn:xpix+nn,ypix-nn:ypix+nn,zpix-znn:zpix+znn);
                        Imedian=median(val(:));%need to decide on this criteria
                        if Imedian>0
                            [intensity(localidx),local(localidx,1:3)]=roicentroid([i,j,k,val(:)],ds,[xorg,yorg,zorg]);
                            chisq(localidx)=Imedian/(intensity(localidx)/nvoxel);
                        else
                            local(localidx,1:3)=[nan,nan,nan];
                            intensity(localidx)=0;
                            chisq(localidx)=1;
                        end
                        localidx=localidx+1;
                    end
                end
            end
            % remove invalid points
            invalid=(intensity==0);
            local(invalid,:)=[];
            intensity(invalid)=[];
            chisq(invalid)=[];
            %x,y,z,photon-count
            localdata.val{cidx,fidx}=[local(:,[1,2,3]),intensity,chisq];
        end
    end
    fprintf(1,'finished\n');
    beep;beep;
    hObject.String='Localise ALL';
    updateimage(handles);
    msgbox(sprintf('Found %g localisation for Channel%g Frame%g.',numel(intensity),cidx,fidx),'Localisation Finished');
catch exception
    errordlg(exception.message,'Calculation Error','modal');
end

% --- Executes on button press in BUTTON_EXPORTALL.
function BUTTON_EXPORTALL_Callback(~, ~, handles)
global localdata;
[pathname,filename,~]=fileparts(handles.FIGURE_IMG2LOC.Name);
if isempty(pathname)
    pathname=pwd;
end
filename=cat(2,filename,'_local');
% ask for one file to open
[filename,pathname,~] = uiputfile({'*.csv','Bruker export ascii file (*.csv)';...
    '*.*','All Files (*.*)'},...
    'Select Localisation ASCII File',cat(2,pathname,filesep,filename));
if ischar(filename)
    for fidx=1:size(localdata.val,2)
        [~,fname,ext]=fileparts(filename);
        fname=cat(2,pathname,fname,'_frame',num2str(fidx),ext);
        threshold=handles.VAL_THRESHOLD.Value;
        locsize=cellfun(@(x)size(x,1),localdata.val);
        tabval=zeros(sum(locsize(:,fidx)),numel(localdata.colname));
        nch=size(localdata.val,1);
        startidx=1;
        for cidx=1:nch
            if ~isempty(localdata.val{cidx,fidx})
                endidx=size(localdata.val{cidx,fidx},1)+startidx-1;
                tabval(startidx:endidx,strcmp(localdata.colname,'probe'))=cidx-1;
                tabval(startidx:endidx,strcmp(localdata.colname,'x'))=localdata.val{cidx,fidx}(:,1)*1000;
                tabval(startidx:endidx,strcmp(localdata.colname,'y'))=localdata.val{cidx,fidx}(:,2)*1000;
                tabval(startidx:endidx,strcmp(localdata.colname,'z'))=localdata.val{cidx,fidx}(:,3)*1000;
                tabval(startidx:endidx,strcmp(localdata.colname,'psf-photon-count'))=localdata.val{cidx,fidx}(:,4);
                tabval(startidx:endidx,strcmp(localdata.colname,'psfx'))=ones(size(localdata.val{cidx,fidx},1),1)*0.025;
                tabval(startidx:endidx,strcmp(localdata.colname,'psfy'))=ones(size(localdata.val{cidx,fidx},1),1)*0.025;
                tabval(startidx:endidx,strcmp(localdata.colname,'psfz'))=ones(size(localdata.val{cidx,fidx},1),1)*0.05;
                tabval(startidx:endidx,strcmp(localdata.colname,'background11'))=ones(size(localdata.val{cidx,fidx},1),1)*threshold;
                tabval(startidx:endidx,strcmp(localdata.colname,'background12'))=ones(size(localdata.val{cidx,fidx},1),1)*threshold;
                tabval(startidx:endidx,strcmp(localdata.colname,'background21'))=ones(size(localdata.val{cidx,fidx},1),1)*threshold;
                tabval(startidx:endidx,strcmp(localdata.colname,'background22'))=ones(size(localdata.val{cidx,fidx},1),1)*threshold;
                tabval(startidx:endidx,strcmp(localdata.colname,'chisq'))=localdata.val{cidx,fidx}(:,5);
                tabval(startidx:endidx,strcmp(localdata.colname,'valid'))=ones(size(localdata.val{cidx,fidx},1),1);
                tabval(startidx:endidx,strcmp(localdata.colname,'accuracy'))=ones(size(localdata.val{cidx,fidx},1),1)*15;
                startidx=endidx+1;
            end
        end
        %csvwrite(filename,tabval);
        headerfmt=[repmat('%s,',1,numel(localdata.colname)-1),'%s\n'];
        filefmt=[repmat('%0.6g,',1,numel(localdata.colname)-1),'%0.6g\n'];
        fileID = fopen(fname,'w');
        fprintf(fileID,headerfmt,localdata.colname{:});
        fprintf(fileID,filefmt,tabval');
        fclose(fileID);
    end
    msgbox('CSV files exported');
else
    msgbox('File Exporting Cancelled');
end
%-------------------------------------------------------------------------
% --- Executes when user attempts to close FIGURE_IMG2LOC.
function FIGURE_IMG2LOC_CloseRequestFcn(hObject, ~, ~)
global tiffobj rawdata; %#ok<NUSED>
clear global rawdata;
if ~isempty(tiffobj)
    tiffobj.close;
end
clear global tiffobj;
delete(hObject);

%=========================================================================
function updateimage(handles)
global rawdata localdata;
cidx=handles.SLIDER_C.Value;
zidx=handles.SLIDER_Z.Value;
fidx=handles.SLIDER_F.Value;
threshold=handles.VAL_THRESHOLD.Value;
val=squeeze(rawdata.val(cidx,:,:,zidx,fidx));
val(val<threshold)=nan;
imshow(val',[handles.SLIDER_IMGMIN.Value,handles.SLIDER_IMGMAX.Value],'Parent',handles.PANEL_IMAGE,'XData',rawdata.x,'YData',rawdata.y);
axis(handles.PANEL_IMAGE,'on');
if ~isempty(localdata.val{cidx,fidx})
    coloridx={'y','r','b','g'};
    if isempty(handles.PANEL_LOCALISATION.Children)
        plot3(handles.PANEL_LOCALISATION,...
            localdata.val{cidx,fidx}(:,1),localdata.val{cidx,fidx}(:,2),localdata.val{cidx,fidx}(:,3),...
            'LineStyle','none','Marker','.','Color',coloridx{cidx},...
            'MarkerSize',1);
        set(handles.PANEL_LOCALISATION,'Color',[0,0,0]);
        xlim(handles.PANEL_LOCALISATION,[rawdata.x(1),rawdata.x(end)]);
        xlabel(handles.PANEL_LOCALISATION,'X');
        ylim(handles.PANEL_LOCALISATION,[rawdata.y(1),rawdata.y(end)]);
        ylabel(handles.PANEL_LOCALISATION,'Y');
        zlim(handles.PANEL_LOCALISATION,[rawdata.z(1),rawdata.z(end)]);
        zlabel(handles.PANEL_LOCALISATION,'Z');
        view(handles.PANEL_LOCALISATION,[0 0 -90]);
        daspect(handles.PANEL_LOCALISATION,[1 1 1]);
    else
        handles.PANEL_LOCALISATION.Children.XData=localdata.val{cidx,fidx}(:,1);
        handles.PANEL_LOCALISATION.Children.YData=localdata.val{cidx,fidx}(:,2);
        handles.PANEL_LOCALISATION.Children.ZData=localdata.val{cidx,fidx}(:,3);
        view(handles.PANEL_LOCALISATION,[0 0 -90]);
    end
else
    if ~isempty(handles.PANEL_LOCALISATION.Children)
        handles.PANEL_LOCALISATION.Children.XData=[];
        handles.PANEL_LOCALISATION.Children.YData=[];
        handles.PANEL_LOCALISATION.Children.ZData=[];
    end
end

function [Itotal,coord]=roicentroid(val,dr,orig)
Itotal=sum(val(:,4));
coord=sum(bsxfun(@times,val(:,1:3),val(:,4)),1)/Itotal.*dr+orig;

%{
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
%}
