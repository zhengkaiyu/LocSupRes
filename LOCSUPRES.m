function varargout = LOCSUPRES(varargin)
% LOCSUPRES MATLAB code for LOCSUPRES.fig
%      LOCSUPRES, by itself, creates a new LOCSUPRES or raises the existing
%      singleton*.
%
%      H = LOCSUPRES returns the handle to a new LOCSUPRES or the handle to
%      the existing singleton*.
%
%      LOCSUPRES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOCSUPRES.M with the given input arguments.
%
%      LOCSUPRES('Property','Value',...) creates a new LOCSUPRES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LOCSUPRES_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LOCSUPRES_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LOCSUPRES

% Last Modified by GUIDE v2.5 02-Nov-2017 18:18:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LOCSUPRES_OpeningFcn, ...
    'gui_OutputFcn',  @LOCSUPRES_OutputFcn, ...
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

% --- Executes just before LOCSUPRES is made visible.
function LOCSUPRES_OpeningFcn(hObject, ~, handles, varargin)
% Choose default command line output for LOCSUPRES
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialise(handles);

% --- Outputs from this function are returned to the command line.
function varargout = LOCSUPRES_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%==GUI functions==========================================================
% --- Executes on button press in BUTTON_OPEN.
function BUTTON_OPEN_Callback(~, ~, handles) %#ok<*DEFNU>
global DATA;
[pathname,~,~]=fileparts(handles.LOCSUPRES.Name);
if isempty(pathname)
    pathname='./';
end
% ask for one file to open
[filename,pathname,~] = uigetfile({'*.laf','Localisation analysis file (*.laf)';...
    '*.*','All Files (*.*)'},...
    'Select Saved Localisation Analysis File',...
    'MultiSelect','off',pathname);
% if files selected
if pathname~=0
    temp = load(cat(2,pathname,filename),'-mat'); % load file
    DATA=temp.DATA;
    % update datainfo table
    display_datainfo( handles.TABLE_DATAINFO );
    % update probe menu
    fname=fieldnames(DATA.datainfo);
    fidx=find(cellfun(@(x)~isempty(x),regexp(fname,'probe\w_name')));
    probe_list=cell(DATA.datainfo.n_probe,1);
    for idx=1:DATA.datainfo.n_probe
        probe_list{idx}=DATA.datainfo.(fname{fidx(idx)});
    end
    handles.MENU_PROBE.String=probe_list;
    handles.MENU_PROBE.Value=1;
    handles.MENU_SORTCLUSTERTABLE.Value=1;
    handles.LOCSUPRES.Name=cat(2,pathname,filename);
    msgbox(sprintf('%s successfully loaded\n',filename),'Open Localisation Analysis File','modal');
end

% --- Executes on button press in BUTTON_SAVE.
function BUTTON_SAVE_Callback(~, ~, handles)
global DATA; %#ok<NUSED>
[pathname,~,~]=fileparts(handles.LOCSUPRES.Name);
if isempty(pathname)
    pathname='./';
end
[filename,pathname,~]=uiputfile({'*.laf','Localisation analysis file (*.laf)';...
    '*.*','All Files (*.*)'},...
    'Select Saved Localisation Analysis File',pathname);
if pathname~=0
    filename=cat(2,pathname,filename);
    version='7.3';
    save(filename,'DATA','-mat',cat(2,'-v',version));
    handles.LOCSUPRES.Name=filename;
    msgbox(sprintf('%s saved in ver %s\n',filename,version),'Save File','modal');
else
    % user interuption
    msgbox(sprintf('File save cancelled\n'),'Save File','modal');
end

% --- Executes on button press in BUTTON_LOADRAW.
function BUTTON_LOADRAW_Callback(~, ~, handles)
global DATA DEFAULT_COLOUR;
[pathname,~,~]=fileparts(handles.LOCSUPRES.Name);
if isempty(pathname)
    pathname='./';
end
% select file from storage
[filename,data_pathname,~]=uigetfile({'*.*','All Files (*.*)';...
    '*.csv','Bruker exported localisation data (*.csv)'},...
    'Select Raw Data File','MultiSelect','off',pathname);
if data_pathname~=0     %if files selected
    filename=cat(2,data_pathname,filename);
    % attempt to open file
    [data_file,~]=fopen(filename,'r');
    if data_file>=3 % successfully opened
        % read header
        buffer=fgetl(data_file);
        header=regexp(buffer,'[\s,\s]','split');
        for header_idx=1:numel(header)
            DATA.datainfo.(cat(2,'header',num2str(header_idx)))=header{header_idx};
        end
        fclose(data_file);% close file
        % read data segment
        val=csvread(filename,1,0);
        % image_id column
        DATA.datainfo.imageidcol=find(strcmp('image-ID',header));
        image_ids=unique(val(:,DATA.datainfo.imageidcol));
        DATA.datainfo.n_image=numel(image_ids);
        % cycle column
        DATA.datainfo.cyclecol=find(strcmp('cycle',header));
        cycles=unique(val(:,DATA.datainfo.cyclecol));
        DATA.datainfo.n_cycle=numel(cycles);
        % zstep column
        DATA.datainfo.zstepcol=find(strcmp('z-step',header));
        zsteps=unique(val(:,DATA.datainfo.zstepcol));
        DATA.datainfo.n_zstep=numel(zsteps);
        % frame column
        DATA.datainfo.framecol=find(strcmp('frame',header));
        frames=unique(val(:,DATA.datainfo.framecol));
        DATA.datainfo.n_frame=numel(frames);
        % accum column
        DATA.datainfo.accumcol=find(strcmp('accum',header));
        accums=unique(val(:,DATA.datainfo.accumcol));
        DATA.datainfo.n_accum=numel(accums);
        % probe column
        DATA.datainfo.probecol=find(strcmp('probe',header));
        probes=unique(val(:,DATA.datainfo.probecol));
        DATA.datainfo.n_probe=numel(probes);
        
        % valid column (column 28, last column)
        DATA.datainfo.validcol=find(strcmp('valid',header));
        valid=find(val(:,DATA.datainfo.validcol)==1);
        DATA.datainfo.n_valid=numel(valid);
        
        % other useful column designation
        DATA.datainfo.psfpccol=find(strcmp('psf-photon-count',header));
        DATA.datainfo.bg11col=find(strcmp('background11',header));
        DATA.datainfo.bg12col=find(strcmp('background12',header));
        DATA.datainfo.bg21col=find(strcmp('background21',header));
        DATA.datainfo.bg22col=find(strcmp('background22',header));
        DATA.datainfo.chisqcol=find(strcmp('chisq',header));
        DATA.datainfo.loglhcol=find(strcmp('log-likelihood',header));
        DATA.datainfo.accuracycol=find(strcmp('accuracy',header));
        
        % convert position from nm to um
        DATA.datainfo.psfxcol=find(strcmp('psfx',header));
        DATA.datainfo.psfycol=find(strcmp('psfy',header));
        DATA.datainfo.psfzcol=find(strcmp('psfz',header));
        DATA.datainfo.xcol=find(strcmp('x',header));
        DATA.datainfo.ycol=find(strcmp('y',header));
        DATA.datainfo.zcol=find(strcmp('z',header));
        X_pos=val(valid,DATA.datainfo.xcol)/1000;
        Y_pos=val(valid,DATA.datainfo.ycol)/1000;
        Z_pos=val(valid,DATA.datainfo.zcol)/1000;
        
        % set dimension scales
        p_scale=probes;
        X_scale=min(X_pos):DATA.datainfo.dX:max(X_pos);
        Y_scale=min(Y_pos):DATA.datainfo.dY:max(Y_pos);
        Z_scale=min(Z_pos):DATA.datainfo.dZ:max(Z_pos);
        T_scale=1;
        % get dimension sizes
        p_size=numel(p_scale);
        X_size=numel(X_scale);
        Y_size=numel(Y_scale);
        Z_size=numel(Z_scale);
        T_size=numel(T_scale);
        DATA.datainfo.data_dim=[p_size,X_size,Y_size,Z_size,T_size];
        [~,DATA.datainfo.dataname,~]=fileparts(filename);
        DATA.datainfo.p=p_scale';
        DATA.datainfo.X=X_scale';
        DATA.datainfo.Y=Y_scale';
        DATA.datainfo.Z=Z_scale';
        DATA.datainfo.T=T_scale';
        DATA.datainfo.probe_colours=DEFAULT_COLOUR(1:p_size);
        probe_list=cell(p_size,1);
        DATA.datainfo.min_cluster_size=repmat(DATA.datainfo.min_cluster_size(1),1,DATA.datainfo.n_probe);
        DATA.datainfo.localdist=repmat(DATA.datainfo.localdist(1),1,DATA.datainfo.n_probe);
        DATA.dataval=val(valid,:);
        % convert position from nm to um
        DATA.dataval(:,[DATA.datainfo.psfxcol,DATA.datainfo.psfycol,DATA.datainfo.psfzcol,DATA.datainfo.xcol,DATA.datainfo.ycol,DATA.datainfo.zcol])=DATA.dataval(:,[DATA.datainfo.psfxcol,DATA.datainfo.psfycol,DATA.datainfo.psfzcol,DATA.datainfo.xcol,DATA.datainfo.ycol,DATA.datainfo.zcol])/1000;
        for probeidx=1:p_size
            probe_list{probeidx}=cat(2,'probe',num2str(p_scale(probeidx)));
            DATA.datainfo.(cat(2,'probe',num2str(p_scale(probeidx)),'_name'))=probe_list{probeidx};
            DATA.probe(probeidx).cluster=[];
            probe=DATA.dataval(:,DATA.datainfo.probecol)==DATA.datainfo.p(probeidx);
            DATA.probe(probeidx).location=DATA.dataval(probe,[DATA.datainfo.xcol,DATA.datainfo.ycol,DATA.datainfo.zcol]);
        end
        DATA.datainfo.last_change=datestr(now);
        
        % update datainfo table
        display_datainfo( handles.TABLE_DATAINFO );
        % update probe menu
        handles.MENU_PROBE.String=probe_list;
        handles.MENU_PROBE.Value=1;
        handles.MENU_SORTCLUSTERTABLE.Value=1;
        handles.LOCSUPRES.Name=filename;
        msgbox(sprintf('%s successfully loaded\n',filename),'Load Exported Localisation File','modal');
    else
        % failed to open file
        errordlg(sprintf('%s failed to open\n',filename));
    end
else
    % failed to open file
    errordlg(sprintf('fileopen cancelled\n'));
end

% --- Executes when entered data in editable cell(s) in TABLE_DATAINFO.
function TABLE_DATAINFO_CellEditCallback(hObject, eventdata, handles)
global DATA;
temp=get(hObject,'Data');
fname=temp(:,1);
switch fname{eventdata.Indices(1)}
    case 'probe_colours'
        DATA.datainfo.probe_colours=eventdata.NewData;
    case {'dX','dY','dZ'}
        DATA.datainfo.(fname{eventdata.Indices(1)})=str2double(eventdata.NewData);
    case 'scaled'
        DATA.datainfo.(fname{eventdata.Indices(1)})=(str2double(eventdata.NewData)==1);
    case 'localdist'
        val=str2num(eventdata.NewData); %#ok<*ST2NM>
        if numel(val)~=DATA.datainfo.n_probe
            val=repmat(val(1),1,DATA.datainfo.n_probe);
        end
        DATA.datainfo.(fname{eventdata.Indices(1)})=val;
    case 'min_cluster_size'
        val=str2num(eventdata.NewData);
        if numel(val)~=DATA.datainfo.n_probe
            val=repmat(val(1),1,DATA.datainfo.n_probe);
        end
        DATA.datainfo.(fname{eventdata.Indices(1)})=val;
    case {'synaptic_range','per_synapse'}
        DATA.datainfo.(fname{eventdata.Indices(1)})=str2double(eventdata.NewData);
    case {'slice_int','isoval'}
        val=str2num(eventdata.NewData);
        if numel(val)~=DATA.datainfo.n_probe
            val=repmat(val(1),1,DATA.datainfo.n_probe);
        end
        DATA.datainfo.(fname{eventdata.Indices(1)})=val;
    case {'probe0_name','probe1_name','probe2_name','probe3_name'}
        DATA.datainfo.(fname{eventdata.Indices(1)})=eventdata.NewData;
        probefname=fieldnames(DATA.datainfo);
        fidx=find(cellfun(@(x)~isempty(x),regexp(probefname,'probe\w_name')));
        probe_list=cell(DATA.datainfo.n_probe,1);
        for idx=1:DATA.datainfo.n_probe
            probe_list{idx}=DATA.datainfo.(probefname{fidx(idx)});
        end
        handles.MENU_PROBE.String=probe_list;
    case {'accum_threshold','psfx_threshold','psfy_threshold','psfz_threshold','snr_threshold','chisq_threshold','loglike_threshold','accuracy_threshold'}
        DATA.datainfo.(fname{eventdata.Indices(1)})=str2double(eventdata.NewData);
    otherwise
        errordlg(sprintf('%s cannot be changed',fname{eventdata.Indices(1)}),'Strict Variable','modal');
end
% update datainfo table
display_datainfo( hObject );

% --- Executes on selection change in MENU_PROBE.
function MENU_PROBE_Callback(hObject, ~, handles)
pidx=get(hObject,'Value');
% update cluster info table
display_clusterinfo(pidx,handles);

% --- Executes on button press in BUTTON_SCATTER.
function BUTTON_SCATTER_Callback(~, ~, handles)
showlocation(handles)

% --- Executes on button press in BUTTON_IMAGE.
function BUTTON_IMAGE_Callback(~, ~, handles)
showimage(handles)

% --- Executes on selection change in MENU_ANALYSOR.
function MENU_ANALYSOR_Callback(~, ~, ~)

% --- Executes on button press in BUTTON_CALCULATE.
function BUTTON_CALCULATE_Callback(~, ~, handles)
analysor=handles.MENU_ANALYSOR.String{handles.MENU_ANALYSOR.Value};
eval(cat(2,analysor,'(handles)'));

% --- Executes when selected cell(s) is changed in TABLE_CLUSTERINFO.
function TABLE_CLUSTERINFO_CellSelectionCallback(hObject, eventdata, handles)
global DATA;
probeidx=handles.MENU_PROBE.Value;
rowidx=unique(eventdata.Indices(:,1));
hObject.UserData=rowidx;
colidx=unique(eventdata.Indices(:,2));
currenthold=handles.PANEL_CLUSTER.NextPlot;
switch currenthold
    case 'replaceall'
        cla(handles.PANEL_CLUSTER);
        handles.PANEL_CLUSTER.NextPlot='add';
end
for idx=1:numel(rowidx)
    if colidx==1
        clusteridx=hObject.Data(rowidx(idx),1);
        centroid=DATA.probe(probeidx).cluster(clusteridx).centroid;
        plot3(handles.PANEL_CLUSTER,centroid(:,1),centroid(:,2),centroid(:,3),...
            'LineStyle','-','LineWidth',2,...
            'MarkerSize',5,'Marker','s',...
            'Color',DATA.datainfo.probe_colours(probeidx));
        showcluster(probeidx, clusteridx, handles);
    elseif colidx>10% selected on nn column
        temp=char(cellfun(@(x)char(x),(regexp(hObject.ColumnName(colidx(1)),'nnc\d_id','match')),'UniformOutput',false));
        if ~isempty(temp)
            nn_idx=find(DATA.datainfo.p==str2double(temp(4)));
            nnclusteridx=hObject.Data(rowidx(idx),colidx(1));
            centroid=DATA.probe(nn_idx).cluster(nnclusteridx).centroid;
            plot3(handles.PANEL_CLUSTER,centroid(:,1),centroid(:,2),centroid(:,3),...
                'LineStyle','-','LineWidth',2,...
                'MarkerSize',8,'Marker','s',...
                'Color',DATA.datainfo.probe_colours(nn_idx));
            showcluster(nn_idx, nnclusteridx, handles);
        else
            % show synapse
            temp=char(cellfun(@(x)char(x),(regexp(hObject.ColumnName(colidx(1)),'_synapse(','match')),'UniformOutput',false));
            if ~isempty(temp)
                clusteridx=hObject.Data(rowidx(idx),1);
                centroid=DATA.probe(probeidx).cluster(clusteridx).centroid_synapse;
                plot3(handles.PANEL_CLUSTER,centroid(:,1),centroid(:,2),centroid(:,3),...
                    'LineStyle','none','LineWidth',4,...
                    'MarkerSize',5,'Marker','+',...
                    'Color','w');
            end
        end
    end
end
handles.PANEL_CLUSTER.NextPlot=currenthold;

% --- Executes on button press in TOGGLE_HOLDPANELSCATTER3D.
function TOGGLE_HOLDPANELSCATTER3D_Callback(hObject, ~, handles)
if hObject.Value
    handles.PANEL_SCATTER3D.NextPlot='add';
else
    handles.PANEL_SCATTER3D.NextPlot='replacechild';
end

% --- Executes on button press in BUTTON_CLEARPANELSCATTER3D.
function BUTTON_CLEARPANELSCATTER3D_Callback(~, ~, handles)
cla(handles.PANEL_SCATTER3D);
handles.PANEL_SCATTER3D.Color=[0,0,0];
handles.PANEL_SCATTER3D.XColor=[0.7,0.7,0.7];
handles.PANEL_SCATTER3D.YColor=[0.7,0.7,0.7];
handles.PANEL_SCATTER3D.ZColor=[0.7,0.7,0.7];
handles.PANEL_SCATTER3D.GridColor=[0.5,0.5,0.5];
handles.PANEL_SCATTER3D.MinorGridColor=[0.5,0.5,0.5];
handles.PANEL_SCATTER3D.XGrid='on';
xlabel(handles.PANEL_SCATTER3D,'X');
handles.PANEL_SCATTER3D.YGrid='on';
ylabel(handles.PANEL_SCATTER3D,'Y');
handles.PANEL_SCATTER3D.ZGrid='on';
zlabel(handles.PANEL_SCATTER3D,'Z');
handles.PANEL_SCATTER3D.NextPlot='replacechild';
handles.TOGGLE_HOLDPANELSCATTER3D.Value=false;

% --- Executes on button press in TOGGLE_HOLDPANELCLUSTER.
function TOGGLE_HOLDPANELCLUSTER_Callback(hObject, ~, handles)
if hObject.Value
    handles.PANEL_CLUSTER.NextPlot='add';
else
    handles.PANEL_CLUSTER.NextPlot='replacechild';
end

% --- Executes on button press in TOGGLE_LINKAXES.
function TOGGLE_LINKAXES_Callback(hObject, ~, handles)
if hObject.Value
    addtarget(hObject.UserData,handles.PANEL_CLUSTER);
else
    removetarget(hObject.UserData,handles.PANEL_CLUSTER);
end

% --- Executes on button press in BUTTON_CLEARPANELCLUSTER.
function BUTTON_CLEARPANELCLUSTER_Callback(~, ~, handles)
cla(handles.PANEL_CLUSTER);
handles.PANEL_CLUSTER.Color=[0,0,0];
handles.PANEL_CLUSTER.XColor=[0.7,0.7,0.7];
handles.PANEL_CLUSTER.YColor=[0.7,0.7,0.7];
handles.PANEL_CLUSTER.ZColor=[0.7,0.7,0.7];
handles.PANEL_CLUSTER.GridColor=[0.5,0.5,0.5];
handles.PANEL_CLUSTER.MinorGridColor=[0.5,0.5,0.5];
handles.PANEL_CLUSTER.XGrid='on';
xlabel(handles.PANEL_CLUSTER,'X');
handles.PANEL_CLUSTER.YGrid='on';
ylabel(handles.PANEL_CLUSTER,'Y');
handles.PANEL_CLUSTER.ZGrid='on';
zlabel(handles.PANEL_CLUSTER,'Z');
handles.PANEL_CLUSTER.NextPlot='replacechild';
handles.TOGGLE_HOLDPANELCLUSTER.Value=false;

% --- Executes on selection change in MENU_SORTCLUSTERTABLE.
function MENU_SORTCLUSTERTABLE_Callback(hObject, ~, handles)
info=handles.TABLE_CLUSTERINFO.Data;
info=sortrows(info,hObject.Value);%sort by volume
handles.TABLE_CLUSTERINFO.Data=info;

% --- Executes on button press in BUTTON_STORECLUSTER.
function BUTTON_STORECLUSTER_Callback(~, ~, handles)
%save selected cluster centroid position so one could track them over files
global DATA;
probeidx=handles.MENU_PROBE.Value;
selectedrow=handles.TABLE_CLUSTERINFO.UserData;
if isempty(selectedrow)
    errordlg('No cluster selected in the current probe','Cluster find','modal');
else
    clusterindex=handles.TABLE_CLUSTERINFO.Data(:,1);
    storedcluster.index=num2cell(clusterindex(selectedrow));
    storedcluster.shape={DATA.probe(probeidx).cluster(clusterindex(selectedrow)).shape};
    storedcluster.centroid={DATA.probe(probeidx).cluster(clusterindex(selectedrow)).centroid};
    storedcluster.centroid_synapse={DATA.probe(probeidx).cluster(clusterindex(selectedrow)).centroid_synapse};
    if isfield(DATA.probe(probeidx).cluster(clusterindex(selectedrow(1))),'nnc0_id')
        storedcluster.nnc0_id={DATA.probe(probeidx).cluster(clusterindex(selectedrow)).nnc0_id};
    end
    if isfield(DATA.probe(probeidx).cluster(clusterindex(selectedrow(1))),'nnc1_id')
        storedcluster.nnc1_id={DATA.probe(probeidx).cluster(clusterindex(selectedrow)).nnc1_id}; %#ok<STRNU>
    end
    filename=cat(2,pwd,filesep,'storedcluster.mat');
    version='7.3';
    save(filename,'storedcluster','-mat',cat(2,'-v',version));
    msgbox(cat(2,'Cluster: ',sprintf('%d, ',clusterindex(selectedrow)),' stored.'),'Save File','modal');
end

% --- Executes on button press in BUTTON_FINDCLUSTER.
function BUTTON_FINDCLUSTER_Callback(~, ~, handles)
%load stored cluster and locate current cluster in the same locality
global DATA;
filename=cat(2,pwd,filesep,'storedcluster.mat');
temp=load(filename,'-mat');
storedcluster=temp.storedcluster;
if isempty(DATA.probe(1).cluster)
    errordlg('No cluster in the current probe','Cluster find','modal');
else
    % find clusters
    current_centroid=cell2mat({DATA.probe(1).cluster.centroid}');
    nstoredcluster=numel(storedcluster.centroid);
    current_clusteridx=zeros(nstoredcluster,1);
    rowidx=zeros(nstoredcluster,1);
    tableidx=handles.TABLE_CLUSTERINFO.Data(:,1);
    for clusteridx=1:nstoredcluster
        % find closes cluster to stored cluster location
        [~,current_clusteridx(clusteridx)]=min(sum(bsxfun(@minus,current_centroid, storedcluster.centroid{clusteridx}).^2,2));
        rowidx(clusteridx)=find(tableidx==current_clusteridx(clusteridx));
    end
    handles.TABLE_CLUSTERINFO.UserData=rowidx;
    %matchidx=[storedcluster.index';num2cell(current_clusteridx)'];
    %msgbox(cat(2,'Most likely cluster pair found are: ',sprintf('\n%5d -> %5d',matchidx{:})),'locate stored cluster');
    matchidx=cellfun(@(x,y)sprintf('%5d ---> %5d',x,y),storedcluster.index,num2cell(current_clusteridx),'UniformOutput',false);
    set(0,'DefaultUicontrolBackgroundColor',[0.3,0.3,0.3]);
    set(0,'DefaultUicontrolForegroundColor','k');
    hfig=figure('NumberTitle','off','Name','Locate Stored Cluster','MenuBar','none','Toolbar','none','Position',[10,10,200,500],'Color','k','Resize','off');
    hlist=uicontrol(hfig,'Style','listbox','String',matchidx,'Position',[10,10,180,480],'BackgroundColor','k','ForegroundColor','w','Enable','on','Max',2);
    hlist.Value=[];
end

% --- Executes on button press in BUTTON_MERGECLUSTER.
function BUTTON_MERGECLUSTER_Callback(hObject, eventdata, handles)
% hObject    handle to BUTTON_MERGECLUSTER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes when user attempts to close LOCSUPRES.
function LOCSUPRES_CloseRequestFcn(hObject, ~, ~)
% ask for confirmation
button = questdlg({'Are you sure you want to quit?',...
    'Make sure you have saved your work'},...
    'Check and Confirm','Yes','No','No');
switch button
    case 'Yes'
        % close main GUI
        delete(hObject);
    case 'No'
        %cancel closure
        msgbox(sprintf('%s\n','Return to work'),'Cancel closing','modal');
end

%==User functions=========================================================
function initialise(handles)
global DATA DEFAULT_COLOUR;
% set default colour scheme to black background and white font for dark
% room usage
set(0,'DefaultUicontrolBackgroundColor','k');
set(0,'DefaultUicontrolForegroundColor','w');
feature('accel','on');
% clear command window
clc;
% stop for debugging if error
dbstop if error;
% But don't bother with warnings
warning off all;
% get current file path
funcpath=mfilename('fullpath');
% move to base directory as we know
cd(funcpath(1:end-9));
% find all subdirectory
path=cat(2,pwd,filesep);
% add all subdirectory for libraries
addpath(genpath(path));
DEFAULT_COLOUR='brgymcwk';
% set up global variable DATA
DATA.datainfo.n_image=[];
DATA.datainfo.n_image=[];
DATA.datainfo.n_cycle=[];
DATA.datainfo.n_zstep=[];
DATA.datainfo.n_frame=[];
DATA.datainfo.n_accum=[];
DATA.datainfo.n_valid=[];
% data value parameters
DATA.datainfo.dp=1;
DATA.datainfo.dX=0.02;
DATA.datainfo.dY=0.02;
DATA.datainfo.dZ=0.05;
DATA.datainfo.dT=1;
DATA.datainfo.p=0;
DATA.datainfo.X=(0:DATA.datainfo.dX:20)'; % 20um default super rest fov
DATA.datainfo.Y=(0:DATA.datainfo.dY:20)'; % 20um default super rest fov
DATA.datainfo.Z=(-1:DATA.datainfo.dZ:1)'; % 2um default super rest fov
DATA.datainfo.T=1;
p_size=numel(DATA.datainfo.p);
X_size=numel(DATA.datainfo.X);
Y_size=numel(DATA.datainfo.Y);
Z_size=numel(DATA.datainfo.Z);
T_size=numel(DATA.datainfo.T);
DATA.datainfo.data_dim=[p_size,X_size,Y_size,Z_size,T_size];
DATA.datainfo.dataname=[];
DATA.datainfo.scaled=false;
DATA.datainfo.probe_colours='r';
% calculation parameters
DATA.datainfo.localdist=[1 1 1];
DATA.datainfo.min_cluster_size=[50 50 50];
DATA.datainfo.synaptic_range=0.2;
DATA.datainfo.per_synapse=0.5;
DATA.datainfo.slice_int=[5,5,5];
DATA.datainfo.isoval=[3,3,3];
DATA.datainfo.accum_threshold=1;
DATA.datainfo.psfx_threshold=0.125;
DATA.datainfo.psfy_threshold=0.125;
DATA.datainfo.psfz_threshold=0.5;
DATA.datainfo.snr_threshold=1.5;
DATA.datainfo.chisq_threshold=0.05;
DATA.datainfo.loglike_threshold=0;
DATA.datainfo.accuracy_threshold=10;
DATA.datainfo.last_change=datestr(now);
DATA.dataval=[];
DATA.probe(1).cluster=[];
% setup GUI components
handles.TABLE_DATAINFO.Data=[];
handles.TABLE_CLUSTERINFO.Data=[];
handles.TOGGLE_HOLDPANELSCATTER3D.Value=false;
handles.TOGGLE_HOLDPANELCLUSTER.Value=false;
handles.MENU_ANALYSOR.String={'reduce_dubious_localisation',...
    'make_cluster',...
    'analyse_cluster_neighbour',...
    'analyse_site_neighbour',...
    'intercluster_distribution',...
    'identify_synapse',...
    'analyse_synapse_nn_site',...
    'analyse_synapse_subcluster_site'};
handles.MENU_ANALYSOR.Value=1;
handles.MENU_.Value=1;
% SET PANEL_SCATTER3D PARAMETERS
BUTTON_CLEARPANELSCATTER3D_Callback([],[],handles);
% SET PANEL_CLUSTER PARAMETERS
BUTTON_CLEARPANELCLUSTER_Callback([],[],handles);
handles.TOGGLE_LINKAXES.UserData=linkprop([handles.PANEL_SCATTER3D,handles.PANEL_CLUSTER],...
    {'CameraPosition','CameraUpVector','XLim','YLim','ZLim'});
removetarget(handles.TOGGLE_LINKAXES.UserData,handles.PANEL_CLUSTER);

function display_datainfo( output_to )
%DISP_DATAINFO display file indata_formation in tables or text field
%   specify which data data_format and the file index
global DATA;
f_name=fieldnames(DATA.datainfo);
%content=cell(length(f_name),2);%add two base field and minus 5 dims
f_idx=0;o_idx=1;
while f_idx<length(f_name)
    f_idx=f_idx+1;
    switch f_name{f_idx}
        %not to display
        case {'p','X','Y','Z','T',...
                'header1','header2','header3','header4','header5','header6','header7','header8','header9','header10',...
                'header11','header12','header13','header14','header15','header16','header17','header18','header19','header20',...
                'header21','header22','header23','header24','header25','header26','header27','header28','header29','header30',...
                'header31','header32','header33','header34','header35','header36','header37','header38','header39','header40',...
                'imageidcol','cyclecol','zstepcol','framecol','accumcol','probecol','psfpccol','psfxcol','psfycol','psfzcol',...
                'xcol','ycol','zcol','bg11col','bg12col','bg21col','bg22col','chisqcol','loglhcol','accuracycol','validcol'}
        otherwise
            f_val=DATA.datainfo.(f_name{f_idx});%field_value
            if isnumeric(f_val)
                content{o_idx,1}=f_name{f_idx};%field name
                if numel(f_val)>10
                    content{o_idx,2}=sprintf('matrix of size %d x %d',size(f_val,1),size(f_val,2));
                else
                    content{o_idx,2}=f_val;
                end
                o_idx=o_idx+1;
            else
                if islogical(f_val)
                    content{o_idx,1}=f_name{f_idx};%field name
                    content{o_idx,2}=f_val;
                else
                    content{o_idx,1}=f_name{f_idx};%field name
                    if ishandle(f_val)
                        % problem with matlab root object
                        if f_val~=0
                            content{o_idx,2}=sprintf('handle of %s type',f_val.Type);
                        else
                            content{o_idx,2}=f_val;
                        end
                    else
                        if iscell(f_val)
                            content{o_idx,2}=char(f_val)';
                        else
                            if isa(f_val,'matlab.graphics.axis.Axes')
                                if f_val.isvalid
                                    content{o_idx,2}=f_val.Tag;
                                else
                                    content{o_idx,2}=[];
                                end
                            else
                                if isstruct(f_val)
                                    content{o_idx,2}=sprintf('structured data');
                                else
                                    content{o_idx,2}=f_val;
                                end
                            end
                        end
                    end
                end
                o_idx=o_idx+1;
            end
    end
end
info=cellfun(@(x)num2str(x),content,'UniformOutput',false);
if isempty(output_to)
    %default output to pipe
    cellfun(@(x)fprintf(1,'%s\n',x),info');
else
    output_to.Data=info;
end

function showlocation(handles)
global DATA;
% get probe index
probeidx=handles.MENU_PROBE.Value;
% plot
plot3(handles.PANEL_SCATTER3D,DATA.probe(probeidx).location(:,1),DATA.probe(probeidx).location(:,2),DATA.probe(probeidx).location(:,3),...
    'LineStyle','none','Marker','.','Color',DATA.datainfo.probe_colours(probeidx),...
    'MarkerSize',1);
xlim(handles.PANEL_SCATTER3D,[DATA.datainfo.X(1),DATA.datainfo.X(end)]);xlabel('X');
ylim(handles.PANEL_SCATTER3D,[DATA.datainfo.Y(1),DATA.datainfo.Y(end)]);ylabel('Y');
zlim(handles.PANEL_SCATTER3D,[DATA.datainfo.Z(1),DATA.datainfo.Z(end)]);zlabel('Z');
view(handles.PANEL_SCATTER3D,2);
daspect(handles.PANEL_SCATTER3D,[1 1 1]);

function showimage(handles)
global DATA;
% get probe index
probeidx=handles.MENU_PROBE.Value;
X_scale=DATA.datainfo.X;
X_res=DATA.datainfo.dX;
Y_scale=DATA.datainfo.Y;
Y_res=DATA.datainfo.dY;
Z_scale=DATA.datainfo.Z;
Z_res=DATA.datainfo.dZ;
val=zeros(DATA.datainfo.data_dim(2),DATA.datainfo.data_dim(3),DATA.datainfo.data_dim(4));
for zslice=1:DATA.datainfo.data_dim(4)
    % go through each z slices
    currentzidx=find(DATA.probe(probeidx).location(:,3)>(Z_scale(zslice)-Z_res/2)&DATA.probe(probeidx).location(:,3)<=(Z_scale(zslice)+Z_res/2));
    [temp,~,~,binx,~]=histcounts2(DATA.probe(probeidx).location(currentzidx,1),DATA.probe(probeidx).location(currentzidx,2),[X_scale-X_res/2;X_scale(end)+X_res/2],[Y_scale-Y_res/2;Y_scale(end)+Y_res/2]);
    if DATA.datainfo.scaled
        for idx=1:numel(currentzidx)
            if binx(idx)>0
                %scaling using column 15
            end
        end
    end
    val(:,:,zslice)=temp;
end
DATA.datainfo.slice_int=[5,5,5];
DATA.datainfo.isoval=[3,3,3];
binsize=DATA.datainfo.slice_int;
% calculate binning
dim_size=DATA.datainfo.data_dim(2:4);
% work out new data size
binsize(binsize>dim_size)=dim_size(binsize>dim_size); % make sure bin<=dim
newsize=floor(dim_size./binsize);% can only have full number bins
offset_size=newsize.*binsize;
%binning
temp=[binsize;newsize];
reshapesize=temp(:);
temp=val(1:offset_size(1),1:offset_size(2),1:offset_size(3));
temp=reshape(temp,reshapesize');
temp=mean(mean(mean(temp,1),3),5);
val=reshape(temp,newsize);
X_scale=mean(reshape(DATA.datainfo.X(1:offset_size(1)),[binsize(1),newsize(1)]),1);
Y_scale=mean(reshape(DATA.datainfo.Y(1:offset_size(2)),[binsize(2),newsize(2)]),1);
Z_scale=mean(reshape(DATA.datainfo.Z(1:offset_size(3)),[binsize(3),newsize(3)]),1);
[y,x,z]=meshgrid(X_scale,Y_scale,Z_scale);
val=smooth3(val,'gaussian',[3,3,1],0.6);
% plot
p = patch(handles.PANEL_SCATTER3D,isosurface(x,y,z,(val>DATA.datainfo.isoval(probeidx)*mean(val(:))),0.9,'noshare'));
isonormals(y,x,z,val,p);
p.FaceColor = DATA.datainfo.probe_colours(probeidx);
p.EdgeColor = 'none';
p.FaceAlpha=0.3;
view(handles.PANEL_SCATTER3D,2);
daspect(handles.PANEL_SCATTER3D,[1 1 1]);

function display_clusterinfo( probeidx, handles )
global DATA;
output_to=handles.TABLE_CLUSTERINFO;
n_cluster=numel(DATA.probe(probeidx).cluster);
if n_cluster>0
    colname=fieldnames(DATA.probe(probeidx).cluster);
    info(:,1)=(1:1:n_cluster)';
    colidx=2;
    for fidx=1:numel(colname)
        switch colname{fidx}
            case 'shape'
                info(:,colidx)=cell2mat(arrayfun(@(x)numel(x.index),DATA.probe(probeidx).cluster,'UniformOutput',false))';
                colname{fidx}='n_site';
                colidx=colidx+1;
            case 'area'
                info(:,colidx)=[DATA.probe(probeidx).cluster.area]';
                colidx=colidx+1;
            case 'volume'
                info(:,colidx)=[DATA.probe(probeidx).cluster.volume]';
                colidx=colidx+1;
            case 'centroid'
                info(:,colidx:colidx+2)=cell2mat({DATA.probe(probeidx).cluster.centroid}');
                colidx=colidx+3;
            case 'Dist'
                info(:,colidx)=cellfun(@(x)max(x),{DATA.probe(probeidx).cluster.Dist})';
                info(:,colidx+1)=cellfun(@(x)mean(x),{DATA.probe(probeidx).cluster.Dist})';
                info(:,colidx+2)=cellfun(@(x)median(x),{DATA.probe(probeidx).cluster.Dist})';
                colidx=colidx+3;
            case {'nnc0_id','nnc1_id','nnc2_id','nnc3_id','nnc0_dist','nnc1_dist','nnc2_dist','nnc3_dist'}
                info(:,colidx)=[DATA.probe(probeidx).cluster.(colname{fidx})]';
                colidx=colidx+1;
            case {'nns0_id','nns1_id','nns2_id','nns3_id','nns0_dist','nns1_dist','nns2_dist','nns3_dist'}
                info(:,colidx)=[DATA.probe(probeidx).cluster.(colname{fidx})]';
                colidx=colidx+1;
            case {'nns0_pt_rmin','nns0_pt_rmean','nns0_pt_rmedian','nns0_pt_Vf','nns0_pt_density','nns0_pt_Asurf',...
                    'nns1_pt_rmin','nns1_pt_rmean','nns1_pt_rmedian','nns1_pt_Vf','nns1_pt_density','nns1_pt_Asurf',...
                    'nns2_pt_rmin','nns2_pt_rmean','nns2_pt_rmedian','nns2_pt_Vf','nns2_pt_density','nns2_pt_Asurf',...
                    'nns0_sh_rmin','nns0_sh_rmean','nns0_sh_rmedian','nns0_sh_Vf','nns0_sh_density','nns0_sh_Asurf',...
                    'nns1_sh_rmin','nns1_sh_rmean','nns1_sh_rmedian','nns1_sh_Vf','nns1_sh_density','nns1_sh_Asurf',...
                    'nns2_sh_rmin','nns2_sh_rmean','nns2_sh_rmedian','nns2_sh_Vf','nns2_sh_density','nns2_sh_Asurf'}
                marker=cell2mat(cellfun(@(x)~isempty(x),{DATA.probe(probeidx).cluster.(colname{fidx})},'UniformOutput',false));
                info(marker,colidx)=[DATA.probe(probeidx).cluster.(colname{fidx})]';
                colidx=colidx+1;
            case 'centroid_synapse'
                info(:,colidx:colidx+2)=cell2mat({DATA.probe(probeidx).cluster.centroid_synapse}');
                colidx=colidx+3;
            case 'loc_synapse'
                %info(:,colidx:colidx+2)=cell2mat({DATA.probe(probeidx).cluster.loc_synapse}');
                %colidx=colidx+3;
            case 'mean_dist_synapse'
                info(:,colidx:colidx+2)=cell2mat({DATA.probe(probeidx).cluster.mean_dist_synapse}');
                colidx=colidx+3;
            case {'synapse','view'}
            otherwise
                
        end
    end
    % add centroid colnames
    temp=find(strcmp(colname,'centroid'));
    if ~isempty(temp)
        colname(temp+3:end+2)=colname(temp+1:end);
        colname(temp:temp+2)={'C(x)','C(y)','C(z)'};
    end
    temp=find(strcmp(colname,'centroid_synapse'));
    if ~isempty(temp)
        colname(temp+3:end+2)=colname(temp+1:end);
        colname(temp:temp+2)={'C_synapse(x)','C_synapse(y)','C_synapse(z)'};
    end
    % add Dist colnames
    temp=find(strcmp(colname,'Dist'));
    if ~isempty(temp)
        colname(temp+3:end+2)=colname(temp+1:end);
        colname(temp:temp+2)={'Mean Dist','Median Dist','Max Dist'};
    end
    % add Dist colnames
    temp=find(strcmp(colname,'mean_dist_synapse'));
    if ~isempty(temp)
        colname(temp+3:end+2)=colname(temp+1:end);
        colname(temp:temp+2)={'Mean Syn Dist','Median Syn Dist','Max Syn Dist'};
    end
    % remove loc_synapse (deprecated fields)
    temp=find(strcmp(colname,'loc_synapse'));
    if ~isempty(temp)
        %colname(temp+3:end+2)=colname(temp+1:end);
        %colname(temp:temp+2)={'loc_synapse(x)','loc_synapse(y)','loc_synapse(z)'};
        colname(temp)=[];
    end
    % remove structured field synapse
    temp=find(strcmp(colname,'synapse'));
    if ~isempty(temp)
        colname(temp)=[];
    end
    % remove view field
    temp=find(strcmp(colname,'view'));
    if ~isempty(temp)
        colname(temp)=[];
    end
else
    info=[];
    colname=[];
end
output_to.Data=info;
output_to.ColumnName=colname;
handles.MENU_SORTCLUSTERTABLE.String=colname;
handles.MENU_SORTCLUSTERTABLE.Value=1;

function showcluster( probeidx, cidx, handles )
global DATA;
if isempty(handles)
    plot(DATA.probe(probeidx).cluster(cidx).shape,'Facecolor',DATA.datainfo.probe_colours(probeidx),...
        'EdgeColor','w',...
        'EdgeAlpha',0.2,...
        'FaceAlpha',0.5)
    view(gca,2);
    daspect(gca,[1 1 1]);
    xlim('auto');ylim('auto');zlim('auto');
else
    plot(DATA.probe(probeidx).cluster(cidx).shape,'Facecolor',DATA.datainfo.probe_colours(probeidx),...
        'EdgeColor','w',...
        'EdgeAlpha',0.2,...
        'FaceAlpha',0.5,...
        'Parent',handles.PANEL_CLUSTER);
    view(handles.PANEL_CLUSTER,2);
    daspect(handles.PANEL_CLUSTER,[1 1 1]);
    xlim('auto');ylim('auto');zlim('auto');
end

function figure_keypress(~,eventkey)
switch eventkey.Key
    case {'f9'}
        % export trace
        %export_panel(findobj(handle,'Type','Axes'));
        export_panel(gca);
    case {'f1','f2','f3','f4'}
        % plot scatter for synapse analysis subplot
        hplot=gca;
        switch hplot.Tag
            case {'synapse_cluster','synapse_ncluster'}
                switch hplot.Children(str2double(eventkey.Key(2))).Visible
                    case 'off'
                        hplot.Children(str2double(eventkey.Key(2))).Visible='on';
                    case 'on'
                        hplot.Children(str2double(eventkey.Key(2))).Visible='off';
                end
        end
    case {'f11'}
        % link axes
        hplot=gca;
        switch hplot.Tag
            case 'synapse_ncluster'
                temp=findobj(hplot.Parent,'Type','Axes');% find all subplot
                hplotidx=find(cellfun(@(x)isempty(x.String),{temp.Title}));% find all scatter
                
                linkprop(temp(fliplr(hplotidx)),{'CameraPosition','CameraUpVector','XLim','YLim','ZLim'});
        end
end

%==User Functions for analysor============================================
function reduce_dubious_localisation(handles)
global DATA;
probeidx=handles.MENU_PROBE.Value;
probe=find(DATA.dataval(:,DATA.datainfo.probecol)==DATA.datainfo.p(probeidx));
fh=figure('Name',sprintf('Parameter Distribution for %s',DATA.datainfo.(cat(2,'probe',num2str(DATA.datainfo.p(probeidx)),'_name'))),...
    'NumberTitle','off',...
    'MenuBar','none',...
    'ToolBar','figure',...
    'Position',[0,0,900,600],...
    'Keypressfcn',@figure_keypress); %#ok<NASGU>
invalid_accum=DATA.dataval(probe,DATA.datainfo.accumcol)>DATA.datainfo.accum_threshold;%accumulation
subplot(2,4,1);
histogram(DATA.dataval(probe,DATA.datainfo.accumcol),linspace(0,5,10));
hold on;line([DATA.datainfo.accum_threshold,DATA.datainfo.accum_threshold],get(gca,'YLim'),'LineStyle','--','LineWidth',3,'Color','r');title('accum');hold off;
subplot(2,4,2);
histogram(DATA.dataval(probe,DATA.datainfo.psfxcol),linspace(-0.5,0.5,50));
hold on;line([DATA.datainfo.psfx_threshold,DATA.datainfo.psfx_threshold],get(gca,'YLim'),'LineStyle','--','LineWidth',3,'Color','r');
line(-[DATA.datainfo.psfx_threshold,DATA.datainfo.psfx_threshold],get(gca,'YLim'),'LineStyle','--','LineWidth',3,'Color','r');title('psfx');hold off;
subplot(2,4,3);histogram(DATA.dataval(probe,DATA.datainfo.psfycol),linspace(-0.5,0.5,50));
hold on;line([DATA.datainfo.psfy_threshold,DATA.datainfo.psfy_threshold],get(gca,'YLim'),'LineStyle','--','LineWidth',3,'Color','r');
line(-[DATA.datainfo.psfy_threshold,DATA.datainfo.psfy_threshold],get(gca,'YLim'),'LineStyle','--','LineWidth',3,'Color','r');title('psfy');hold off;
subplot(2,4,4);histogram(DATA.dataval(probe,DATA.datainfo.psfzcol),linspace(-1,1,50));
hold on;line([DATA.datainfo.psfz_threshold,DATA.datainfo.psfz_threshold],get(gca,'YLim'),'LineStyle','--','LineWidth',3,'Color','r');
line(-[DATA.datainfo.psfz_threshold,DATA.datainfo.psfz_threshold],get(gca,'YLim'),'LineStyle','--','LineWidth',3,'Color','r');title('psfz');hold off;
invalid_x=abs(DATA.dataval(probe,DATA.datainfo.psfxcol))>DATA.datainfo.psfx_threshold;%250nm good psfx
invalid_y=abs(DATA.dataval(probe,DATA.datainfo.psfycol))>DATA.datainfo.psfy_threshold;%250nm good psfy
invalid_z=abs(DATA.dataval(probe,DATA.datainfo.psfzcol))>DATA.datainfo.psfz_threshold;%1050nm good psfz
sn_ratio=DATA.dataval(probe,DATA.datainfo.psfpccol)./sum(DATA.dataval(probe,[DATA.datainfo.bg11col,DATA.datainfo.bg12col,DATA.datainfo.bg21col,DATA.datainfo.bg22col]),2);%psf-photon-count/background11
subplot(2,4,5);
histogram(sn_ratio,linspace(min(sn_ratio),max(sn_ratio),50));
hold on;line([DATA.datainfo.snr_threshold,DATA.datainfo.snr_threshold],get(gca,'YLim'),'LineStyle','--','LineWidth',3,'Color','r');title('snr');hold off;
invalid_snr=sn_ratio<DATA.datainfo.snr_threshold;
chisq=DATA.dataval(probe,DATA.datainfo.chisqcol)./DATA.dataval(probe,DATA.datainfo.psfpccol);
subplot(2,4,6);
histogram(chisq,linspace(min(chisq),max(chisq),50));
hold on;line([DATA.datainfo.chisq_threshold,DATA.datainfo.chisq_threshold],get(gca,'YLim'),'LineStyle','--','LineWidth',3,'Color','r');title('chisq');hold off;
invalid_chisq=chisq>=DATA.datainfo.chisq_threshold;
loglike=DATA.dataval(probe,DATA.datainfo.loglhcol)./DATA.dataval(probe,DATA.datainfo.psfpccol);
subplot(2,4,7);histogram(loglike,linspace(-2,0,50));
hold on;line([DATA.datainfo.loglike_threshold,DATA.datainfo.loglike_threshold],get(gca,'YLim'),'LineStyle','--','LineWidth',3,'Color','r');title('loglike');hold off;
invalid_loglike=loglike>DATA.datainfo.loglike_threshold;
accuracy=DATA.dataval(probe,DATA.datainfo.accuracycol);
subplot(2,4,8);
histogram(accuracy,linspace(min(accuracy),max(accuracy),50));
hold on;line([DATA.datainfo.accuracy_threshold,DATA.datainfo.accuracy_threshold],get(gca,'YLim'),'LineStyle','--','LineWidth',3,'Color','r');title('accuracy');hold off;
invalid_accuracy=accuracy<=DATA.datainfo.accuracy_threshold;
% combine all
invalid=invalid_accum|invalid_x|invalid_y|invalid_z|invalid_snr|invalid_chisq|invalid_loglike|invalid_accuracy;
fh=figure('Name','Comparison of clean up operation',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'ToolBar','figure',...
    'Position',[0,0,900,600],...
    'Keypressfcn',@figure_keypress);
% plot scatter to compare
subplot(1,2,1);set(gca,'Color',[0,0,0]);
plot3(DATA.dataval(probe(~invalid),DATA.datainfo.xcol),DATA.dataval(probe(~invalid),DATA.datainfo.ycol),DATA.dataval(probe(~invalid),DATA.datainfo.zcol),...
    'LineStyle','none','Marker','.','Color',DATA.datainfo.probe_colours(probeidx),...
    'MarkerSize',1);
set(gca,'Color',[0,0,0]);
title('Cleaned (close figure to chose)');
xlim([DATA.datainfo.X(1),DATA.datainfo.X(end)]);xlabel('X');
ylim([DATA.datainfo.Y(1),DATA.datainfo.Y(end)]);ylabel('Y');
zlim([DATA.datainfo.Z(1),DATA.datainfo.Z(end)]);zlabel('Z');
view(gca,2);
daspect(gca,[1 1 1]);
subplot(1,2,2);
plot3(DATA.dataval(probe,DATA.datainfo.xcol),DATA.dataval(probe,DATA.datainfo.ycol),DATA.dataval(probe,DATA.datainfo.zcol),...
    'LineStyle','none','Marker','.','Color',DATA.datainfo.probe_colours(probeidx),...
    'MarkerSize',1);
set(gca,'Color',[0,0,0]);
title('Original');
xlim([DATA.datainfo.X(1),DATA.datainfo.X(end)]);xlabel('X');
ylim([DATA.datainfo.Y(1),DATA.datainfo.Y(end)]);ylabel('Y');
zlim([DATA.datainfo.Z(1),DATA.datainfo.Z(end)]);zlabel('Z');
view(gca,2);
daspect(gca,[1 1 1]);
uiwait(fh);
% ask if you want to do this
answer=questdlg('Do you want to apply this operation? It is irreversible!!','reduce_dubious_localisation','Yes','No','Reset','Yes');
switch answer
    case 'Yes'
        DATA.probe(probeidx).location=DATA.dataval(probe(~invalid),[DATA.datainfo.xcol,DATA.datainfo.ycol,DATA.datainfo.zcol]);
    case 'Reset'
        DATA.probe(probeidx).location=DATA.dataval(probe,[DATA.datainfo.xcol,DATA.datainfo.ycol,DATA.datainfo.zcol]);
    case 'No'
        
end

function make_cluster(handles)
global DATA;
try
    % make clusters using clusterdata function
    pidx=handles.MENU_PROBE.Value;
    %probe=find(DATA.dataval(:,DATA.datainfo.probecol)==DATA.datainfo.p(pidx));
    %position=DATA.dataval(probe,[DATA.datainfo.xcol,DATA.datainfo.ycol,DATA.datainfo.zcol]);
    position=DATA.probe(pidx).location;
    if size(position,1)<8e4
        savemem='off';
    else
        savemem='on';
    end
    %create waitbar for user attention
    waitbar_handle = waitbar(0,'Please wait...','Progress Bar','Calculating...',...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)',...
        'WindowStyle','normal',...
        'Color',[0.2,0.2,0.2]);
    setappdata(waitbar_handle,'canceling',0);
    T = clusterdata(position,'criterion','distance','cutoff',DATA.datainfo.localdist(pidx)/2,'distance','euclidean','linkage','median','savememory',savemem);
    [clustersize,~,~]=histcounts(T,0.5:1:max(T)+0.5);
    validclusterid=find(clustersize>DATA.datainfo.min_cluster_size(pidx));
    xlim(handles.PANEL_CLUSTER,[DATA.datainfo.X(1),DATA.datainfo.X(end)]);xlabel('X');
    ylim(handles.PANEL_CLUSTER,[DATA.datainfo.Y(1),DATA.datainfo.Y(end)]);ylabel('Y');
    zlim(handles.PANEL_CLUSTER,[DATA.datainfo.Z(1),DATA.datainfo.Z(end)]);zlabel('Z');
    n_cluster=numel(validclusterid);
    DATA.probe(pidx).cluster=[];
    for clusteridx=1:n_cluster
        %output some progress so we know it is doing things
        if getappdata(waitbar_handle,'canceling')
            fprintf('make cluster calculation cancelled\n');
            delete(waitbar_handle);       % DELETE the waitbar; don't try to CLOSE it.
            return;
        end
        clusterid=validclusterid(clusteridx);
        particle_idx=(T==clusterid);
        DATA.probe(pidx).cluster(clusteridx).index=find(particle_idx);
        P_cluster=[position(particle_idx,1),position(particle_idx,2),position(particle_idx,3)];
        %k=boundary(P_cluster,1);
        %DATA.probe(pidx).cluster(clusteridx).surface=k;
        temp=alphaShape(P_cluster,1,'HoleThreshold',0.2);
        pc=criticalAlpha(temp,'one-region');
        temp.Alpha=pc;
        DATA.probe(pidx).cluster(clusteridx).shape=temp;
        DATA.probe(pidx).cluster(clusteridx).volume=volume(temp);
        DATA.probe(pidx).cluster(clusteridx).area=surfaceArea(temp);
        [~,C,~,D] = kmeans(P_cluster,1);
        DATA.probe(pidx).cluster(clusteridx).centroid=C;
        DATA.probe(pidx).cluster(clusteridx).Dist=D;
        %showcluster(pidx, clusteridx, handles);
        
        % Report current estimate in the waitbar's message field
        done=clusteridx/n_cluster;
        waitbar(done,waitbar_handle,sprintf('%3.1f%%',100*done));
    end
    delete(waitbar_handle);       % DELETE the waitbar; don't try to CLOSE it.
    view(handles.PANEL_CLUSTER,2);
    daspect(handles.PANEL_CLUSTER,[1 1 1]);
    % sort by centroid
    [~,ind]=sortrows(cell2mat({DATA.probe(pidx).cluster.centroid}'),[1,2,3]);
    DATA.probe(pidx).cluster=DATA.probe(pidx).cluster(ind);
    % update cluster info table
    display_clusterinfo(pidx,handles);
    msgbox(sprintf('%s clustering successfully analysed.\nIt has %g clusters.\nUse localdist and min_cluster_size are variable parameters.',DATA.datainfo.(cat(2,'probe',num2str(DATA.datainfo.p(pidx)),'_name')),n_cluster),'Cluster Analysis','modal');
catch exception
    if exist('waitbar_handle','var')&&ishandle(waitbar_handle)
        delete(waitbar_handle);
    end
    errordlg(exception.message,'Calculation Error','modal');
end

function analyse_cluster_neighbour(handles)
% find pairing using centroid positions and minimum distance requirement
% using k nearest neighbour function
global DATA;
try
    % update probe menu
    fname=fieldnames(DATA.datainfo);
    fidx=find(cellfun(@(x)~isempty(x),regexp(fname,'probe\w_name')));
    probe_list=cell(DATA.datainfo.n_probe,1);
    for idx=1:DATA.datainfo.n_probe
        probe_list{idx}=DATA.datainfo.(fname{fidx(idx)});
    end
    set(0,'DefaultUicontrolBackgroundColor',[0.3,0.3,0.3]);
    set(0,'DefaultUicontrolForegroundColor','k');
    [s,v] = listdlg('PromptString','Select neighbouring probe:',...
        'SelectionMode','multiple',...
        'InitialValue',1:1:numel(probe_list),...
        'ListString',probe_list);
    set(0,'DefaultUicontrolBackgroundColor','k');
    set(0,'DefaultUicontrolForegroundColor','w');
    if v%validated by click on OK
        currentprobe=handles.MENU_PROBE.Value;
        currentcentroid=cell2mat({DATA.probe(currentprobe).cluster.centroid}');
        n_cluster=numel(DATA.probe(currentprobe).cluster);
        %create waitbar for user attention
        waitbar_handle = waitbar(0,'Please wait...','Progress Bar','Calculating...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)',...
            'WindowStyle','normal',...
            'Color',[0.2,0.2,0.2]);
        setappdata(waitbar_handle,'canceling',0);
        for probeidx=1:numel(s)%go through each probe including self if selected
            if ~isempty(DATA.probe(s(probeidx)).cluster)
                othercentroid=cell2mat({DATA.probe(s(probeidx)).cluster.centroid}');
                if currentprobe==s(probeidx)%nearest neighbour of same kind
                    [IDX,D] = knnsearch(othercentroid,currentcentroid,...
                        'K',2,'IncludeTies',false,'NSMethod','kdtree','Distance','euclidean');
                    IDX(:,1)=[];D(:,1)=[];
                else%nearest neighbour of a different kind
                    [IDX,D] = knnsearch(othercentroid,currentcentroid,...
                        'K',1,'IncludeTies',false,'NSMethod','kdtree','Distance','euclidean');
                end
                for clusterid=1:n_cluster
                    %output some progress so we know it is doing things
                    if getappdata(waitbar_handle,'canceling')
                        fprintf('search cluster neighbour calculation cancelled\n');
                        delete(waitbar_handle);       % DELETE the waitbar; don't try to CLOSE it.
                        return;
                    end
                    DATA.probe(currentprobe).cluster(clusterid).(cat(2,'nnc',num2str(DATA.datainfo.p(s(probeidx))),'_id'))=IDX(clusterid);
                    DATA.probe(currentprobe).cluster(clusterid).(cat(2,'nnc',num2str(DATA.datainfo.p(s(probeidx))),'_dist'))=D(clusterid);
                    % Report current estimate in the waitbar's message field
                    done=clusterid/n_cluster;
                    waitbar(done,waitbar_handle,sprintf('%3.1f%%',100*done));
                end
            end
        end
        delete(waitbar_handle);       % DELETE the waitbar; don't try to CLOSE it.
        display_clusterinfo(currentprobe,handles);
        msgbox(sprintf('%s nearest neighbour cluster search successfully analysed.\n',DATA.datainfo.(cat(2,'probe',num2str(DATA.datainfo.p(s(probeidx))),'_name'))),'Cluster Analysis','modal');
    else
        errordlg(sprintf('search neighbour cancelled\n'));
    end
catch exception
    if exist('waitbar_handle','var')&&ishandle(waitbar_handle)
        delete(waitbar_handle);
    end
    errordlg(exception.message,'Calculation Error','modal');
end

function analyse_site_neighbour(handles)
% find pairing using centroid positions and minimum distance requirement
% using k nearest neighbour function
global DATA;
try
    % update probe menu
    fname=fieldnames(DATA.datainfo);
    fidx=find(cellfun(@(x)~isempty(x),regexp(fname,'probe\w_name')));
    probe_list=cell(DATA.datainfo.n_probe,1);
    for idx=1:DATA.datainfo.n_probe
        probe_list{idx}=DATA.datainfo.(fname{fidx(idx)});
    end
    set(0,'DefaultUicontrolBackgroundColor',[0.3,0.3,0.3]);
    set(0,'DefaultUicontrolForegroundColor','k');
    [s,v] = listdlg('PromptString','Select neighbouring probe:',...
        'SelectionMode','multiple',...
        'InitialValue',1:1:numel(probe_list),...
        'ListString',probe_list);
    set(0,'DefaultUicontrolBackgroundColor','k');
    set(0,'DefaultUicontrolForegroundColor','w');
    if v%validated by click on OK
        plotcount=0;
        currentprobe=handles.MENU_PROBE.Value;
        selected_cluster=handles.TABLE_CLUSTERINFO.Data(handles.TABLE_CLUSTERINFO.UserData,1);
        currentcentroid=cell2mat({DATA.probe(currentprobe).cluster(selected_cluster).centroid}');
        n_cluster=numel(selected_cluster);
        prox_dist=max(DATA.datainfo.localdist);
        probeidx=1;
        for clusterid=1:n_cluster
            fh=figure('Name',sprintf('Nearest Probe Site Distance to probe %s cluster%g Centroid',probe_list{currentprobe},selected_cluster(clusterid)),...
                'NumberTitle','off',...
                'MenuBar','none',...
                'ToolBar','figure',...
                'Position',[0,0,900,600],...
                'Keypressfcn',@figure_keypress);
            for probeidx=1:numel(s)%go through each probe including self if selected
                othersite=DATA.probe(s(probeidx)).location;
                d_len=bsxfun(@minus,currentcentroid(clusterid,:),othersite);
                [az,el,rad]=cart2sph(d_len(:,1),d_len(:,2),d_len(:,3));
                proximity=rad<=prox_dist;
                rad=rad(proximity);
                az=rad2deg(az(proximity));el=rad2deg(el(proximity));
                subplot(2,numel(s),probeidx,'Parent',fh);
                histogram2(az,el,-180:5:180,-90:5:90,'FaceColor',DATA.datainfo.probe_colours(s(probeidx)));
                view([0,90]);axis('equal');
                xlabel('az (^o)','FontSize',8);
                ylabel('el (^o)','FontSize',8);
                title(sprintf('Probe %s',probe_list{probeidx}));
                subplot(2,numel(s),probeidx+numel(s),'Parent',fh);
                histogram(rad,linspace(0,prox_dist,25),'FaceColor',DATA.datainfo.probe_colours(s(probeidx)));
                xlabel('r (\mum)','FontSize',8);
                title(sprintf('median = %f',median(rad)));
                %radial_dist{probeidx,clusterid}=[rad(proximity),az(proximity),el(proximity)];
            end
            plotcount=plotcount+1;
        end
        %display_clusterinfo(currentprobe,handles);
        msgbox(sprintf('%s nearest neighbour site search successfully analysed.\n',probe_list{s(probeidx)}),'Cluster Analysis','modal');
    else
        errordlg(sprintf('search neighbour cancelled\n'));
    end
catch exception
    if exist('waitbar_handle','var')&&ishandle(waitbar_handle)
        delete(waitbar_handle);
    end
    errordlg(exception.message,'Calculation Error','modal');
end

function intercluster_distribution(handles)
% calculate the distribution of cluster sites to its centroid and the
% distribution of the shortest distance of the sites of neighbouring
% cluster to current cluster
global DATA;
try
    temp=handles.TABLE_CLUSTERINFO.Data;
    probeidx=handles.MENU_PROBE.Value;
    selectedrow=handles.TABLE_CLUSTERINFO.UserData;
    n_site=numel(selectedrow);
    fh=figure('Name','Nearest Neighbour Cluster Distance to Their Centroid',...
        'NumberTitle','off',...
        'MenuBar','none',...
        'ToolBar','figure',...
        'Position',[0,0,900,600],...
        'Keypressfcn',@figure_keypress);
    %create waitbar for user attention
    waitbar_handle = waitbar(0,'Please wait...','Progress Bar','Calculating...',...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)',...
        'WindowStyle','normal',...
        'Color',[0.2,0.2,0.2]);
    setappdata(waitbar_handle,'canceling',0);
    for ridx=1:n_site
        %output some progress so we know it is doing things
        if getappdata(waitbar_handle,'canceling')
            fprintf('make cluster calculation cancelled\n');
            delete(waitbar_handle);       % DELETE the waitbar; don't try to CLOSE it.
            return;
        end
        rowidx=selectedrow(ridx);
        pc_idx=temp(rowidx,1);
        pc_cluster=DATA.probe(probeidx).cluster(pc_idx);
        % distance distribution of individual clusters
        for pidx=1:DATA.datainfo.n_probe
            if isfield(pc_cluster,cat(2,'nnc',num2str(pidx-1),'_id'))
                if probeidx~=pidx
                    clusterid=pc_cluster.(cat(2,'nnc',num2str(pidx-1),'_id'));
                    po_cluster=DATA.probe(pidx).cluster(clusterid);
                    sh1=subplot(n_site,DATA.datainfo.n_probe*2,(ridx-1)*DATA.datainfo.n_probe*2+pidx,'Parent',fh);
                    hist(sh1,po_cluster.Dist,20);
                    sh1.YTickLabel=[];
                    pairedDist{pidx,ridx} = pdist2(DATA.probe(pidx).location(po_cluster.index,:),DATA.probe(probeidx).location(pc_cluster.index,:),'euclidean','Smallest',1); %#ok<AGROW>
                    sh2=subplot(n_site,DATA.datainfo.n_probe*2,(ridx-1)*DATA.datainfo.n_probe*2+pidx+DATA.datainfo.n_probe,'Parent',fh);
                    hist(sh2,pairedDist{pidx,ridx},20);
                    sh2.YTickLabel=[];
                    ylabel(sh1,sprintf('cluster%g',clusterid),'FontSize',8);
                    if ridx==n_site
                        nnc_dist=[DATA.probe(probeidx).cluster.(cat(2,'nnc',num2str(pidx-1),'_dist'))];
                        centroid_dist{pidx}=nnc_dist(temp(selectedrow,1))'; %#ok<AGROW>
                        xlabel(sh1,'r (\mum)','FontSize',8);
                        xlabel(sh2,'r (\mum)','FontSize',8);
                    end
                    if ridx==1
                        title(sh1,sprintf('nn %s cluster',eval(cat(2,'DATA.datainfo.probe',num2str(pidx-1),'_name'))),'FontSize',8);
                        title(sh2,sprintf('d_{shortest} to nn %s sites',eval(cat(2,'DATA.datainfo.probe',num2str(pidx-1),'_name'))),'FontSize',8);
                    end
                else
                    sh1=subplot(n_site,DATA.datainfo.n_probe*2,(ridx-1)*DATA.datainfo.n_probe*2+pidx,'Parent',fh);
                    hist(sh1,pc_cluster.Dist,20);
                    sh1.YTickLabel=[];
                    clusterid=pc_cluster.(cat(2,'nnc',num2str(pidx-1),'_id'));
                    po_cluster=DATA.probe(pidx).cluster(clusterid);
                    % only find distance to the next neighbour of the same kind
                    pairedDist{pidx,ridx} = pdist2(DATA.probe(pidx).location(po_cluster.index,:),DATA.probe(probeidx).location(pc_cluster.index,:),'euclidean','Smallest',1); %#ok<AGROW>
                    sh2=subplot(n_site,DATA.datainfo.n_probe*2,(ridx-1)*DATA.datainfo.n_probe*2+pidx+DATA.datainfo.n_probe,'Parent',fh);
                    hist(sh2,pairedDist{pidx,ridx},20);
                    sh2.YTickLabel=[];
                    ylabel(sh1,sprintf('cluster%g',pc_idx),'FontSize',8);
                    if ridx==n_site
                        %centroid_dist{pidx}=[DATA.probe(probeidx).cluster.(cat(2,'nnc',num2str(pidx-1),'_dist'))]'; %#ok<AGROW>
                        nnc_dist=[DATA.probe(probeidx).cluster.(cat(2,'nnc',num2str(pidx-1),'_dist'))];
                        centroid_dist{pidx}=nnc_dist(temp(selectedrow,1))'; %#ok<AGROW>
                        xlabel(sh1,'r (\mum)','FontSize',8);
                        xlabel(sh2,'r (\mum)','FontSize',8);
                    end
                    if ridx==1
                        title(sh1,sprintf('selected %s cluster',eval(cat(2,'DATA.datainfo.probe',num2str(pidx-1),'_name'))),'FontSize',8);
                        title(sh2,sprintf('d_{shortest} to nn %s sites',eval(cat(2,'DATA.datainfo.probe',num2str(pidx-1),'_name'))),'FontSize',8);
                    end
                end
                if pidx==0
                    ylabel(sh1,num2str(pc_idx),'FontSize',8);
                end
            end
        end
        % Report current estimate in the waitbar's message field
        done=ridx/n_site;
        waitbar(done,waitbar_handle,sprintf('%3.1f%%',100*done));
    end
    fh=figure('Name','Nearest Neighbour Cluster Distance to Their Centroid Summary',...
        'NumberTitle','off',...
        'MenuBar','none',...
        'ToolBar','figure',...
        'Position',[0,0,900,600],...
        'Keypressfcn',@figure_keypress);
    fname=fieldnames(DATA.datainfo);
    fidx=find(cellfun(@(x)~isempty(x),regexp(fname,'probe\w_name')));
    plotcount=1;
    for idx=1:DATA.datainfo.n_probe
        if idx<=size(pairedDist,1)
            temp=cell2mat(pairedDist(idx,:))';
            if ~isempty(temp)
                [N,edges] = histcounts(temp,max(round(numel(temp)/20),10));
                sh=subplot(2,1,1,'Parent',fh);
                line(sh,edges(1:end-1),N,'Marker','o','MarkerSize',3,'LineStyle','-','LineWidth',1,'Color',DATA.datainfo.probe_colours(idx),'MarkerFaceColor',DATA.datainfo.probe_colours(idx));
                probe_list{plotcount}=DATA.datainfo.(fname{fidx(idx)}); %#ok<AGROW>
                plotcount=plotcount+1;
            end
            if ~isempty(centroid_dist{idx})
                [N,edges] = histcounts(centroid_dist{idx},max(round(numel(centroid_dist{idx})/20),5));
                sh=subplot(2,1,2,'Parent',fh);
                line(sh,edges(1:end-1),N,'Marker','o','MarkerSize',3,'LineStyle','-','LineWidth',1,'Color',DATA.datainfo.probe_colours(idx),'MarkerFaceColor',DATA.datainfo.probe_colours(idx));
            end
        end
    end
    sh=subplot(2,1,1,'Parent',fh);
    xlabel(sh,'r (\mum)','FontSize',8);
    title(sh,sprintf('distribution of mean nn site distance between nn clusters'));
    legend(sh,probe_list);
    sh=subplot(2,1,2,'Parent',fh);
    xlabel(sh,'r (\mum)','FontSize',8);
    title(sh,sprintf('distribution of centroid distance between nn clusters'));
    legend(sh,probe_list);
    delete(waitbar_handle);       % DELETE the waitbar; don't try to CLOSE it.
catch exception
    if exist('waitbar_handle','var')&&ishandle(waitbar_handle)
        delete(waitbar_handle);
    end
    errordlg(exception.message,'Calculation Error','modal');
end

function identify_synapse(handles)
global DATA;
try
    probeidx=handles.MENU_PROBE.Value;
    %selectedrow=handles.TABLE_CLUSTERINFO.UserData;
    % update probe menu
    fname=fieldnames(DATA.datainfo);
    fidx=find(cellfun(@(x)~isempty(x),regexp(fname,'probe\w_name')));
    probe_list=cell(DATA.datainfo.n_probe,1);
    for idx=1:DATA.datainfo.n_probe
        probe_list{idx}=DATA.datainfo.(fname{fidx(idx)});
    end
    set(0,'DefaultUicontrolBackgroundColor',[0.3,0.3,0.3]);
    set(0,'DefaultUicontrolForegroundColor','k');
    [s,v] = listdlg('PromptString','Select neighbouring probe:',...
        'SelectionMode','single',...
        'InitialValue',numel(probe_list),...
        'ListString',probe_list);
    set(0,'DefaultUicontrolBackgroundColor','k');
    set(0,'DefaultUicontrolForegroundColor','w');
    if v%validated by click on OK
        n_site=numel(DATA.probe(probeidx).cluster);
        pdistrec=cell(n_site,1);
        distrec=cell(n_site,1);
        %create waitbar for user attention
        waitbar_handle = waitbar(0,'Please wait...','Progress Bar','Calculating...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)',...
            'WindowStyle','normal',...
            'Color',[0.2,0.2,0.2]);
        setappdata(waitbar_handle,'canceling',0);
        for pc_idx=1:n_site
            po_idx=DATA.probe(probeidx).cluster(pc_idx).(cat(2,'nnc',num2str(s-1),'_id'));
            %pc_cluster=DATA.probe(probeidx).cluster(pc_idx);
            %po_cluster=DATA.probe(s).cluster(po_idx);
            %pc_pts=pc_cluster.shape.Points;
            %po_pts=po_cluster.shape.Points;
            pc_pts=DATA.probe(probeidx).location(DATA.probe(probeidx).cluster(pc_idx).index,:);
            po_pts=DATA.probe(s).location(DATA.probe(s).cluster(po_idx).index,:);
            [d,Ind]=pdist2(po_pts,pc_pts,'euclidean','Smallest',1);
            synapse_po=Ind(d<=DATA.datainfo.synaptic_range(1));
            if ~isempty(synapse_po)
                if (numel(synapse_po)/size(pc_pts,1)>=DATA.datainfo.per_synapse)
                    [d,Ind]=pdist2(pc_pts,po_pts,'euclidean','Smallest',1);
                    synapse_pc=Ind(d<=DATA.datainfo.synaptic_range(1));
                    [~,Centre,~,Dist] = kmeans([pc_pts(synapse_pc,:);po_pts(synapse_po,:)],1,'Distance','sqeuclidean');
                    meanptc_dist=[mean(Dist),median(Dist),max(Dist)];
                    pdistrec{pc_idx,1}=d;
                    viewvec=mean([DATA.probe(probeidx).cluster(pc_idx).centroid-Centre;DATA.probe(s).cluster(po_idx).centroid-Centre],1);
                    DATA.probe(probeidx).cluster(pc_idx).view=viewvec;
                    %[~,pind]=min(d);
                    %loc=(pc_pts(pind(1),:)+po_pts(Ind(pind(1)),:))/2;%half way of closes site pair
                    loc=[nan,nan,nan];
                    distrec{pc_idx,1}=Dist';
                else
                    Centre=[nan,nan,nan];
                    loc=[nan,nan,nan];
                    meanptc_dist=[nan,nan,nan];
                end
            else
                Centre=[nan,nan,nan];
                loc=[nan,nan,nan];
                meanptc_dist=[nan,nan,nan];
            end
            DATA.probe(probeidx).cluster(pc_idx).('centroid_synapse')=Centre;
            DATA.probe(probeidx).cluster(pc_idx).('loc_synapse')=loc;
            DATA.probe(probeidx).cluster(pc_idx).('mean_dist_synapse')=meanptc_dist;
            % Report current estimate in the waitbar's message field
            done=pc_idx/n_site;
            waitbar(done,waitbar_handle,sprintf('%3.1f%%',100*done));
        end
        delete(waitbar_handle);       % DELETE the waitbar; don't try to CLOSE it.
        fh=figure('Name',sprintf('%s synapse search against %s',probe_list{probeidx},probe_list{s}),...
            'NumberTitle','off',...
            'MenuBar','none',...
            'ToolBar','figure',...
            'Position',[0,0,900,600],...
            'Keypressfcn',@figure_keypress);
        temp=[pdistrec{:}];
        subplot(2,1,1,'Parent',fh);
        histogram(temp(:),50);
        xlabel('r (\mum)','FontSize',8);
        title(sprintf('%s to %s site synapse pairwise distances',probe_list{probeidx},probe_list{s}));
        temp=[distrec{:}]';
        subplot(2,1,2,'Parent',fh);
        histogram(temp(:),50);
        xlabel('r (\mum)','FontSize',8);
        title(sprintf('%s and %s to synapse centroid distances',probe_list{probeidx},probe_list{s}));
        display_clusterinfo(probeidx,handles);
        msgbox(sprintf('%s synapse search against %s successfully analysed.\n',probe_list{probeidx},probe_list{s}),'Cluster Analysis','modal');
    end
catch exception
    if exist('waitbar_handle','var')&&ishandle(waitbar_handle)
        delete(waitbar_handle);
    end
    errordlg(exception.message,'Calculation Error','modal');
end

function analyse_synapse_nn_site(handles)
global DATA;
try
    % update probe menu
    fname=fieldnames(DATA.datainfo);
    fidx=find(cellfun(@(x)~isempty(x),regexp(fname,'probe\w_name')));
    probe_list=cell(DATA.datainfo.n_probe,1);
    for idx=1:DATA.datainfo.n_probe
        probe_list{idx}=DATA.datainfo.(fname{fidx(idx)});
    end
    set(0,'DefaultUicontrolBackgroundColor',[0.3,0.3,0.3]);
    set(0,'DefaultUicontrolForegroundColor','k');
    [s,v] = listdlg('PromptString','Select neighbouring probe:',...
        'SelectionMode','multiple',...
        'InitialValue',1:1:numel(probe_list),...
        'ListString',probe_list);
    set(0,'DefaultUicontrolBackgroundColor','k');
    set(0,'DefaultUicontrolForegroundColor','w');
    if v%validated by click on OK
        [pathname,~,~]=fileparts(handles.LOCSUPRES.Name);
        if isempty(pathname)
            pathname=pwd;
        end
        pathname=cat(2,pathname,filesep,'synapse_nn_site');
        if ~isdir(pathname)
            dos(cat(2,'mkdir "',pathname,'"'));
        end
        plotcount=0;
        currentprobe=handles.MENU_PROBE.Value;
        selected_cluster=handles.TABLE_CLUSTERINFO.Data(handles.TABLE_CLUSTERINFO.UserData,1);
        n_cluster=numel(selected_cluster);
        currentsynapse=cell2mat({DATA.probe(currentprobe).cluster(selected_cluster).centroid_synapse}');
        prox_dist=max(DATA.datainfo.localdist);
        set(0,'DefaultUicontrolBackgroundColor',[0.3,0.3,0.3]);
        set(0,'DefaultUicontrolForegroundColor','k');
        options.Interpreter = 'tex';
        options.WindowStyle = 'modal';
        answer = inputdlg('Enter analysis sphere radius (\mum)[default use max(localdist)]:',...
            'Synapse sphere radius', 1,{num2str(prox_dist)},options);
        set(0,'DefaultUicontrolBackgroundColor','k');
        set(0,'DefaultUicontrolForegroundColor','w');
        if ~isempty(prox_dist)
            prox_dist=max(str2double(answer),0.1);
        end
        for clusterid=1:n_cluster
            fh=figure('Name',sprintf('Nearest Probe Site Distance to cluster%g synapse',selected_cluster(clusterid)),...
                'NumberTitle','off',...
                'MenuBar','none',...
                'ToolBar','figure',...
                'Position',[0,0,900,600],...
                'Color',[0.5,0.5,0.5],...
                'Keypressfcn',@figure_keypress);
            sph=subplot(2,numel(s)*3,[(numel(s)*2+1):numel(s)*3,(numel(s)*5+1):2*(numel(s)*3)],'Parent',fh);hold all;
            sph.Tag='synapse_cluster';
            plot3(handles.PANEL_CLUSTER,currentsynapse(clusterid,1),currentsynapse(clusterid,2),currentsynapse(clusterid,3),...
                'LineStyle','none','LineWidth',4,...
                'MarkerSize',5,'Marker','+',...
                'Color','w','Parent',sph);
            synapse_centre=currentsynapse(clusterid,:);
            for probeidx=1:numel(s)%go through each probe including self if selected
                othersite=DATA.probe(s(probeidx)).location;
                d_len=bsxfun(@minus,othersite,currentsynapse(clusterid,:));
                [az,el,rad]=cart2sph(d_len(:,1),d_len(:,2),d_len(:,3));
                proximity=rad<=prox_dist;
                rad=rad(proximity);
                az=rad2deg(az(proximity));el=rad2deg(el(proximity));
                %viewaz(probeidx)=mean(az);viewel(probeidx)=mean(el);
                if isempty(DATA.probe(s(probeidx)).cluster)
                    vf_shape=0;rou_shape=0;A_surf_shape=0;
                    az_shape=[];el_shape=[];rad_shape=nan;
                else
                    if s(probeidx)==currentprobe
                        %selected probe
                        checkcluster=selected_cluster(clusterid);
                    else
                        fname=cat(2,'nnc',num2str(s(probeidx)-1),'_id');
                        if isfield(DATA.probe(currentprobe).cluster(selected_cluster(clusterid)),fname)
                            checkcluster=DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(fname);
                        end
                    end
                    if isfield(DATA.probe(s(probeidx)).cluster(checkcluster),'shape')
                        d_len_shape=bsxfun(@minus,DATA.probe(s(probeidx)).cluster(checkcluster).shape.Points,currentsynapse(clusterid,:));
                        [az_shape,el_shape,rad_shape]=cart2sph(d_len_shape(:,1),d_len_shape(:,2),d_len_shape(:,3));
                        az_shape=rad2deg(az_shape);el_shape=rad2deg(el_shape);
                        V_cluster=DATA.probe(s(probeidx)).cluster(checkcluster).volume;
                        vf_shape=V_cluster/(4/3*pi*prox_dist^3);
                        rou_shape=size(DATA.probe(s(probeidx)).cluster(checkcluster).shape.Points,1)/V_cluster;
                        A_surf_shape=DATA.probe(s(probeidx)).cluster(checkcluster).area;
                    else
                        az_shape=[];el_shape=[];rad_shape=nan;
                        vf_shape=0;rou_shape=0;A_surf_shape=0;
                    end
                end
                %calculate volume information
                temp=alphaShape(d_len(proximity,1:3),Inf,'HoleThreshold',0,'RegionThreshold',0);
                pc=criticalAlpha(temp,'one-region');
                if isempty(pc)
                    vf_max=0;rou_max=0;A_surf_max=0;
                else
                    temp.Alpha=pc;
                    V_cluster_max=volume(temp);
                    vf_max=V_cluster_max/(4/3*pi*prox_dist^3);
                    rou_max=numel(rad)/V_cluster_max;
                    A_surf_max=surfaceArea(temp);
                end
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.(cat(2,'probe',num2str(s(probeidx))))=temp;
                subplot(2,numel(s)*3,(probeidx-1)*2+1,'Parent',fh);
                histogram2(az,el,-180:10:180,-90:10:90,'FaceColor','flat','DisplayStyle','tile','EdgeColor','none');
                view([0,90]);axis('equal');xlabel('az (^o)','FontSize',8);ylabel('el (^o)','FontSize',8);
                title(sprintf('Probe\n%s\n(point)',probe_list{probeidx}),'Color',DATA.datainfo.probe_colours(s(probeidx)));
                subplot(2,numel(s)*3,(probeidx-1)*2+2,'Parent',fh);
                histogram2(az_shape,el_shape,-180:10:180,-90:10:90,'FaceColor','flat','DisplayStyle','tile','EdgeColor','none');
                view([0,90]);axis('equal');xlabel('az (^o)','FontSize',8);ylabel('el (^o)','FontSize',8);
                title(sprintf('Probe\n%s\n(shape)',probe_list{probeidx}),'Color',DATA.datainfo.probe_colours(s(probeidx)));
                subplot(2,numel(s)*3,numel(s)*3+(probeidx-1)*2+1,'Parent',fh);
                histogram(rad,linspace(0,prox_dist,25),'FaceColor',DATA.datainfo.probe_colours(s(probeidx)));
                xlabel('r (\mum)','FontSize',8);
                title({cat(2,'r_{min} = ',sprintf('%4.3f',min(rad)),'\mum');...
                    cat(2,'r_{mean} = ',sprintf('%4.3f',mean(rad)),'\mum');...
                    cat(2,'r_{median} = ',sprintf('%4.3f',median(rad)),'\mum');...
                    cat(2,'V_f = ',sprintf('%4.2f',vf_max*100),'%');...
                    cat(2,'\rho = ',sprintf('%4.0f',rou_max),'\mum^{-3}');...
                    cat(2,'A_{surf} = ',sprintf('%4.2f',A_surf_max),'\mum^2')},...
                    'Interpreter','tex');
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_pt_rmin'))=min(rad);
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_pt_rmean'))=mean(rad);
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_pt_rmedian'))=median(rad);
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_pt_Vf'))=vf_max*100;
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_pt_density'))=rou_max;
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_pt_Asurf'))=A_surf_max;
                %{
                title({'parameter = point|shape';...
                    cat(2,'r_{min} = ',sprintf('%4.3f | %4.3f',min(rad),min(rad_shape)),'\mum');...
                    cat(2,'r_{mean} = ',sprintf('%4.3f | %4.3f',mean(rad),mean(rad_shape)),'\mum');...
                    cat(2,'r_{median} = ',sprintf('%4.3f | %4.3f',median(rad),median(rad_shape)),'\mum');...
                    cat(2,'V_f = ',sprintf('%4.2f | %4.2f',vf_max*100,vf_shape*100),'%');...
                    cat(2,'\rho = ',sprintf('%4.2f | %4.2f',rou_max,rou_shape),'\mum^{-3}');...
                    cat(2,'A_{surface} = ',sprintf('%4.2f | %4.2f',A_surf_max,A_surf_shape),'\mum^2')},...
                    'Interpreter','tex');
                %}
                subplot(2,numel(s)*3,numel(s)*3+(probeidx-1)*2+2,'Parent',fh);
                histogram(rad_shape,linspace(0,prox_dist,25),'FaceColor',DATA.datainfo.probe_colours(s(probeidx)));
                xlabel('r (\mum)','FontSize',8);
                title({cat(2,'r_{min} = ',sprintf('%4.3f',min(rad_shape)),'\mum');...
                    cat(2,'r_{mean} = ',sprintf('%4.3f',mean(rad_shape)),'\mum');...
                    cat(2,'r_{median} = ',sprintf('%4.3f',median(rad_shape)),'\mum');...
                    cat(2,'V_f = ',sprintf('%4.2f',vf_shape*100),'%');...
                    cat(2,'\rho = ',sprintf('%4.0f',rou_shape),'\mum^{-3}');...
                    cat(2,'A_{surf} = ',sprintf('%4.2f',A_surf_shape),'\mum^2')},...
                    'Interpreter','tex');
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_sh_rmin'))=min(rad_shape);
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_sh_rmean'))=mean(rad_shape);
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_sh_rmedian'))=median(rad_shape);
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_sh_Vf'))=vf_shape*100;
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_sh_density'))=rou_shape;
                DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nns',num2str(s(probeidx)-1),'_sh_Asurf'))=A_surf_shape;
                % export raw angle and distance data
                filename=sprintf('%s%scluster%d_%s.dat',pathname,filesep,selected_cluster(clusterid),probe_list{probeidx});
                fid=fopen(filename,'w');
                fprintf(fid,'%4.4g,%4.4g,%4.4g\n',[az';el';rad']);
                fclose(fid);
                subplot(2,numel(s)*3,[(numel(s)*2+1):numel(s)*3,(numel(s)*5+1):2*(numel(s)*3)],'Parent',fh);
                if isfield(DATA.probe(s(probeidx)).cluster,'shape')
                    % has cluster defined
                    if s(probeidx)==currentprobe
                        % if it is selected one
                        plot(DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).shape,...
                            'Facecolor',DATA.datainfo.probe_colours(currentprobe),...
                            'EdgeColor','k',...
                            'EdgeAlpha',0.2,...
                            'FaceAlpha',0.4,...
                            'Parent',sph);
                        %viewvec(probeidx,:)=DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).centroid-synapse_centre;
                    else
                        %if ~isfield(DATA.probe(currentprobe).cluster(selected_cluster(clusterid)),cat(2,'nnc',num2str(s(probeidx)-1),'_id'))
                        pco_idx=DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nnc',num2str(s(probeidx)-1),'_id'));
                        plot(DATA.probe(s(probeidx)).cluster(pco_idx).shape,...
                            'Facecolor',DATA.datainfo.probe_colours(s(probeidx)),...
                            'EdgeColor','k',...
                            'EdgeAlpha',0.2,...
                            'FaceAlpha',0.4,...
                            'Parent',sph);
                    end
                else
                    %temporary shape
                    %{
                    temp.Points=bsxfun(@plus,temp.Points,currentsynapse(clusterid,:));
                        plot(temp,...
                            'Facecolor',DATA.datainfo.probe_colours(s(probeidx)),...
                            'EdgeColor','k',...
                            'EdgeAlpha',0.2,...
                            'FaceAlpha',0.4,...
                            'Parent',sph);
                    %}
                    % just scatter points
                    
                    probesite=DATA.probe(s(probeidx)).location(proximity,:);
                    plot3(probesite(:,1),probesite(:,2),probesite(:,3),...
                        'LineStyle','none','Marker','o','Color',DATA.datainfo.probe_colours(s(probeidx)),...
                        'MarkerSize',3,'Parent',sph);
                    
                    %viewaz(probeidx)=nan;viewel(probeidx)=nan;
                end
                daspect(sph,[1 1 1]);
                xlim(sph,'auto');ylim(sph,'auto');zlim(sph,'auto');
            end
            plotcount=plotcount+1;
            sph.Color=[0.5,0.5,0.5];
            sph.XColor=[0.7,0.7,0.7];
            sph.YColor=[0.7,0.7,0.7];
            sph.ZColor=[0.7,0.7,0.7];
            sph.GridColor=[0.5,0.5,0.5];
            sph.MinorGridColor=[0.5,0.5,0.5];
            sph.XGrid='on';xlabel(sph,'X');
            sph.YGrid='on';ylabel(sph,'Y');
            sph.ZGrid='on';zlabel(sph,'Z');
            axis(sph,'equal');
            if isfield(DATA.probe(currentprobe).cluster(selected_cluster(clusterid)),'view')
                normvec=DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).view;
                [az,el,~]=cart2sph(normvec(1),normvec(2),normvec(3));
                view(sph,[rad2deg(az)-90,rad2deg(el)-90]);
                %{
            transl=eye(4);%transl(4,1:3)=-synapse_centre;
            az=-az;el=-el;
            rotz=[cos(az),sin(az),0,0;...
                -sin(az),cos(az),0,0;...
                0,0,1,0;...
                0,0,0,1];
            roty=[cos(el),0,-sin(el),0;...
                0,1,0,0;...
                sin(el),0,cos(el),0;...
                0,0,0,1];
            transM=roty*rotz*transl;
            tform=affine3d(transM);
            [x1,y1,z1] = transformPointsForward(tform,DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe1.Points(:,1),DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe1.Points(:,2),DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe1.Points(:,3));
            [x2,y2,z2] = transformPointsForward(tform,DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe2.Points(:,1),DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe2.Points(:,2),DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe2.Points(:,3));
            [x3,y3,z3] = transformPointsForward(tform,DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe3.Points(:,1),DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe3.Points(:,2),DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe3.Points(:,3));
            [x4,y4,z4] = transformPointsForward(tform,DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe4.Points(:,1),DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe4.Points(:,2),DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.probe4.Points(:,3));
            figure(10);plot3(x1,y1,z1,'o','Color',DATA.datainfo.probe_colours(1));hold all;
            plot3(x2,y2,z2,'o','Color',DATA.datainfo.probe_colours(2));
            plot3(x3,y3,z3,'.','Color',DATA.datainfo.probe_colours(3));
            plot3(x4,y4,z4,'.','Color',DATA.datainfo.probe_colours(4));
            ax=gca;axis(ax,'equal');xlabel(ax,'X');ylabel(ax,'Y');zlabel(ax,'Z');
                %}
            end
            %view(sph,DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).view);
        end
        % update cluster info table
        display_clusterinfo(currentprobe,handles);
        msgbox(sprintf('cluster %s synapse nearest neighbour site search successfully analysed.\n',sprintf('%d, ',selected_cluster)),'Cluster Analysis','modal');
    else
        errordlg(sprintf('synapse nearest site analysis cancelled\n'));
    end
catch exception
    if exist('waitbar_handle','var')&&ishandle(waitbar_handle)
        delete(waitbar_handle);
    end
    errordlg(exception.message,'Calculation Error','modal');
end

function analyse_synapse_subcluster_site(handles)
global DATA;
try
    % update probe menu
    fname=fieldnames(DATA.datainfo);
    fidx=find(cellfun(@(x)~isempty(x),regexp(fname,'probe\w_name')));
    probe_list=cell(DATA.datainfo.n_probe,1);
    for idx=1:DATA.datainfo.n_probe
        probe_list{idx}=DATA.datainfo.(fname{fidx(idx)});
    end
    set(0,'DefaultUicontrolBackgroundColor',[0.3,0.3,0.3]);
    set(0,'DefaultUicontrolForegroundColor','k');
    [s,v] = listdlg('PromptString','Select neighbouring probe:',...
        'SelectionMode','multiple',...
        'InitialValue',1:1:numel(probe_list),...
        'ListString',probe_list);
    set(0,'DefaultUicontrolBackgroundColor','k');
    set(0,'DefaultUicontrolForegroundColor','w');
    if v%validated by click on OK
        [pathname,~,~]=fileparts(handles.LOCSUPRES.Name);
        if isempty(pathname)
            pathname=pwd;
        end
        plotcount=0;
        prox_dist=max(DATA.datainfo.localdist);
        currentprobe=handles.MENU_PROBE.Value;
        selected_cluster=handles.TABLE_CLUSTERINFO.Data(handles.TABLE_CLUSTERINFO.UserData,1);
        n_cluster=numel(selected_cluster);
        currentsynapse=cell2mat({DATA.probe(currentprobe).cluster(selected_cluster).centroid_synapse}');
        for clusterid=1:n_cluster
            fh=figure('Name',sprintf('Sub-cluster analysis for cluster%g synapse',selected_cluster(clusterid)),...
                'NumberTitle','off',...
                'MenuBar','none',...
                'ToolBar','figure',...
                'Position',[0,0,900,600],...
                'Color',[0.5,0.5,0.5],...
                'Keypressfcn',@figure_keypress);
            sph=subplot(2,numel(s)+2,[numel(s)+1,numel(s)+2,2*(numel(s)+1)+1,2*(numel(s)+1)+2],'Parent',fh);hold all;
            sph.Tag='synapse_ncluster';
            plot3(handles.PANEL_CLUSTER,currentsynapse(clusterid,1),currentsynapse(clusterid,2),currentsynapse(clusterid,3),...
                'LineStyle','none','LineWidth',4,...
                'MarkerSize',5,'Marker','+',...
                'Color','w','Parent',sph);
            % work out orientation
            for probeidx=1:numel(s)%go through each probe including self if selected
                probecluster=DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.(cat(2,'probe',num2str(s(probeidx))));
                if ~isempty(probecluster.Points)
                    clustertree=linkage(probecluster.Points,'ward','euclidean','savememory','on');
                    %c = cluster(clustertree,'cutoff',0.3,'depth',1,'criterion','inconsistent');
                    c = cluster(clustertree,'cutoff',2,'depth',3,'criterion','inconsistent');
                    %c = cluster(clustertree,'maxclust',5);
                    subplot(2,numel(s)+2,probeidx,'Parent',fh);
                    dendrogram(clustertree);
                    xlabel('index','FontSize',8);
                    ylabel('level','FontSize',8);
                    title(sprintf('Probe %s',probe_list{probeidx}));
                    subplot(2,numel(s)+2,probeidx+numel(s)+2,'Parent',fh);
                    scatter3(probecluster.Points(:,1)+currentsynapse(clusterid,1),...
                        probecluster.Points(:,2)+currentsynapse(clusterid,2),...
                        probecluster.Points(:,3)+currentsynapse(clusterid,3),...
                        10,c,'filled');
                    view(DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).view);
                    xlim([currentsynapse(clusterid,1)-prox_dist,currentsynapse(clusterid,1)+prox_dist]);
                    ylim([currentsynapse(clusterid,2)-prox_dist,currentsynapse(clusterid,2)+prox_dist]);
                    zlim([currentsynapse(clusterid,3)-prox_dist,currentsynapse(clusterid,3)+prox_dist]);
                    axis('equal');
                    % export raw angle and distance data
                    %filename=sprintf('%s%scluster%d_%s.dat',pathname,filesep,selected_cluster(clusterid),probe_list{probeidx});
                    %fid=fopen(filename,'w');
                    %fprintf(fid,'%4.4g,%4.4g,%4.4g\n',[az';el';rad']);
                    %fclose(fid);
                    subplot(2,numel(s)+2,[numel(s)+1,numel(s)+2,2*(numel(s)+1)+1,2*(numel(s)+1)+2],'Parent',fh);
                    if ~isfield(DATA.probe(currentprobe).cluster(selected_cluster(clusterid)),cat(2,'nnc',num2str(s(probeidx)-1),'_id'))
                        probesite=DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).synapse.(cat(2,'probe',num2str(s(probeidx)))).Points;
                        plot3(probesite(:,1)+currentsynapse(clusterid,1),...
                            probesite(:,2)+currentsynapse(clusterid,2),...
                            probesite(:,3)+currentsynapse(clusterid,3),...
                            'LineStyle','none','Marker','o','Color',DATA.datainfo.probe_colours(s(probeidx)),...
                            'MarkerSize',3,'Parent',sph);
                    elseif s(probeidx)==currentprobe
                        plot(DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).shape,...
                            'Facecolor',DATA.datainfo.probe_colours(currentprobe),...
                            'EdgeColor','k',...
                            'EdgeAlpha',0.2,...
                            'FaceAlpha',0.4,...
                            'Parent',sph);
                    else
                        pco_idx=DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).(cat(2,'nnc',num2str(s(probeidx)-1),'_id'));
                        plot(DATA.probe(s(probeidx)).cluster(pco_idx).shape,...
                            'Facecolor',DATA.datainfo.probe_colours(s(probeidx)),...
                            'EdgeColor','k',...
                            'EdgeAlpha',0.2,...
                            'FaceAlpha',0.4,...
                            'Parent',sph);
                    end
                    daspect(sph,[1 1 1]);
                    xlim(sph,'auto');ylim(sph,'auto');zlim(sph,'auto');
                end
            end
            plotcount=plotcount+1;
            sph.Color=[0.5,0.5,0.5];
            sph.XColor=[0.7,0.7,0.7];
            sph.YColor=[0.7,0.7,0.7];
            sph.ZColor=[0.7,0.7,0.7];
            sph.GridColor=[0.5,0.5,0.5];
            sph.MinorGridColor=[0.5,0.5,0.5];
            sph.XGrid='on';xlabel(sph,'X');
            sph.YGrid='on';ylabel(sph,'Y');
            sph.ZGrid='on';zlabel(sph,'Z');
            axis(sph,'equal');
            view(sph,DATA.probe(currentprobe).cluster(selected_cluster(clusterid)).view);
        end
        msgbox(sprintf('cluster %d synapse synapse sub-cluster successfully analysed.\n',selected_cluster),'Sub-Cluster Analysis','modal');
    else
        errordlg(sprintf('synapse sub-cluster analysis cancelled\n'));
    end
catch exception
    if exist('waitbar_handle','var')&&ishandle(waitbar_handle)
        delete(waitbar_handle);
    end
    errordlg(exception.message,'Calculation Error','modal');
end
