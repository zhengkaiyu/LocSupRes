function [ status, message ] = export_panel( panel_handle_inarg )
%EXPORT_PANEL export different type of plots from matlab
%   export plot type
%   trace,histogram,surface,scatter

%% function complete

% assume worst
status=false;message='';
% default to current path
rootpath=cat(2,pwd,filesep);
PathName=rootpath;
try
    for panel_idx=1:numel(panel_handle_inarg)
        panel_handle=panel_handle_inarg(panel_idx);
        panel_tag=cat(2,panel_handle.Tag,num2str(panel_idx));
        if ishandle(panel_handle)
            % --- axes has child plot ---
            %{
            subplots=findobj(panel_handle,'Type','axes');
            if ~isempty(subplots)
                [status,message]=export_panel(panel_handle.Children);
            end
            %}
            % --- export histograms ---
            histplot=findobj(panel_handle,'Type','histogram');% find trace plot
            if ~isempty(histplot)
                % has trace line plot
                % --- Calculate Data ---
                trace_x=get(histplot,'BinEdges');% get xdata
                trace_y=get(histplot,'Values');% get ydata
                if size(trace_x,1)==1%single line to export
                    % turn single trace into cell data for easy processing
                    trace_x={trace_x};
                end
                if size(trace_y,1)==1
                    trace_y={trace_y};
                end
                if ~isempty(trace_y)
                    % create data cell holder
                    data=cell(numel(trace_x)*2,1);
                    data(1:2:end)=flipud(trace_x);data(2:2:end)=flipud(trace_y);
                    % --- Ask for output ---
                    button=questdlg('Where do you want the hist plot data exported to?>','Export plots','Clipboard','File','Cancel','Clipboard');
                    switch button
                        case {'Cancel',''}
                            % if user cancelled action
                            message=sprintf('%s\n','exporting hist cancelled');
                        case 'Clipboard'
                            % send data to clipboard
                            data2clip(cellfun(@(x)x',data,'UniformOutput',false));
                            message=sprintf('%g histograms exported to %s\n',numel(trace_x),'Clipboard');
                            status=true;
                        case 'File'
                            % work out automatic filename
                            tag=cat(2,rootpath,'hist_',panel_tag,'_',datestr(now,'yyyymmddHHMMSS'));
                            % ask for save file confirmation
                            [FileName,PathName,~]=uiputfile('*.dat','Export Traces',tag);
                            if ischar(FileName)
                                % if user happy
                                FileName=cat(2,PathName,FileName);% get filename
                                % open output file for writing
                                fid = fopen(FileName,'w');
                                % output data into file
                                for i=1:1:numel(data)
                                    fprintf(fid,'%d ',data{i});
                                    fprintf(fid,'\n');
                                end
                                fclose(fid);% close file
                                message=sprintf('%g histograms exported to %s\n',numel(trace_x),FileName);
                                % update saved path
                                rootpath=PathName;
                                status=true;
                            else
                                % if user cancelled action
                                message=sprintf('%s\n','exporting histogram cancelled');
                            end
                    end
                else
                    % no line object found
                    message=sprintf('\n%s\n','No histograms found');
                end
            else
                % no line object found
                message=sprintf('\n%s\n','No histograms found');
            end
            % --- export histograms2d ---
            histplot=findobj(panel_handle,'Type','histogram2');% find trace plot
            if ~isempty(histplot)
                % has trace line plot
                % --- Calculate Data ---
                trace_x=get(histplot,'XBinEdges');% get xdata
                trace_y=get(histplot,'YBinEdges');% get xdata
                map=get(histplot,'Values');% get ydata
                if size(trace_x,1)==1%single line to export
                    % turn single trace into cell data for easy processing
                    trace_x={trace_x};
                end
                if size(trace_y,1)==1
                    trace_y={trace_y};
                end
                if ~isempty(map)
                    % --- Ask for output ---
                    button=questdlg('Where do you want the 2D hist plot data exported to?>','Export plots','Clipboard','File','Cancel','Clipboard');
                    switch button
                        case {'Cancel',''}
                            % if user cancelled action
                            message=sprintf('%s\n','exporting 2D hist cancelled');
                        case 'Clipboard'
                            % send data to clipboard
                            data2clip(map);
                            message=sprintf('2D histogram exported to %s\n','Clipboard');
                            status=true;
                        case 'File'
                            % work out automatic filename
                            tag=cat(2,rootpath,'hist2d_',panel_tag,'_',datestr(now,'yyyymmddHHMMSS'));
                            choice = questdlg(...
                                sprintf('Please choose the image export format.\nASCII format will output MxN matrix of Z-values.\nTIFF will export a MxNx3 coloured image file.'),...
                                'Choose Export Format','ASCII','TIFF','TIFF');
                            switch choice
                                case 'ASCII'
                                    [FileName,PathName,~] = uiputfile('*.asc','Export hist2D as',tag);
                                    if ischar(FileName)
                                        FileName=cat(2,PathName,FileName);
                                        save(FileName,'hist2D','-ascii');
                                        message=sprintf('%s\nmap exported in ascii format to %s\n',message,FileName);
                                        % update saved path
                                        rootpath=PathName;
                                        status=true;
                                    else
                                        message=sprintf('%s%s\n',message,'saveing map cancelled');
                                    end
                                case 'TIFF'
                                    [FileName,PathName,~] = uiputfile('*.tiff','Export surface mesh as',tag);
                                    if ischar(FileName)
                                        FileName=cat(2,PathName,FileName);
                                        databit=8;
                                        colours=flipud(map');
                                        colours=colours/max(colours(:))*2^databit;
                                        % construct tiff file
                                        tifobj = Tiff(FileName,'w');
                                        tagstruct.ImageLength=size(colours,1);
                                        tagstruct.ImageWidth=size(colours,2);
                                        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                                        tagstruct.BitsPerSample = databit;
                                        tagstruct.SamplesPerPixel = 1;
                                        %tagstruct.RowsPerStrip = 1;
                                        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                                        tagstruct.Software = 'MATLAB';
                                        tagstruct.Copyright = 'FIMAS';
                                        tifobj.setTag(tagstruct);
                                        tifobj.write(uint8(colours));
                                        % close tiff file construct
                                        tifobj.close();
                                        message=sprintf('%s\nmap exported in TIFF format to %s\n',message,FileName);
                                        % update saved path
                                        rootpath=PathName;
                                        status=true;
                                    else
                                        message=sprintf('%s\n%s\n',message,'saveing map cancelled');
                                    end
                                otherwise
                                    % if user cancelled action
                                    message=sprintf('%s\n','exporting trace cancelled');
                            end
                    end
                else
                    % no line object found
                    message=sprintf('\n%s\n','No 2D hist found');
                end
            else
                % no line object found
                message=sprintf('\n%s\n','No 2D hist found');
            end
            
            % --- export traces ---
            trace=findobj(panel_handle,'Type','line');% find trace plot
            if ~isempty(trace)
                % has trace line plot
                % --- Calculate Data ---
                trace_x=get(trace,'XData');% get xdata
                trace_y=get(trace,'YData');% get ydata
                if size(trace_x,1)==1%single line to export
                    % turn single trace into cell data for easy processing
                    trace_x={trace_x};
                end
                if size(trace_y,1)==1
                    trace_y={trace_y};
                end
                % create data cell holder
                data=cell(numel(trace_x)*2,1);
                data(1:2:end)=flipud(trace_x);data(2:2:end)=flipud(trace_y);
                % --- Ask for output ---
                button=questdlg('Where do you want the line plot data exported to?>','Export plots','Clipboard','File','Cancel','Clipboard');
                switch button
                    case {'Cancel',''}
                        % if user cancelled action
                        message=sprintf('%s\n','exporting trace cancelled');
                    case 'Clipboard'
                        % send data to clipboard
                        data2clip(cellfun(@(x)x',data,'UniformOutput',false));
                        message=sprintf('%g traces exported to %s\n',numel(trace_x),'Clipboard');
                        status=true;
                    case 'File'
                        % work out automatic filename
                        tag=cat(2,rootpath,'trace_',panel_tag,'_',datestr(now,'yyyymmddHHMMSS'));
                        % ask for save file confirmation
                        [FileName,PathName,~]=uiputfile('*.dat','Export Traces',tag);
                        if ischar(FileName)
                            % if user happy
                            FileName=cat(2,PathName,FileName);% get filename
                            % open output file for writing
                            fid = fopen(FileName,'w');
                            % output data into file
                            for i=1:1:numel(data)
                                fprintf(fid,'%d ',data{i});
                                fprintf(fid,'\n');
                            end
                            fclose(fid);% close file
                            message=sprintf('%g traces exported to %s\n',numel(trace_x),FileName);
                            % update saved path
                            rootpath=PathName;
                            status=true;
                        else
                            % if user cancelled action
                            message=sprintf('%s\n','exporting trace cancelled');
                        end
                end
            else
                % no line object found
                message=sprintf('\n%s\n','No traces found');
            end
            
            % --- export surfaces ---
            hsurface=findobj(panel_handle,'Type','surface');% find surface plot
            if ~isempty(hsurface)
                % has surf plot data
                map=get(hsurface,'ZData'); %get surface data
                button=questdlg('Where do you want the surface plot data exported to?>','Export plots','Clipboard','File','Cancel','Clipboard');
                switch button
                    case {'Cancel',''}
                        % if user cancelled action
                        message=sprintf('%s\nexporting trace cancelled\n',message);
                    case 'Clipboard'
                        % send data to clipboard
                        data2clip(map);
                        message=sprintf('%s\nmap data exported to Clipboard\n',message);
                        status=true;
                    case 'File'
                        tag=cat(2,rootpath,'surface_',panel_tag,'_',datestr(now,'yyyymmddHHMMSS'));
                        choice = questdlg(...
                            sprintf('Please choose the image export format.\nASCII format will output MxN matrix of Z-values.\nTIFF will export a MxNx3 coloured image file.'),...
                            'Choose Export Format','ASCII','TIFF','TIFF');
                        switch choice
                            case 'ASCII'
                                [FileName,PathName,~] = uiputfile('*.asc','Export surface mesh as',tag);
                                if ischar(FileName)
                                    FileName=cat(2,PathName,FileName);
                                    save(FileName,'map','-ascii');
                                    message=sprintf('%s\nmap exported in ascii format to %s\n',message,FileName);
                                    % update saved path
                                    rootpath=PathName;
                                    status=true;
                                else
                                    message=sprintf('%s%s\n',message,'saveing map cancelled');
                                end
                            case 'TIFF'
                                [FileName,PathName,~] = uiputfile('*.tiff','Export surface mesh as',tag);
                                if ischar(FileName)
                                    FileName=cat(2,PathName,FileName);
                                    databit=8;
                                    colours=get(hsurface,'CData');
                                    colours=colours/max(colours(:))*2^databit;
                                    % construct tiff file
                                    tifobj = Tiff(FileName,'w');
                                    tagstruct.ImageLength=size(colours,1);
                                    tagstruct.ImageWidth=size(colours,2);
                                    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                                    tagstruct.BitsPerSample = databit;
                                    tagstruct.SamplesPerPixel = 1;
                                    %tagstruct.RowsPerStrip = 1;
                                    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                                    tagstruct.Software = 'MATLAB';
                                    tagstruct.Copyright = 'FIMAS';
                                    tifobj.setTag(tagstruct);
                                    tifobj.write(uint8(colours));
                                    % close tiff file construct
                                    tifobj.close();
                                    message=sprintf('%s\nmap exported in TIFF format to %s\n',message,FileName);
                                    % update saved path
                                    rootpath=PathName;
                                    status=true;
                                else
                                    message=sprintf('%s\n%s\n',message,'saveing map cancelled');
                                end
                            otherwise
                                % if user cancelled action
                                message=sprintf('%s\n','exporting trace cancelled');
                        end
                end
            else
                % didn't find surf plot data
                message=sprintf('\n%s%s\n',message,'No surface data found');
            end
            % --- export rgbimage ---
            himage=findobj(panel_handle,'Type','image');% find surface plot
            if ~isempty(himage)
                % has surf plot data
                map=get(himage,'CData'); %get surface data
                button=questdlg('Where do you want the rgbimage plot data exported to?>','Export plots','Clipboard','File','Cancel','Clipboard');
                switch button
                    case {'Cancel',''}
                        % if user cancelled action
                        message=sprintf('%s\nexporting trace cancelled\n',message);
                    case 'Clipboard'
                        % send data to clipboard
                        data2clip(map);
                        message=sprintf('%s\nmap data exported to Clipboard\n',message);
                        status=true;
                    case 'File'
                        tag=cat(2,rootpath,'image_',panel_tag,'_',datestr(now,'yyyymmddHHMMSS'));
                        choice = questdlg(...
                            sprintf('Please choose the image export format.\nASCII format will output MxNx3 matrix of C-values.\nTIFF will export a MxNx3 coloured image file.'),...
                            'Choose Export Format','MAT','TIFF','TIFF');
                        switch choice
                            case 'MAT'
                                [FileName,PathName,~] = uiputfile('*.mat','Export image as',tag);
                                if ischar(FileName)
                                    FileName=cat(2,PathName,FileName);
                                    save(FileName,'map','-mat');
                                    message=sprintf('%s\nmap exported in matlab format to %s\n',message,FileName);
                                    % update saved path
                                    rootpath=PathName;
                                    status=true;
                                else
                                    message=sprintf('%s%s\n',message,'saveing map cancelled');
                                end
                            case 'TIFF'
                                [FileName,PathName,~] = uiputfile('*.tiff','Export image as',tag);
                                if ischar(FileName)
                                    FileName=cat(2,PathName,FileName);
                                    dr=himage.UserData;
                                    databit=8;
                                    colours=map;
                                    colours=colours/max(colours(:))*2^databit;
                                    % construct tiff file
                                    tifobj = Tiff(FileName,'w');
                                    tagstruct.ImageLength=size(colours,1);
                                    tagstruct.ImageWidth=size(colours,2);
                                    tagstruct.Photometric = Tiff.Photometric.RGB;
                                    tagstruct.BitsPerSample = databit;
                                    tagstruct.SamplesPerPixel = size(colours,3);%rgb
                                    %tagstruct.ResolutionUnit=Tiff.ResolutionUnit.Centimeter;
                                    tagstruct.ResolutionUnit=1;
                                    tagstruct.XResolution=(1/dr);
                                    tagstruct.YResolution=(1/dr);
                                    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                                    tagstruct.Software = 'MATLAB';
                                    tagstruct.Copyright = 'FIMAS';
                                    tifobj.setTag(tagstruct);
                                    tifobj.write(uint8(colours));
                                    % close tiff file construct
                                    tifobj.close();
                                    message=sprintf('%s\nmap exported in TIFF format to %s\n',message,FileName);
                                    % update saved path
                                    rootpath=PathName;
                                    status=true;
                                else
                                    message=sprintf('%s\n%s\n',message,'saveing map cancelled');
                                end
                            otherwise
                                % if user cancelled action
                                message=sprintf('%s\n','exporting trace cancelled');
                        end
                end
            else
                % didn't find surf plot data
                message=sprintf('\n%s%s\n',message,'No surface data found');
            end
            % --- export scatter plot ---
            scatter=findobj(panel_handle,'Type','hggroup','-or','Type','scatter');% find surface plot
            if ~isempty(scatter)
                trace_x=(get(scatter,'XData'));
                trace_y=(get(scatter,'YData'));
                % create data cell holder
                data=[trace_x,trace_y];
                data=cellfun(@(x)x',data,'UniformOutput',false)';
                button=questdlg('Where do you want the line plot data exported to?>','Export plots','Clipboard','ASCII','MATLAB','Clipboard');
                switch button
                    case {''}
                        % if user cancelled action
                        message=sprintf('%s\nexporting trace cancelled\n',message);
                        
                    case 'Clipboard'
                        % send data to clipboard
                        data2clip(data);
                        message=sprintf('%s\nscatter points exported to Clipboard\n',message);
                        status=true;
                    case 'ASCII'
                        tag=cat(2,rootpath,'scatter_',panel_tag,'_',datestr(now,'yyyymmddHHMMSS'));
                        [FileName,PathName,~]=uiputfile('*.dat','Export scatter',tag);
                        if ischar(FileName)
                            % get filename
                            FileName=cat(2,PathName,FileName);
                            fid = fopen(FileName,'w');
                            % export data to ascii
                            for i=1:1:numel(data)
                                fprintf(fid,'%d ',data{i});
                                fprintf(fid,'\n');
                            end
                            fclose(fid);%close file
                            message=sprintf('%s\nscatter plot exported in ascii format to %s\n',message,FileName);
                            % update saved path
                            rootpath=PathName;
                            status=true;
                        else
                            % if user cancelled action
                            message=sprintf('\n%s%s\n',message,'exporting trace cancelled');
                        end
                    case {'MATLAB'}
                        tag=cat(2,rootpath,'scatter_',panel_tag,'_',datestr(now,'yyyymmddHHMMSS'));
                        [FileName,PathName,~]=uiputfile('*.mat','Export scatter',tag);
                        if ischar(FileName)
                            % get filename
                            FileName=cat(2,PathName,FileName);
                            save(FileName,'data','-mat');
                            message=sprintf('%s\nscatter plot exported in mat format to %s\n',message,FileName);
                            % update saved path
                            rootpath=PathName;
                            status=true;
                        else
                            % if user cancelled action
                            message=sprintf('\n%s%s\n',message,'exporting trace cancelled');
                        end
                end
            else
                % didn't find scatter plot data
                message=sprintf('%s%s\n',message,'No scatter data found');
            end
        else
            % couldn't find panel for some reason
            message=sprintf('Panel currently does not exist\n');
        end
    end
catch exception
    % error handling
    message=exception.message;
    errordlg(sprintf('Error Exporting %s',message),'Error Exporting Figure');
end