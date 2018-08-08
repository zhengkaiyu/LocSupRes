function [ output_args ] = combine_projection(varargin)
%COMBINE_PROJECTION Summary of this function goes here
%   Detailed explanation goes here

%------------------------------------------------------------------------
% reading options
if nargin>0
    for arginidx=1:2:(numel(nargin))
        fname=varargin{arginidx};
        fval=varargin{arginidx+1};
        switch fname
            case 'pathname'
                if ischar(fval)
                    pathname=fval;
                else
                    errordlg(sprintf('%s need to be a string\n',fname));
                end
            case 'probe_colours'
                if ischar(fval)
                    probe_colours=fval;
                else
                    errordlg(sprintf('%s need r/g/b in the first three colours\n',fname));
                end
            case 'synapsepair'
                if isnumeric(fval)
                    eval(cat(2,fname,'=',num2str(fval)));
                else
                    errordlg(sprintf('%s need to be a numeric value\n',fname));
                end
            case {'minres','maxres','dres'}
                if isnumeric(fval)
                    eval(cat(2,fname,'=',num2str(fval)));
                else
                    errordlg(sprintf('%s need to be a numeric value\n',fname));
                end
            case 'FWHM_ratio'
                if isnumeric(fval)
                    fval=max(min(fval,1),0);% FWHM percentage of max bounded by (0,1) default to 50%
                    eval(cat(2,fname,'=',num2str(fval)));
                else
                    errordlg(sprintf('%s need to be a numeric value\n',fname));
                end
            case 'nsig'
                if isnumeric(fval)
                    fval=max(min(fval,5),0.5);% signifant figure in from 1 to 5
                    eval(cat(2,fname,'=',num2str(fval)));
                else
                    errordlg(sprintf('%s need to be a numeric value between 0.5-5\n',fname));
                end
            case 'radius'
                if isnumeric(fval)
                    fval=max(min(fval,5),0.1);% spherical cutoff in from 0.1 to 5um
                    eval(cat(2,fname,'=',num2str(fval)));
                else
                    errordlg(sprintf('%s need to be a numeric value >0\n',fname));
                end
            otherwise
                errordlg(sprintf('unknown argument option ''%s''\n',fname));
        end
    end
end
%------------------------------------------------------------------------
% assign default values if not specified
if ~exist('pathname','var')
    % default folder
    pathname='D:\Temporary\Janosh\analysed';
end
% default colour coding
if ~exist('probe_colours','var')
    probe_colours='brgymc';
    % rgb index in probe colour for rgb image formation
    rgbidx=[find(probe_colours=='r'),find(probe_colours=='g'),find(probe_colours=='b')];
end
if ~exist('synapsepair','var')
    % synapse probe pair index, default to probe 2 and 3
    synapsepair=[2,3];
end
if ~exist('minres','var')
    % minimum resolution
    minres=0.005;
end
if ~exist('maxres','var')
    % maximum resolution
    maxres=0.5;
end
if ~exist('dres','var')
    % delta resolution
    dres=0.005;
end
if ~exist('FWHM_ratio','var')
    % FWHM percentage of max bounded by (0,1) default to 50%
    FWHM_ratio=max(min(0.5,1),0);
end
if ~exist('nsig','var')
    % significant figure bound by (0.5,5) default to 1
    nsig=max(min(1,5),0.5);
end
if ~exist('radius','var')
    radius=max(min(1,5),0.1);% spherical cutoff in from 0.1 to 5um
end
%------------------------------------------------------------------------
% ask for folder
[filename,pathname,~]=uigetfile({'*.mat','MAT-files (*.mat)';'*.*','All Files (*.*)'},...
    'Select exported mat file folder','MultiSelect','on',pathname);
% check folder
if ischar(pathname)
    % if single file selection convert from str2cell
    if ~iscell(filename)
        filename={filename};
    end
    % initialise rz_tot, the rz matrix, containing three probes
    rz_tot=cell(1,3);
    % fill in rz_tot with data files containing rz_tot data
    for fileidx=1:numel(filename)
        % load mat file
        temp=load(cat(2,pathname,filename{fileidx}),'-mat');
        % check file
        if isfield(temp,'data')
            % data field exist
            if sum(size(temp.data)==[2,3])==2
                % correct size 2x3
                rz_tot{1}=[rz_tot{1};cell2mat(temp.data(:,1)')];
                rz_tot{2}=[rz_tot{2};cell2mat(temp.data(:,2)')];
                rz_tot{3}=[rz_tot{3};cell2mat(temp.data(:,3)')];
            else
                errordlg(sprintf('File %s data size is not 2x3 (rz X 3probe)\n',filename{fileidx}));
            end
        else
            errordlg(sprintf('File %s do not contain the correct data field\n',filename{fileidx}));
        end
    end
    %--------------------
    if isempty(rz_tot{1})
        %no data loaded
        errordlg(sprintf('no data loaded from files\n'));
    else
        % initialise dr,dz,zminmax
        dr=minres;dz=dr;zminmax=[];
        % combined data plot
        fh_tot=figure('Name',sprintf('Nearest Probe Site Distance to synapse collective (F9=export|pageup/pagedown=adjust dr)'),...
            'NumberTitle','off',...
            'MenuBar','none',...
            'ToolBar','figure',...
            'Position',[0,0,900,600],...
            'Color',[0.5,0.5,0.5],...
            'Tag','Synapse_tot',...
            'Keypressfcn',@figure_keypress);
        % set userdata for interactive plotting
        fh_tot.UserData.dr=dr;
        fh_tot.UserData.zminmax=zminmax;
        
        %--------------------
        % plot cylindrical r_distance histogram
        subplot(2,3,1);cla;hold on;
        for probeidx=1:numel(rz_tot)
            [n,e]=histcounts(rz_tot{probeidx}(:,1),min(rz_tot{probeidx}(:,1)):dr:max(rz_tot{probeidx}(:,1)));
            if ~isempty(find(probeidx==synapsepair))
                deltaA=dr*radius;
            else
                deltaA=(dr*((radius^2-e(1:end-1).^2).^0.5+(radius^2-e(2:end).^2).^0.5));
            end
            n=n./deltaA;
            e=e(1:end-1)+dr/2;
            plot(e,n,'LineWidth',2,'Color',probe_colours(probeidx));
        end
        grid minor;xlabel('r_{cylindrical}');ylabel('\rho_{localisation}(\mum^{-2})');
        
        %--------------------
        % plot cylindrical z_distance histogram
        subplot(2,3,2);cla;hold on;
        for probeidx=1:numel(rz_tot)
            [n,e]=histcounts(rz_tot{probeidx}(:,2),min(rz_tot{probeidx}(:,2)):dz:max(rz_tot{probeidx}(:,2)));
            if ~isempty(find(probeidx==synapsepair))
                % only record synapse pair probe values for estimating
                % synaptic edges
                sigma=nsig*std(rz_tot{probeidx}(:,2));
                zminmax=[zminmax,[-sigma;sigma]];
                deltaA=dz*radius;
            else
                % correction for distribution because of spherical cutoff
                % in cylindrical coordinate for glt1
                deltaA=(dz*((radius^2-e(1:end-1).^2).^0.5+(radius^2-e(2:end).^2).^0.5));
            end
            n=n./deltaA;
            e=e(1:end-1)+dz/2;
            plot(e,n,'LineWidth',2,'Color',probe_colours(probeidx));
        end
        grid minor;xlabel('z_{cylindrical}');ylabel('\rho_{localisation}(\mum^{-2})');
        % work out estimated synapse edges
        zminmax=[min(zminmax(:)),max(zminmax(:))];
        fh_tot.UserData.zminmax=zminmax;
        
        %--------------------
        % plot 2D scatter of all localisations with estimated synapse edges
        subplot(2,3,3);cla;
        % plot synapse edges
        temp=scatter([0;0],zminmax,300,'k','.');hold on;
        temp.Tag='synapse_edge';
        % plot all probes
        for probeidx=1:numel(rz_tot)
            scatter(rz_tot{probeidx}(:,1),rz_tot{probeidx}(:,2),5,probe_colours(probeidx),'.');
        end
        axis 'square';grid minor;xlabel('r');ylabel('z');
        title(sprintf('F1/F2/F3/F4 to toggle scatter'));
        set(gca,'Tag','Synapse_tot_scatter');
        
        %--------------------
        % plot histogram of radial rz distance i.e. distance to the origin
        subplot(2,3,4);cla;hold on;
        for probeidx=1:numel(rz_tot)
            % sign preserving distances
            rzdist=sqrt(sum(rz_tot{probeidx}.^2,2)).*sign(rz_tot{probeidx}(:,1));
            [n,e]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
            deltaA=pi*dr*(2*abs(e(1:end-1))+dr);
            e=e(1:end-1)+dr/2;
            n=n./deltaA;
            plot(e,n,'LineWidth',2,'Color',probe_colours(probeidx));
        end
        grid minor;xlabel('Dist_{rz}');ylabel('\rho_{localisation}(\mum^{-2})');
        
        %--------------------
        % plot distance to the nearest synapse edge
        subplot(2,3,5);hold on;
        for probeidx=1:numel(rz_tot)
            syndist=[];rzdist=[];
            syndist(:,1)=sqrt(sum(bsxfun(@minus,rz_tot{probeidx},[0,zminmax(1)]).^2,2));
            syndist(:,2)=sqrt(sum(bsxfun(@minus,rz_tot{probeidx},[0,zminmax(2)]).^2,2));
            syndist(:,3)=sqrt(sum(rz_tot{probeidx}.^2,2));
            rzdist=min(syndist,[],2).*sign(rz_tot{probeidx}(:,1));
            [n,e]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
            e=e(1:end-1)+dr/2;
            plot(e,n,'LineWidth',2,'Color',probe_colours(probeidx));
        end
        title(sprintf('w_{synapse} from cluster %g s.f. = %gnm',nsig,1e3*abs(diff(zminmax))));
        grid minor;xlabel('Dist_{synapse}');ylabel('N_{localisation}');
        
        %--------------------
        % plot rgb image of the synapse
        subplot(2,3,6);set(gca,'Tag','synrgbimg');
        scalemin=min([min(rz_tot{1},[],1);min(rz_tot{2},[],1);min(rz_tot{3},[],1)],[],1);
        scalemax=max([max(rz_tot{1},[],1);max(rz_tot{2},[],1);max(rz_tot{3},[],1)],[],1);
        [rz_ch{1},~]=hist3(rz_tot{1},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
        [rz_ch{2},~]=hist3(rz_tot{2},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
        [rz_ch{3},~]=hist3(rz_tot{3},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
        rgbimg=cat(3,rz_ch{rgbidx(1)}./max(rz_ch{rgbidx(1)}(:)),rz_ch{rgbidx(2)}./max(rz_ch{rgbidx(2)}(:)),rz_ch{rgbidx(3)}./max(rz_ch{rgbidx(3)}(:)));
        himage=image(imrotate(rgbimg,90));axis square;
        set(himage,'UserData',dr);
        % print out dr value
        title(sprintf('dr=%g',dr));
    end
else
    %action cancelled
    errordlg('action cancelled');
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
                    case {'Synapse_tot_scatter'}
                        switch hplot.Children(str2double(eventkey.Key(2))).Visible
                            case 'off'
                                hplot.Children(str2double(eventkey.Key(2))).Visible='on';
                            case 'on'
                                hplot.Children(str2double(eventkey.Key(2))).Visible='off';
                        end
                end
            case {'home','end','pagedown','pageup'}
                hplot=gcf;
                switch hplot.Tag
                    case 'Synapse_tot'
                        dr=hplot.UserData.dr;
                        zminmax=hplot.UserData.zminmax;
                        switch eventkey.Key
                            case 'home'
                                dr=minres;
                            case 'end'
                                dr=maxres;
                            case 'pageup'
                                dr=max(dr-dres,minres);
                            case 'pagedown'
                                dr=min(dr+dres,maxres);
                        end
                        hplot.UserData.dr=dr;
                        
                        %----------------------------------
                        % get scatter data in reverse to plot order
                        subplot(2,3,3);temp=gca;
                        rz_tot{3}=[temp.Children(1).XData(:),temp.Children(1).YData(:)];
                        rz_tot{2}=[temp.Children(2).XData(:),temp.Children(2).YData(:)];
                        rz_tot{1}=[temp.Children(3).XData(:),temp.Children(3).YData(:)];
                        % recalcute scaling
                        scalemin=min([min(rz_tot{1},[],1);min(rz_tot{2},[],1);min(rz_tot{3},[],1)],[],1);
                        scalemax=max([max(rz_tot{1},[],1);max(rz_tot{2},[],1);max(rz_tot{3},[],1)],[],1);
                        dz=dr;
                        
                        %--------------------
                        % plot cylindrical r_distance histogram
                        subplot(2,3,1);lineplot=get(gca);
                        for pidx=1:numel(rz_tot)
                            [n,e]=histcounts(rz_tot{pidx}(:,1),min(rz_tot{pidx}(:,1)):dr:max(rz_tot{pidx}(:,1)));
                            if ~isempty(find(pidx==synapsepair))
                                deltaA=dr*radius;
                            else
                                deltaA=(dr*((radius^2-e(1:end-1).^2).^0.5+(radius^2-e(2:end).^2).^0.5));
                            end
                            n=n./deltaA;
                            e=e(1:end-1)+dr/2;
                            set(lineplot.Children(numel(rz_tot)+1-pidx),'XData',e,'YData',n);
                        end
                        
                        %--------------------
                        % plot cylindrical z_distance histogram
                        subplot(2,3,2);lineplot=get(gca);
                        for pidx=1:numel(rz_tot)
                            [n,e]=histcounts(rz_tot{pidx}(:,2),min(rz_tot{pidx}(:,2)):dz:max(rz_tot{pidx}(:,2)));
                            if ~isempty(find(pidx==synapsepair))
                                deltaA=dz*radius;
                            else
                                deltaA=(dz*((radius^2-e(1:end-1).^2).^0.5+(radius^2-e(2:end).^2).^0.5));
                            end
                            n=n./deltaA;
                            e=e(1:end-1)+dz/2;
                            set(lineplot.Children(numel(rz_tot)+1-pidx),'XData',e,'YData',n);
                        end
                        
                        %--------------------
                        % alter synapse edge
                        %subplot(2,3,3);temp=gca;
                        %set(temp.Children(4),'YData',zminmax);
                        %--------------------
                        % plot histogram of radial rz distance i.e. distance to the origin
                        subplot(2,3,4); lineplot=get(gca);
                        for pidx=1:numel(rz_tot)
                            % sign preserving distances
                            rzdist=sqrt(sum(rz_tot{pidx}.^2,2)).*sign(rz_tot{pidx}(:,1));
                            [n,e]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
                            deltaA=pi*dr*(2*abs(e(1:end-1))+dr);
                            e=e(1:end-1)+dr/2;
                            n=n./deltaA;
                            set(lineplot.Children(numel(rz_tot)+1-pidx),'XData',e,'YData',n);
                        end
                        
                        %--------------------
                        % plot distance to the nearest synapse edge
                        subplot(2,3,5);lineplot=get(gca);
                        for pidx=1:numel(rz_tot)
                            syndist=[];rzdist=[];
                            %{
                            % above
                            pt_sec1=rz_tot{pidx}(rz_tot{pidx}(:,2)>=max(zminmax),:);
                            dist_sec1=sqrt(sum(bsxfun(@minus,pt_sec1,[0,max(zminmax)]).^2,2)).*sign(pt_sec1(:,1));
                            % central
                            dist_sec2=rz_tot{pidx}((rz_tot{pidx}(:,2)<max(zminmax)&rz_tot{pidx}(:,2)>min(zminmax)),1);
                            % below
                            pt_sec3=rz_tot{pidx}(rz_tot{pidx}(:,2)<=min(zminmax),:);
                            dist_sec3=sqrt(sum(bsxfun(@minus,pt_sec3,[0,min(zminmax)]).^2,2)).*sign(pt_sec3(:,1));
                            % combined
                            rzdist=[dist_sec1;dist_sec2;dist_sec3];
                            %}
                            
                            syndist(:,1)=sqrt(sum(bsxfun(@minus,rz_tot{pidx},[0,zminmax(1)]).^2,2));
                            syndist(:,2)=sqrt(sum(bsxfun(@minus,rz_tot{pidx},[0,zminmax(2)]).^2,2));
                            syndist(:,3)=sqrt(sum(rz_tot{pidx}.^2,2));
                            rzdist=min(syndist,[],2).*sign(rz_tot{pidx}(:,1));
                            
                            [n,e]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
                            %deltaA=pi*dr*(2*abs(e(1:end-1))+dr);
                            %n=n./deltaA;
                            e=e(1:end-1)+dr/2;
                            set(lineplot.Children(numel(rz_tot)+1-pidx),'XData',e,'YData',n);
                        end
                        
                        %----------------------------------
                        subplot(2,3,6);
                        % calculate rgb images
                        [rz_ch{1},~]=hist3(rz_tot{1},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
                        [rz_ch{2},~]=hist3(rz_tot{2},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
                        [rz_ch{3},~]=hist3(rz_tot{3},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
                        rgbimg=cat(3,rz_ch{rgbidx(1)}./max(rz_ch{rgbidx(1)}(:)),rz_ch{rgbidx(2)}./max(rz_ch{rgbidx(2)}(:)),rz_ch{rgbidx(3)}./max(rz_ch{rgbidx(3)}(:)));
                        himage=get(gca,'Children');
                        set(himage,'CData',imrotate(rgbimg,90));axis tight;
                        % print out dr value
                        title(sprintf('dr=%g',dr));
                end
        end
    end
end