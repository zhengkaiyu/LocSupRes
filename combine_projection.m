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
            case 'probecolour'
                if ischar(fval)
                    probecolour=fval;
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
if ~exist('probecolour','var')
    probecolour='grbymc';
    % rgb index in probe colour for rgb image formation
    rgbidx=[find(probecolour=='r'),find(probecolour=='g'),find(probecolour=='b')];
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
    appendidx=0;
    aligned_cluster_pts=[];
    % fill in rz_tot with data files containing rz_tot data
    for fileidx=1:numel(filename)
        % load mat file
        temp=load(cat(2,pathname,filename{fileidx}),'-mat');
        % check file
        if isfield(temp,'cluster_pts')
            % data field exist
            if isfield(temp,'TM')
                cluster_pts=temp.cluster_pts;
                n_cluster=numel(temp.TM);
                % transform points for each cluster
                for clusteridx=1:n_cluster
                    pmidpt_synapse=temp.TM{clusteridx}(1:3)';
                    pnorm_synapse=temp.TM{clusteridx}(4:6)';
                    pnorm_p2=temp.TM{clusteridx}(7:9)';
                    pnorm_p3=temp.TM{clusteridx}(10:12)';
                    prox_dist=temp.prox_dist;
                    
                    %origin translation to synapse point
                    shift_trans=[[eye(3);-pmidpt_synapse],[0;0;0;1]];
                    %plane rotation to match xz plane to
                    %synapse line parallel to separation plane
                    rot_trans=zeros(4,4);rot_trans(4,4)=1;
                    rot_trans(1:3,1:3)=([pnorm_synapse',pnorm_p3',pnorm_p2']);%rotation to x,y,z
                    line_synpara_TM=shift_trans*rot_trans;
                    %plane rotation to match xz plane to
                    %synapse line parallel to p2 plane
                    rot_trans(1:3,1:3)=([pnorm_p2',pnorm_p3',pnorm_synapse']);
                    line_synortho_TM=shift_trans*rot_trans;
                    
                    %transform points to new local coordinate
                    %in line_synpara
                    aligned_cluster_pts{clusteridx+appendidx}=cellfun(@(x)[x,ones(size(x,1),1)]*line_synpara_TM,cluster_pts(clusteridx,1:3),'UniformOutput',false);
                    %transform points to new local coordinate
                    %in line_synortho is just swap x and z
                    
                    %transform from cart to cylindrical
                    %use angle to separate two groups in term
                    %of distance.
                    [theta,r,z]=cellfun(@(y)cart2pol(y(:,1),y(:,2),y(:,3)),aligned_cluster_pts{clusteridx+appendidx},'UniformOutput',false);
                    %first t,r,z is para
                    %transform
                    % combine rz and crop off over proxdist r/z
                    in_rzvol_pts=cellfun(@(m,n)abs(m)<=prox_dist&abs(n)<=prox_dist,r,z,'UniformOutput',false);
                    % collapse theta dimension
                    % maintain left/right polarity using theta value
                    trans_cluster_pts{clusteridx+appendidx,1}=cellfun(@(id,a,b,c)[((abs(a(id))<=pi/2)*2-1).*b(id),c(id)],in_rzvol_pts,theta,r,z,'UniformOutput',false);%0-pi/2=+ve
                    
                    % second t,r,z is ortho do the same for orthogonal transform
                    [theta,r,z]=cellfun(@(y)cart2pol(y(:,3),y(:,2),y(:,1)),aligned_cluster_pts{clusteridx+appendidx},'UniformOutput',false);
                    in_rzvol_pts=cellfun(@(m,n)abs(m)<=prox_dist&abs(n)<=prox_dist,r,z,'UniformOutput',false);
                    trans_cluster_pts{clusteridx+appendidx,2}=cellfun(@(id,a,b,c)[((abs(a(id))<=pi/2)*2-1).*b(id),c(id)],in_rzvol_pts,theta,r,z,'UniformOutput',false);%0-pi/2=+ve
                end
                appendidx=n_cluster;
            else
                errordlg(sprintf('File %s do not contain the correct TM field\n',filename{fileidx}));
            end
        else
            errordlg(sprintf('File %s do not contain the correct cluster_pts field\n',filename{fileidx}));
        end
    end
    clear temp theta r z in_rzvol_pts;
    %--------------------
    if ~exist('aligned_cluster_pts','var')&&isempty(aligned_cluster_pts{1})
        %no data loaded
        errordlg(sprintf('no data loaded from files\n'));
    else
        % initialise dr,dz,zminmax
        dr=minres;dz=dr;zminmax=[];
        n_probe=size(cluster_pts,2);
        n_cluster=numel(aligned_cluster_pts);
        %------------------------------------
        % combined data plot
        fh_tot=figure('Name',sprintf('Nearest Probe Site Distance to synapse collective (F9=export|pageup/pagedown=adjust dr)'),...
            'NumberTitle','off',...
            'MenuBar','none',...
            'ToolBar','figure',...
            'Position',[0,0,900,600],...
            'Color',[0.5,0.5,0.5],...
            'Tag','Synapse_tot',...
            'Keypressfcn',@figure_keypress);
        fh_tot.UserData.dr=dr;
        fh_tot.UserData.rgbidx=rgbidx;
        fh_tot.UserData.nsig=nsig;
        fh_tot.UserData.probecolour=probecolour;
        fh_tot.UserData.prox_dist=prox_dist;
        %trans_cluster_pts %n_cluster X 2transform (contain 1x3 probe, contain nx2 rzcoord)
        %cluster_pts %n_cluster X n_probe (contain nx3 xyzcoord)
        %-------------------------
        aligedcollect=reshape([aligned_cluster_pts{:}]',n_probe,size(trans_cluster_pts,1))';
        %parallel transform
        transidx=1;%1st transform
        rzcollect=reshape([trans_cluster_pts{:,transidx}]',n_probe,size(trans_cluster_pts,1))';
        for probeidx=1:n_probe
            % all cluster transformed loc for the probe
            rzcoord{probeidx}=cell2mat(rzcollect(:,probeidx));
            % all original cluster xyz loc for the probe aligned
            xyzcoord{probeidx}=cell2mat(aligedcollect(:,probeidx));
        end
        % work out synapse edge using parallel transform
        sigma=[nsig*std(rzcoord{1}(:,2));1.0*std(rzcoord{2}(:,2))];
        zminmax=[-sigma;sigma];
        zminmax=[min(zminmax(:)),max(zminmax(:))];
        scalemin=min(cell2mat(rzcoord'));
        scalemax=max(cell2mat(rzcoord'));
        % scatter all points
        subplot(2,3,3);cla;hold on;
        plot(zminmax,[0,0],'k-','LineWidth',5);
        subplot(2,3,1);cla;hold on;
        scatter([0;0],zminmax,30,'k','s','MarkerFaceColor','k');
        for probeidx=1:n_probe
            subplot(2,3,1);
            scatter(rzcoord{probeidx}(:,1),rzcoord{probeidx}(:,2),5,probecolour(probeidx),'.');
            %preserving r dim so that pre and post side are
            %maintained
            dist=sqrt(sum(rzcoord{probeidx}.^2,2)).*sign(rzcoord{probeidx}(:,1));
            [n,e]=histcounts(dist,-prox_dist:dr:prox_dist);
            c=e(1:end-1)+dr/2;
            deltaA=pi*(2*prox_dist)*dr*(2*abs(c)+dr);
            n=n./deltaA;
            subplot(2,3,3);
            plot(c,n,'LineWidth',2,'Color',probecolour(probeidx));
            
            [rz_ch{probeidx},~]=hist3(rzcoord{probeidx},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
        end
        subplot(2,3,1);
        axis 'square';grid minor;xlabel('r');ylabel('z');
        title(sprintf('w_{synapse} from cluster %g s.f. = %gnm',nsig,1e3*abs(diff(zminmax))));
        set(gca,'Tag','synapse_rz_scatter');
        subplot(2,3,3);
        grid minor;xlabel('dist_{rz}');ylabel('\rho_{localisation}(\mum^{-3})');
        % parallel transform combined rgbimage
        subplot(2,3,2);set(gca,'Tag','synrgbimg');
        rgbimg=cat(3,rz_ch{rgbidx(1)}./max(rz_ch{rgbidx(1)}(:)),rz_ch{rgbidx(2)}./max(rz_ch{rgbidx(2)}(:)),rz_ch{rgbidx(3)}./max(rz_ch{rgbidx(3)}(:)));
        himage=image(imrotate(rgbimg,90));axis square;title(sprintf('dr=%g',dr));
        set(himage,'UserData',dr);
        
        %orthogonal transform
        transidx=2;%1st transform
        rzcollect=reshape([trans_cluster_pts{:,transidx}]',n_probe,n_cluster)';
        for probeidx=1:n_probe
            % all cluster transformed loc for the probe
            rzcoord{probeidx}=cell2mat(rzcollect(:,probeidx));
        end
        % work out synapse edge using parallel transform
        sigma=[nsig*std(rzcoord{1}(:,1));1.0*std(rzcoord{2}(:,1))];
        rminmax=[-sigma;sigma];
        rminmax=[min(rminmax(:)),max(rminmax(:))];
        scalemin=min(cell2mat(rzcoord'));
        scalemax=max(cell2mat(rzcoord'));
        % scatter all points
        subplot(2,3,6);cla;hold on;
        plot(rminmax,[0,0],'k-','LineWidth',5);
        subplot(2,3,4);cla;hold on;
        scatter(rminmax,[0;0],30,'k','s','MarkerFaceColor','k');
        for probeidx=1:n_probe
            subplot(2,3,4);
            scatter(rzcoord{probeidx}(:,1),rzcoord{probeidx}(:,2),5,probecolour(probeidx),'.');
            %preserving z dim so that pre and post side are
            %maintained
            dist=sqrt(sum(rzcoord{probeidx}.^2,2)).*sign(rzcoord{probeidx}(:,1));
            [n,e]=histcounts(dist,-prox_dist:dr:prox_dist);
            c=e(1:end-1)+dr/2;
            deltaA=pi*(2*prox_dist)*dr*(2*abs(c)+dr);
            n=n./deltaA;
            subplot(2,3,6);
            plot(c,n,'LineWidth',2,'Color',probecolour(probeidx));
            
            [rz_ch{probeidx},~]=hist3(rzcoord{probeidx},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
        end
        subplot(2,3,4);
        axis 'square';grid minor;xlabel('r');ylabel('z');
        title(sprintf('w_{synapse} from cluster %g s.f. = %gnm',nsig,1e3*abs(diff(rminmax))));
        set(gca,'Tag','synapse_rz_scatter');
        subplot(2,3,6);
        grid minor;xlabel('dist_{rz}');ylabel('\rho_{localisation}(\mum^{-3})');
        % parallel transform combined rgbimage
        subplot(2,3,5);set(gca,'Tag','synrgbimg');
        rgbimg=cat(3,rz_ch{rgbidx(1)}./max(rz_ch{rgbidx(1)}(:)),rz_ch{rgbidx(2)}./max(rz_ch{rgbidx(2)}(:)),rz_ch{rgbidx(3)}./max(rz_ch{rgbidx(3)}(:)));
        himage=image(imrotate(rgbimg,90));axis square;title(sprintf('dr=%g',dr));
        set(himage,'UserData',dr);
        
        
        
        
        
        
        
        %-----------------
        
        %{
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
            plot(e,n,'LineWidth',2,'Color',probecolour(probeidx));
        end
        grid minor;
        set(gca,'YScale','linear');
        xlabel('r_{cylindrical}');ylabel('\rho_{localisation}(\mum^{-2})');
        
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
            plot(e,n,'LineWidth',2,'Color',probecolour(probeidx));
        end
        grid minor;
        set(gca,'YScale','linear');
        xlabel('z_{cylindrical}');ylabel('\rho_{localisation}(\mum^{-2})');
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
            scatter(rz_tot{probeidx}(:,1),rz_tot{probeidx}(:,2),5,probecolour(probeidx),'.');
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
            plot(e,n,'LineWidth',2,'Color',probecolour(probeidx));
        end
        grid minor;
        set(gca,'YScale','linear');
        xlabel('Dist_{rz}');ylabel('\rho_{localisation}(\mum^{-2})');
        
        %--------------------
        % plot distance to the within synapse width
        subplot(2,3,5);hold on;
        for probeidx=1:numel(rz_tot)
        %{
            syndist=[];rzdist=[];
            syndist(:,1)=sqrt(sum(bsxfun(@minus,rz_tot{probeidx},[0,zminmax(1)]).^2,2));
            syndist(:,2)=sqrt(sum(bsxfun(@minus,rz_tot{probeidx},[0,zminmax(2)]).^2,2));
            syndist(:,3)=sqrt(sum(rz_tot{probeidx}.^2,2));
            rzdist=min(syndist,[],2).*sign(rz_tot{probeidx}(:,1));
            [n,e]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
            e=e(1:end-1)+dr/2;
            plot(e,n,'LineWidth',2,'Color',probecolour(probeidx));
        %}
            % central
            dist=rz_tot{probeidx}((rz_tot{probeidx}(:,2)<max(zminmax)&rz_tot{probeidx}(:,2)>min(zminmax)),1);
            [n,e]=histcounts(dist,min(dist):dr:max(dist));
            e=e(1:end-1)+dr/2;
            plot(e,n,'LineWidth',2,'Color',probecolour(probeidx));
        end
        title(sprintf('w_{synapse} from cluster %g s.f. = %gnm',nsig,1e3*abs(diff(zminmax))));
        grid minor;
        set(gca,'YScale','linear');
        xlabel('Dist_{synapse}');ylabel('N_{localisation}');
        
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
        xlabel('r');ylabel('z');
        set(himage,'UserData',dr);
        % print out dr value
        title(sprintf('dr=%g',dr));
        %}
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
                    case {'synapse_rz_scatter'}
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
                            %{
                                                        syndist=[];rzdist=[];
                             syndist(:,1)=sqrt(sum(bsxfun(@minus,rz_tot{pidx},[0,zminmax(1)]).^2,2));
                            syndist(:,2)=sqrt(sum(bsxfun(@minus,rz_tot{pidx},[0,zminmax(2)]).^2,2));
                            syndist(:,3)=sqrt(sum(rz_tot{pidx}.^2,2));
                            rzdist=min(syndist,[],2).*sign(rz_tot{pidx}(:,1));
                             [n,e]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
                            %deltaA=pi*dr*(2*abs(e(1:end-1))+dr);
                            %n=n./deltaA;
                            e=e(1:end-1)+dr/2;
                            set(lineplot.Children(numel(rz_tot)+1-pidx),'XData',e,'YData',n);
                            %}
                            % central
                            dist_sec=rz_tot{pidx}((rz_tot{pidx}(:,2)<max(zminmax)&rz_tot{pidx}(:,2)>min(zminmax)),1);
                            [n,e]=histcounts(dist_sec,min(dist_sec):dr:max(dist_sec));
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