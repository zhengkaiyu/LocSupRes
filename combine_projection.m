function combine_projection(varargin)
%COMBINE_PROJECTION combine exported 3probe 2d projection scatter data
% input options: pathname = folder containing exported files
%                probecolour = colour code for plotting, default 'grbymc'
%                synapsepair = index of pre and post probe id, default [1,2]
%                minres,maxres,dres = resolution settings, default 0.005,0.5,0.005 um
%                nsig = num of significant figures for estimating synapse, default 1.0
%
% exported file format: [cluster_pts,TM,prox_dist] containing cluster
%                       points before transformation, transformation matrix
%                       for [origin, synapse plane normal, symmetry plane
%                       normal p2, orthogonal plane p3 ] and proximity
%                       distance which is generally 1um (Note this is for
%                       cylindrical coordinates only, spherical coordinate
%                       cutoff has already been applied as sqrt(2)*R before
%                       exporting
% unit um
%------------------------------------------------------------------------
% read argument in options
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
if ~exist('probecolour','var')
    % default colour coding
    probecolour='grbymc';
    % rgb index in probe colour for rgb image formation
    rgbidx=[find(probecolour=='r'),find(probecolour=='g'),find(probecolour=='b')];
end
if ~exist('synapsepair','var')
    % synapse probe pair index, default to probe 2 and 3
    synapsepair=[1,2];
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
% subplot grid size row for transform and col for parameters
plotsize=[3,5];
n_transform=2;
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
                    % TM transform matrix in format
                    %[origin,synapse plane normal,symmetry plane normal p2,orthogonal plane p3]
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
        dr=dres;dz=dr;
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
        aligedcollect=reshape([aligned_cluster_pts{:}]',n_probe,n_cluster)';
        
        %1st transform = parallel transform
        %2nd transform = orthogonal transform
        for transidx=0:n_transform
            % set plots to hold on, except rgb image
            for subplotidx=1:plotsize(2)
                if subplotidx~=2
                    subplot(plotsize(1),plotsize(2),subplotidx+transidx*plotsize(2));cla;hold on;
                end
            end
            if transidx==0
                % original location
                % all original cluster xyz loc for the probe aligned
                for probeidx=1:n_probe
                    xyzcoord{probeidx}=cell2mat(aligedcollect(:,probeidx));
                    % only get the ones in the cube
                    xyzcoord{probeidx}=xyzcoord{probeidx}((abs(xyzcoord{probeidx}(:,1))<=prox_dist)&(abs(xyzcoord{probeidx}(:,2))<=prox_dist)&(abs(xyzcoord{probeidx}(:,3))<=prox_dist),1:3);
                end
                % work out synapse edge using parallel transform
                zsigma=[nsig*std(xyzcoord{synapsepair(1)}(:,3));nsig*std(xyzcoord{synapsepair(2)}(:,3))];% in z axis
                ysigma=[nsig*std(xyzcoord{synapsepair(1)}(:,2));nsig*std(xyzcoord{synapsepair(2)}(:,2))];% in y axis
                zsminmax=max(zsigma(:));
                ysminmax=max(ysigma(:));
                scalemin=min(cell2mat(xyzcoord'));
                scalemax=max(cell2mat(xyzcoord'));
                for probeidx=1:n_probe
                    subplot(plotsize(1),plotsize(2),1+transidx*plotsize(2));
                    scatter3(xyzcoord{probeidx}(:,1),xyzcoord{probeidx}(:,2),xyzcoord{probeidx}(:,3),1,probecolour(probeidx),'.');
                    % rgb colour channel
                    [xy_ch{probeidx},~]=hist3(xyzcoord{probeidx}(:,1:2),{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
                    
                    %preserving x dim so that pre and post side are
                    %maintained!! remember this for subsequent calculations
                    % x distance
                    dist=xyzcoord{probeidx}(:,1);
                    [n,e]=histcounts(dist,[fliplr(0:-dr:-prox_dist),dr:dr:prox_dist]);
                    c=e(1:end-1)+dr/2;
                    %correction for rectangular slab volume
                    deltaA=dr*(2*prox_dist)^2;
                    n=n./deltaA;
                    subplot(plotsize(1),plotsize(2),3+transidx*plotsize(2));
                    plot(c,n,'LineWidth',2,'Color',probecolour(probeidx));
                    
                    % radial distance | xyz distance to origin
                    dist=sqrt(sum(xyzcoord{probeidx}.^2,2)).*sign(sign(xyzcoord{probeidx}(:,1))+0.5);
                    [n,e]=histcounts(dist,[fliplr(0:-dr:-prox_dist),dr:dr:prox_dist]);
                    c=e(1:end-1)+dr/2;
                    % correctin for spherical shell volume, but only half
                    % shell for polarity preservation in x
                    deltaA=0.5*(4/3)*pi*dr*(3*c.^2+3*dr*abs(c)+dr^2);
                    n=n./deltaA;
                    subplot(plotsize(1),plotsize(2),4+transidx*plotsize(2));
                    plot(c,n,'LineWidth',2,'Color',probecolour(probeidx));
                    
                    % edge distance (depending on sminmax
                    % cutoff in r(sminmax)
                    % calculate radial distance to edge from selected points
                    % histogram of edge distance
                    subplot(plotsize(1),plotsize(2),5+transidx*plotsize(2));
                    %plot(c,n,'LineWidth',2,'Color',probecolour(probeidx));
                end
                subplot(plotsize(1),plotsize(2),1+transidx*plotsize(2));
                % plot estimated synapse
                fz=@(t)zsminmax*sin(t);fy=@(t)ysminmax*cos(t);fx=@(t)0*t;
                t = linspace(0,2*pi,101);
                line(fx(t),fy(t),fz(t),'Color','k','LineWidth',5);
                %set plot scatter properties
                axis 'square';grid minor;xlabel('x');ylabel('y');zlabel('z');
                title(sprintf('w_{synapse} from cluster %g s.f. = %g|%gnm',nsig,1e3*zsminmax*2,1e3*ysminmax*2));
                set(gca,'Tag','synapse_xyz_scatter');
                % transform combined rgbimage
                subplot(plotsize(1),plotsize(2),2+transidx*plotsize(2));
                rgbimg=cat(3,xy_ch{rgbidx(1)}./max(xy_ch{rgbidx(1)}(:)),xy_ch{rgbidx(2)}./max(xy_ch{rgbidx(2)}(:)),xy_ch{rgbidx(3)}./max(xy_ch{rgbidx(3)}(:)));
                himage=image(imrotate(rgbimg,90));axis square;title(sprintf('dr=%g',dr));
                set(gca,'XTickLabel',[],'YTickLabel',[]);
                xlabel('x');ylabel('y');
                set(gca,'Tag','synrgbimg');
                
                %set plot r distance plot properties
                subplot(plotsize(1),plotsize(2),3+transidx*plotsize(2));
                % plot synpase size for illustration
                plot([-zsminmax,zsminmax],[100,100],'k-','LineWidth',5);
                set(gca,'YScale','linear');
                axis square;grid minor;xlabel('x_{dist}');ylabel('\rho_{localisation}(\mum^{-3})');
                set(gca,'Tag','synapse_z_hist');
                
                %set plot r distance plot properties
                subplot(plotsize(1),plotsize(2),4+transidx*plotsize(2));
                plot([-ysminmax,ysminmax],[100,100],'k-','LineWidth',5);
                set(gca,'YScale','linear');
                axis square;grid minor;xlabel('xyz_{dist}');ylabel('\rho_{localisation}(\mum^{-3})');
                set(gca,'Tag','synapse_xyz_hist');
                
                %set plot edge distance plot properties
                subplot(plotsize(1),plotsize(2),5+transidx*plotsize(2));
                set(gca,'YScale','linear');
                axis square;grid minor;xlabel('edge_{dist}');ylabel('\rho_{localisation}(\mum^{-3})');
                set(gca,'Tag','synapse_edge_hist');
            else
                % transformed location
                rzcollect=reshape([trans_cluster_pts{:,transidx}]',n_probe,n_cluster)';
                for probeidx=1:n_probe
                    % all cluster transformed loc for the probe
                    rzcoord{probeidx}=cell2mat(rzcollect(:,probeidx)); %#ok<*AGROW>
                end
                scalemin=min(cell2mat(rzcoord'));
                scalemax=max(cell2mat(rzcoord'));
                % going through each probe and calculate and plot
                for probeidx=1:n_probe
                    subplot(plotsize(1),plotsize(2),1+transidx*plotsize(2));
                    scatter(rzcoord{probeidx}(:,1),rzcoord{probeidx}(:,2),1,probecolour(probeidx),'.');
                    
                    % rgb colour channel
                    [rz_ch{probeidx},~]=hist3(rzcoord{probeidx},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
                    
                    %preserving r dim so that pre and post side are
                    %maintained!! remember this for subsequent calculations
                    % r distance
                    dist=rzcoord{probeidx}(:,1);
                    [n,e]=histcounts(dist,[fliplr(0:-dr:-prox_dist),dr:dr:prox_dist]);
                    c=e(1:end-1)+dr/2;
                    %correction for theta collapse in theta/r/z coord, can be derived from cylindrical jacobian
                    %prox_dist is range(z)
                    %0.5 because of polarity preservation, so only half
                    deltaA=0.5*(2*pi*prox_dist)*dr*(dr+2*abs(c));
                    n=n./deltaA;
                    subplot(plotsize(1),plotsize(2),3+transidx*plotsize(2));
                    plot(c,n,'LineWidth',2,'Color',probecolour(probeidx));
                    
                    % radial distance | rz distance to origin
                    dist=sqrt(sum(rzcoord{probeidx}.^2,2)).*sign(sign(rzcoord{probeidx}(:,1))+0.5);
                    [n,e]=histcounts(dist,[fliplr(0:-dr:-prox_dist),dr:dr:prox_dist]);
                    c=e(1:end-1)+dr/2;
                    % correction
                    deltaA=0.5*(4/3)*pi*dr*(3*c.^2+3*dr*abs(c)+dr^2);%spherical
                    n=n./deltaA;
                    subplot(plotsize(1),plotsize(2),4+transidx*plotsize(2));
                    plot(c,n,'LineWidth',2,'Color',probecolour(probeidx));
                    
                    % edge distance (depending on sminmax
                    % cutoff in r(sminmax)
                    % calculate radial distance to edge from selected points
                    % histogram of edge distance
                    subplot(plotsize(1),plotsize(2),5+transidx*plotsize(2));
                    
                    %plot(c,n,'LineWidth',2,'Color',probecolour(probeidx));
                end
                %------------
                % scatter synapse edge points on scatter plot
                subplot(plotsize(1),plotsize(2),1+transidx*plotsize(2));
                % work out synapse edge using parallel transform
                switch transidx
                    case 1
                        sigma=[nsig*std(rzcoord{synapsepair(1)}(:,2));nsig*std(rzcoord{synapsepair(2)}(:,2))];
                        sminmax=max(sigma(:));
                        scatter([0;0],[-sminmax,sminmax],30,'k','s','MarkerFaceColor','k');
                    case 2
                        sigma=[nsig*std(rzcoord{synapsepair(1)}(:,1));nsig*std(rzcoord{synapsepair(2)}(:,1))];
                        sminmax=max(sigma(:));
                        scatter([-sminmax,sminmax],[0,0],30,'k','s','MarkerFaceColor','k');
                end
                %set plot scatter properties
                axis 'square';grid minor;xlabel('r');ylabel('z');
                title(sprintf('w_{synapse} from cluster %g s.f. = %gnm',nsig,1e3*sminmax*2));
                set(gca,'Tag','synapse_rz_scatter');
                
                % transform combined rgbimage
                subplot(plotsize(1),plotsize(2),2+transidx*plotsize(2));
                rgbimg=cat(3,rz_ch{rgbidx(1)}./max(rz_ch{rgbidx(1)}(:)),rz_ch{rgbidx(2)}./max(rz_ch{rgbidx(2)}(:)),rz_ch{rgbidx(3)}./max(rz_ch{rgbidx(3)}(:)));
                himage=image(imrotate(rgbimg,90));axis square;title(sprintf('dr=%g',dr));
                set(gca,'XTickLabel',[],'YTickLabel',[]);
                set(gca,'Tag','synrgbimg');
                
                %set plot r distance plot properties
                subplot(plotsize(1),plotsize(2),3+transidx*plotsize(2));
                % plot synpase size for illustration
                plot([-sminmax,sminmax],[100,100],'k-','LineWidth',5);
                set(gca,'YScale','linear');
                axis square;grid minor;xlabel('r_{dist}');ylabel('\rho_{localisation}(\mum^{-3})');
                set(gca,'Tag','synapse_r_hist');
                
                %set plot r distance plot properties
                subplot(plotsize(1),plotsize(2),4+transidx*plotsize(2));
                plot([-sminmax,sminmax],[100,100],'k-','LineWidth',5);
                set(gca,'YScale','linear');
                axis square;grid minor;xlabel('rz_{dist}');ylabel('\rho_{localisation}(\mum^{-3})');
                set(gca,'Tag','synapse_rz_hist');
                
                %set plot edge distance plot properties
                subplot(plotsize(1),plotsize(2),5+transidx*plotsize(2));
                set(gca,'YScale','linear');
                axis square;grid minor;xlabel('edge_{dist}');ylabel('\rho_{localisation}(\mum^{-3})');
                set(gca,'Tag','synapse_edge_hist');
            end
        end
        %-----------------
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
                    case {'synapse_xyz_scatter','synapse_rz_scatter',...
                            'synapse_r_hist','synapse_rz_hist',...
                            'synapse_x_hist','synapse_xyz_hist'}
                        keyidx=str2double(eventkey.Key(2));
                        if keyidx<=numel(hplot.Children)
                            switch hplot.Children(keyidx).Visible
                                case 'off'
                                    hplot.Children(keyidx).Visible='on';
                                case 'on'
                                    hplot.Children(keyidx).Visible='off';
                            end
                        end
                end
            case {'home','end','pagedown','pageup'}
                hplot=gcf;
                switch hplot.Tag
                    case 'Synapse_tot'
                        dr=hplot.UserData.dr;
                        %sminmax=hplot.UserData.sminmax;
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
                        dz=dr;
                        hplot.UserData.dr=dr;
                        for Tidx=0:n_transform
                            % get plot handle
                            scatterplot=subplot(plotsize(1),plotsize(2),1+Tidx*plotsize(2));
                            rgbimgplot=subplot(plotsize(1),plotsize(2),2+Tidx*plotsize(2));
                            subplot(plotsize(1),plotsize(2),3+Tidx*plotsize(2));
                            lineplot1=get(gca,'Children');lineplot1=flipud(lineplot1);
                            subplot(plotsize(1),plotsize(2),4+Tidx*plotsize(2));
                            lineplot2=get(gca,'Children');lineplot2=flipud(lineplot2);
                            subplot(plotsize(1),plotsize(2),5+Tidx*plotsize(2));
                            lineplot3=get(gca,'Children');lineplot3=flipud(lineplot3);
                            if Tidx==0
                                %get processed scatter points
                                for pidx=1:n_probe
                                    plotidx=n_probe-pidx+2;
                                    xyz_pts{pidx}=[scatterplot.Children(plotidx).XData(:),...
                                        scatterplot.Children(plotidx).YData(:),...
                                        scatterplot.Children(plotidx).ZData(:)];
                                end
                                % recalculate and plot
                                scalemin=min(cell2mat(xyz_pts'));
                                scalemax=max(cell2mat(xyz_pts'));
                                for pidx=1:n_probe
                                    plotidx=pidx;
                                    % rgb colour channel
                                    [xy_ch{pidx},~]=hist3(xyz_pts{pidx}(:,1:2),{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
                                    
                                    %preserving x dim so that pre and post side are
                                    %maintained!! remember this for subsequent calculations
                                    % x distance
                                    dist=xyz_pts{pidx}(:,1);
                                    [n,e]=histcounts(dist,[fliplr(0:-dr:-prox_dist),dr:dr:prox_dist]);
                                    c=e(1:end-1)+dr/2;
                                    %correction for rectangular slab volume
                                    deltaA=dr*(2*prox_dist)^2;
                                    n=n./deltaA;
                                    set(lineplot1(plotidx),'XData',c,'YData',n);
                                    
                                    % radial distance | xyz distance to origin
                                    dist=sqrt(sum(xyz_pts{pidx}.^2,2)).*sign(sign(xyz_pts{pidx}(:,1))+0.5);
                                    [n,e]=histcounts(dist,[fliplr(0:-dr:-prox_dist),dr:dr:prox_dist]);
                                    c=e(1:end-1)+dr/2;
                                    % correctin for spherical shell volume, but only half
                                    % shell for polarity preservation in x
                                    deltaA=0.5*(4/3)*pi*dr*(3*c.^2+3*dr*abs(c)+dr^2);
                                    n=n./deltaA;
                                    set(lineplot2(plotidx),'XData',c,'YData',n);
                                    
                                    % edge distance (depending on sminmax
                                    % cutoff in r(sminmax)
                                    % calculate radial distance to edge from selected points
                                    % histogram of edge distance
                                    
                                    % set(lineplot3(plotidx),'XData',c,'YData',n);
                                end
                                % transform combined rgbimage
                                rgbimg=cat(3,xy_ch{rgbidx(1)}./max(xy_ch{rgbidx(1)}(:)),xy_ch{rgbidx(2)}./max(xy_ch{rgbidx(2)}(:)),xy_ch{rgbidx(3)}./max(xy_ch{rgbidx(3)}(:)));
                                set(rgbimgplot.Children,'CData',imrotate(rgbimg,90));
                                title(rgbimgplot,sprintf('dr=%g',dr));
                                axis(rgbimgplot,'tight');
                            else
                                % get transformed scatter points
                                %rz_pts{:,transidx}=[temp.Children(n_probe-pidx+1).XData(:),temp.Children(n_probe-pidx+1).YData(:)];
                                %get processed scatter points
                                for pidx=1:n_probe
                                    plotidx=n_probe-pidx+2;
                                    rz_pts{pidx}=[scatterplot.Children(plotidx).XData(:),...
                                        scatterplot.Children(plotidx).YData(:)];
                                end
                                % recalculate and plot
                                scalemin=min(cell2mat(rz_pts'));
                                scalemax=max(cell2mat(rz_pts'));
                                for pidx=1:n_probe
                                    plotidx=pidx;
                                    % rgb colour channel
                                    [rz_ch{pidx},~]=hist3(rz_pts{pidx},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
                                    
                                    %preserving r dim so that pre and post side are
                                    %maintained!! remember this for subsequent calculations
                                    % r distance
                                    dist=rz_pts{pidx}(:,1);
                                    [n,e]=histcounts(dist,[fliplr(0:-dr:-prox_dist),dr:dr:prox_dist]);
                                    c=e(1:end-1)+dr/2;
                                    %correction for theta collapse in theta/r/z coord, can be derived from cylindrical jacobian
                                    %prox_dist is range(z)
                                    %0.5 because of polarity preservation, so only half
                                    deltaA=0.5*(2*pi*prox_dist)*dr*(dr+2*abs(c));
                                    n=n./deltaA;
                                    set(lineplot1(plotidx),'XData',c,'YData',n);
                                    
                                    % radial distance | rz distance to origin
                                    dist=sqrt(sum(rz_pts{pidx}.^2,2)).*sign(sign(rz_pts{pidx}(:,1))+0.5);
                                    [n,e]=histcounts(dist,[fliplr(0:-dr:-prox_dist),dr:dr:prox_dist]);
                                    c=e(1:end-1)+dr/2;
                                    % correction
                                    deltaA=0.5*(4/3)*pi*dr*(3*c.^2+3*dr*abs(c)+dr^2);
                                    n=n./deltaA;
                                    set(lineplot2(plotidx),'XData',c,'YData',n);
                                    
                                    % edge distance (depending on sminmax
                                    % cutoff in r(sminmax)
                                    % calculate radial distance to edge from selected points
                                    % histogram of edge distance
                                    % set(lineplot3(plotidx),'XData',c,'YData',n);
                                end
                                % transform combined rgbimage
                                rgbimg=cat(3,rz_ch{rgbidx(1)}./max(rz_ch{rgbidx(1)}(:)),rz_ch{rgbidx(2)}./max(rz_ch{rgbidx(2)}(:)),rz_ch{rgbidx(3)}./max(rz_ch{rgbidx(3)}(:)));
                                set(rgbimgplot.Children,'CData',imrotate(rgbimg,90));
                                title(rgbimgplot,sprintf('dr=%g',dr));
                                axis(rgbimgplot,'tight');
                            end
                        end
                end
        end
    end
end