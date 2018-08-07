function [ output_args ] = combine_projection(pathname,probe_colours)
%COMBINE_PROJECTION Summary of this function goes here
%   Detailed explanation goes here

if ~exist('pathname','var')
    pathname='D:\Temporary\Janosh\analysed';
end
if ~exist('probe_colours','var')
    probe_colours='brg';
    rgbidx=[find(probe_colours=='r'),find(probe_colours=='g'),find(probe_colours=='b')];
end
minres=0.005;%minimum resolution
maxres=0.5;%maximum resolution
dres=0.005;%delta resolution
FWHM_ratio=max(min(0.5,1),0);%FWHM percentage of max bounded by (0,1)

[filename,pathname,~]=uigetfile({'*.mat','MAT-files (*.mat)';'*.*','All Files (*.*)'},...
    'Select exported mat file folder','MultiSelect','on',pathname);
if ischar(pathname)
    if ~iscell(filename)
        filename={filename};
    end
    rz_tot=cell(1,3);
    for fileidx=1:numel(filename)
        temp=load(cat(2,pathname,filename{fileidx}),'-mat');
        rz_tot{1}=[rz_tot{1};cell2mat(temp.data(:,1)')];
        rz_tot{2}=[rz_tot{2};cell2mat(temp.data(:,2)')];
        rz_tot{3}=[rz_tot{3};cell2mat(temp.data(:,3)')];
    end
    dr=minres;dz=dr;zminmax=[];
    
    %combined data plot
    fh_tot=figure('Name',sprintf('Nearest Probe Site Distance to synapse collective (pageup/pagedown adjust dr)'),...
        'NumberTitle','off',...
        'MenuBar','none',...
        'ToolBar','figure',...
        'Position',[0,0,900,600],...
        'Color',[0.5,0.5,0.5],...
        'Tag','Synapse_tot',...
        'Keypressfcn',@figure_keypress);
    fh_tot.UserData.dr=dr;
    fh_tot.UserData.zminmax=zminmax;
    
    subplot(2,3,2);cla;
    [n1,e1]=histcounts(rz_tot{1}(:,1),min(rz_tot{1}(:,1)):dr:max(rz_tot{1}(:,1)));
    e1=e1(1:end-1)+dr/2;
    [n2,e2]=histcounts(rz_tot{2}(:,1),min(rz_tot{2}(:,1)):dr:max(rz_tot{2}(:,1)));
    e2=e2(1:end-1)+dr/2;
    [n3,e3]=histcounts(rz_tot{3}(:,1),min(rz_tot{3}(:,1)):dr:max(rz_tot{3}(:,1)));
    e3=e3(1:end-1)+dr/2;
    plot(e1,n1,'LineWidth',2,'Color',probe_colours(1));hold on;
    plot(e2,n2,'LineWidth',2,'Color',probe_colours(2));
    plot(e3,n3,'LineWidth',2,'Color',probe_colours(3));
    grid minor;xlabel('r');ylabel('N_{loc}');
    
    subplot(2,3,3);cla;
    [n1,e1]=histcounts(rz_tot{1}(:,2),min(rz_tot{1}(:,2)):dz:max(rz_tot{1}(:,2)));
    e1=e1(1:end-1)+dr/2;
    [n2,e2]=histcounts(rz_tot{2}(:,2),min(rz_tot{2}(:,2)):dz:max(rz_tot{2}(:,2)));
    e2=e2(1:end-1)+dr/2;
    [n3,e3]=histcounts(rz_tot{3}(:,2),min(rz_tot{3}(:,2)):dz:max(rz_tot{3}(:,2)));
    e3=e3(1:end-1)+dr/2;
    plot(e1,n1,'LineWidth',2,'Color',probe_colours(1));hold on;
    plot(e2,n2,'LineWidth',2,'Color',probe_colours(2));
    plot(e3,n3,'LineWidth',2,'Color',probe_colours(3));
    psdidx1=find(n2>=max(n2)*FWHM_ratio);
    psdidx2=find(n3>=max(n3)*FWHM_ratio);
    psdidx=[min([psdidx1,psdidx2]),max([psdidx1,psdidx2])];
    zminmax=[e3(psdidx(1));e3(psdidx(2))];
    fh_tot.UserData.zminmax=zminmax;
    grid minor;xlabel('z');ylabel('N_{loc}');
    
    subplot(2,3,1);cla;set(gca,'Tag','Synapse_tot_scatter');
    scatter([0;0],zminmax,200,'k','.');hold on;
    scatter(rz_tot{1}(:,1),rz_tot{1}(:,2),5,probe_colours(1),'.');
    scatter(rz_tot{2}(:,1),rz_tot{2}(:,2),5,probe_colours(2),'.');
    scatter(rz_tot{3}(:,1),rz_tot{3}(:,2),5,probe_colours(3),'.');
    axis 'square';grid minor;xlabel('r');ylabel('z');
    
    subplot(2,3,4);set(gca,'Tag','synrgbimg');
    scalemin=min([min(rz_tot{1},[],1);min(rz_tot{2},[],1);min(rz_tot{3},[],1)],[],1);
    scalemax=max([max(rz_tot{1},[],1);max(rz_tot{2},[],1);max(rz_tot{3},[],1)],[],1);
    [rz_ch{1},~]=hist3(rz_tot{1},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
    [rz_ch{2},~]=hist3(rz_tot{2},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
    [rz_ch{3},~]=hist3(rz_tot{3},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
    rgbimg=cat(3,rz_ch{rgbidx(1)}./max(rz_ch{rgbidx(1)}(:)),rz_ch{rgbidx(2)}./max(rz_ch{rgbidx(2)}(:)),rz_ch{rgbidx(3)}./max(rz_ch{rgbidx(3)}(:)));
    himage=image(imrotate(rgbimg,90));axis square;title(sprintf('dr=%g',dr));
    set(himage,'UserData',dr);
    
    subplot(2,3,5);cla;
    rzdist=sqrt(sum(rz_tot{1}.^2,2)).*sign(rz_tot{1}(:,1));
    [n1,e1]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
    e1=e1(1:end-1)+dr/2;
    rzdist=sqrt(sum(rz_tot{2}.^2,2)).*sign(rz_tot{2}(:,1));
    [n2,e2]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
    e2=e2(1:end-1)+dr/2;
    rzdist=sqrt(sum(rz_tot{3}.^2,2)).*sign(rz_tot{3}(:,1));
    [n3,e3]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
    e3=e3(1:end-1)+dr/2;
    plot(e1,n1,'LineWidth',2,'Color',probe_colours(1));hold on;
    plot(e2,n2,'LineWidth',2,'Color',probe_colours(2));
    plot(e3,n3,'LineWidth',2,'Color',probe_colours(3));
    grid minor;xlabel('rzdist');ylabel('N_{loc}');
    
    subplot(2,3,6);hold on;
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
    title(sprintf('W_{synapse} estimated from cluster FWHM = %gnm',1e3*abs(diff(zminmax))));
    grid minor;xlabel('synapsedist');ylabel('N_{loc}');
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
                        subplot(2,3,1);
                        temp=gca;
                        rz_tot{3}=[temp.Children(1).XData(:),temp.Children(1).YData(:)];
                        rz_tot{2}=[temp.Children(2).XData(:),temp.Children(2).YData(:)];
                        rz_tot{1}=[temp.Children(3).XData(:),temp.Children(3).YData(:)];
                        scalemin=min([min(rz_tot{1},[],1);min(rz_tot{2},[],1);min(rz_tot{3},[],1)],[],1);
                        scalemax=max([max(rz_tot{1},[],1);max(rz_tot{2},[],1);max(rz_tot{3},[],1)],[],1);
                        dz=dr;
                        [rz_ch{1},~]=hist3(rz_tot{1},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
                        [rz_ch{2},~]=hist3(rz_tot{2},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
                        [rz_ch{3},~]=hist3(rz_tot{3},{scalemin(1):dr:scalemax(1),scalemin(2):dz:scalemax(2)});
                        rgbimg=cat(3,rz_ch{rgbidx(1)}./max(rz_ch{rgbidx(1)}(:)),rz_ch{rgbidx(2)}./max(rz_ch{rgbidx(2)}(:)),rz_ch{rgbidx(3)}./max(rz_ch{rgbidx(3)}(:)));
                        subplot(2,3,3);
                        [n1,e1]=histcounts(rz_tot{1}(:,2),min(rz_tot{1}(:,2)):dz:max(rz_tot{1}(:,2)));
                        e1=e1(1:end-1)+dr/2;
                        [n2,e2]=histcounts(rz_tot{2}(:,2),min(rz_tot{2}(:,2)):dz:max(rz_tot{2}(:,2)));
                        e2=e2(1:end-1)+dr/2;
                        [n3,e3]=histcounts(rz_tot{3}(:,2),min(rz_tot{3}(:,2)):dz:max(rz_tot{3}(:,2)));
                        e3=e3(1:end-1)+dr/2;
                        lineplot=get(gca);
                        set(lineplot.Children(1),'XData',e3,'YData',n3);
                        set(lineplot.Children(2),'XData',e2,'YData',n2);
                        set(lineplot.Children(3),'XData',e1,'YData',n1);
                        psdidx1=find(n2>=max(n2)*FWHM_ratio);
                        psdidx2=find(n3>=max(n3)*FWHM_ratio);
                        psdidx=[min([psdidx1,psdidx2]),max([psdidx1,psdidx2])];
                        zminmax=[e3(psdidx(1));e3(psdidx(2))]
                        hplot.UserData.zminmax=zminmax;
                        subplot(2,3,4);
                        himage=get(gca,'Children');
                        set(himage,'CData',imrotate(rgbimg,90));axis tight;
                        set(himage,'UserData',dr);
                        title(sprintf('dr=%g',dr));
                        subplot(2,3,2);
                        [n1,e1]=histcounts(rz_tot{1}(:,1),min(rz_tot{1}(:,1)):dr:max(rz_tot{1}(:,1)));
                        e1=e1(1:end-1)+dr/2;
                        [n2,e2]=histcounts(rz_tot{2}(:,1),min(rz_tot{2}(:,1)):dr:max(rz_tot{2}(:,1)));
                        e2=e2(1:end-1)+dr/2;
                        [n3,e3]=histcounts(rz_tot{3}(:,1),min(rz_tot{3}(:,1)):dr:max(rz_tot{3}(:,1)));
                        e3=e3(1:end-1)+dr/2;
                        lineplot=get(gca);
                        set(lineplot.Children(1),'XData',e3,'YData',n3);
                        set(lineplot.Children(2),'XData',e2,'YData',n2);
                        set(lineplot.Children(3),'XData',e1,'YData',n1);
                        subplot(2,3,5);
                        rzdist=sqrt(sum(rz_tot{1}.^2,2)).*sign(rz_tot{1}(:,1));
                        [n1,e1]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
                        e1=e1(1:end-1)+dr/2;
                        rzdist=sqrt(sum(rz_tot{2}.^2,2)).*sign(rz_tot{2}(:,1));
                        [n2,e2]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
                        e2=e2(1:end-1)+dr/2;
                        rzdist=sqrt(sum(rz_tot{3}.^2,2)).*sign(rz_tot{3}(:,1));
                        [n3,e3]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
                        e3=e3(1:end-1)+dr/2;
                        lineplot=get(gca);
                        set(lineplot.Children(1),'XData',e3,'YData',n3);
                        set(lineplot.Children(2),'XData',e2,'YData',n2);
                        set(lineplot.Children(3),'XData',e1,'YData',n1);
                        subplot(2,3,6);
                        lineplot=get(gca);
                        for pidx=1:numel(rz_tot)
                            syndist=[];rzdist=[];
                            syndist(:,1)=sqrt(sum(bsxfun(@minus,rz_tot{pidx},[0,zminmax(1)]).^2,2));
                            syndist(:,2)=sqrt(sum(bsxfun(@minus,rz_tot{pidx},[0,zminmax(2)]).^2,2));
                            syndist(:,3)=sqrt(sum(rz_tot{pidx}.^2,2));
                            rzdist=min(syndist,[],2).*sign(rz_tot{pidx}(:,1));
                            [n,e]=histcounts(rzdist,min(rzdist):dr:max(rzdist));
                            e=e(1:end-1)+dr/2;
                            set(lineplot.Children(numel(rz_tot)+1-pidx),'XData',e,'YData',n);
                        end
                        title(sprintf('W_{synapse} estimated from cluster FWHM = %gnm',1e3*abs(diff(zminmax))));
                        grid minor;xlabel('synapsedist');ylabel('N_{loc}');
                end
        end
    end
end

