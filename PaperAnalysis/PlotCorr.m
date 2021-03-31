%Load azimuth/elevation maps
Azimuth_map='your path\mean_azimuth_map.npy';
Altitude_map='your path\mean_altitude_map.npy';

Azimuth=readNPY(Azimuth_map);
Altitude=readNPY(Altitude_map);
% Load tif
im_dir='';

% Filenames
filename_tif=['CCF_borders_10um' '.tif'];

    InfoImage.tif=imfinfo([im_dir filesep filename_tif]);
    im_tif=zeros(InfoImage.tif.Height,InfoImage.tif.Width,length(InfoImage.tif),'double');

    TifLink = Tiff([im_dir filesep filename_tif], 'r');
    for i=1:length(InfoImage.tif)
       TifLink.setDirectory(i);
       im_tif(:,:,i)=TifLink.read();
    end
    TifLink.close();

%% Remove outer area
Az=Azimuth;
Az(Az==1000) = nan;
caxis([-30 90]);
% Al=Altitude;
% Al(Al==1000) = nan;
% caxis([-30 30]);

%%%Load 3d brain (top-down view)%%%

%% Overlay azimuth
Az_map=flipdim(Az,1);
% Az_map=imrotate(flip,90);

hold on
Az_m=image(Az_map,'XData',0,'YData',-1100,'AlphaData',0.9);
% im=image(Az_map,'AlphaData',0.8);
% contourcmap('jet');
hold off

% flip=flipdim(Al,2);
% Al_map=imrotate(flip,90);
% 
% hold on
% image(Al_map,'XData',500,'YData',0,'AlphaData',0.8);
% % contourcmap('jet');
% hold off

%% Overlay CCF boundaries 
hold on
% CCF=imrotate(im_tif, 90);
im=image(im_tif,'XData',0,'YData',0);
im.AlphaData=max(im_tif,[],3);
hold off
%% Remove outer area
Az=Azimuth;
Az(Az==1000) = nan;
caxis([-30 90]);
% Al=Altitude;
% Al(Al==1000) = nan;
% caxis([-30 30]);

%%Creat blank matrix 1320x1140
blank=zeros(size(im_tif));
% blank(501:1100,1:600)=Az;
blank(501:1100,1:600)=Az;

% contourcmap('jet',[-30:5:30]); % Al

contourcmap('jet',[-30:5:90]); % Az

%% Getting injection coords
% load data tables
% load injection coords

% HVAs only
inj_order_hva=table2array(data_hva(1,3:end));% exclude first injection
areas_hva=table2array(data_hva(2:end,1));
ml_hva=abs(Inj_coords(inj_order_hva,2))'; %take abs to make all same side
ap_hva=Inj_coords(inj_order_hva,1)';

% Whole brain (excluding HVAs)
inj_order_wb=table2array(data_wholebrain(1,2:end));
areas_wb=table2array(data_wholebrain(2:22,1)); % exclude HVAs
ml_wb=abs(Inj_coords(inj_order_wb,2))'; %take abs to make all same side
ap_wb=Inj_coords(inj_order_wb,1)';

% Whole brain grouped(excluding HVAs)
areas_wb_gr=table2array(data_wb_gr(2:end,1)); % exclude HVAs


%% Getting retinotopy coords
%%Creat blank matrix 1320x1140
allen_map_azi=nan(1320,1140);
% blank(501:1100,1:600)=Az;
azi_im(azi_im==1000) = nan;
allen_map_azi(501:1100,1:600)=azi_im;
figure;
imagesc(allen_map_azi)
title('azimuth map')

%%Creat blank matrix 1320x1140
allen_map_alt=nan(1320,1140);
% blank(501:1100,1:600)=Az;
alt_im(alt_im==1000) = nan;
allen_map_alt(501:1100,1:600)=alt_im;
figure;
imagesc(allen_map_alt)
title('altitude map')

addpath(genpath('your path\GitHub\allenCCF'))
bregma=allenCCFbregma();
atlas_resolution=0.01;

%%
allen_ml=-abs(Inj_coords(2:end,2))/atlas_resolution+bregma(3);% exclude first injection
allen_ap=-Inj_coords(2:end,1)/atlas_resolution+bregma(1);% exclude first injection
allen_ml_hva=-ml_hva/atlas_resolution+bregma(3);
allen_ap_hva=-ap_hva/atlas_resolution+bregma(1);
allen_ml_wb=-ml_wb/atlas_resolution+bregma(3);
allen_ap_wb=-ap_wb/atlas_resolution+bregma(1);

figure;
imagesc(allen_map)
hold on
sz=60;
h1=scatter(allen_ml, allen_ap,sz,'ko');
hold on
sz=20;
h2=scatter(allen_ml_hva, allen_ap_hva,sz,'ro');
hold on
sz=5;
h3=scatter(allen_ml_wb, allen_ap_wb,sz,'bo');
hold on
im=image(im_bound,'XData',0,'YData',0);
im.AlphaData=max(im_bound,[],3);
hold off
legend('all','HVA','whole brain')
caxis([-30 90]);
box off 
axis off equal
title ('azimuth ret')
contourcmap('jet',[-30:5:90]); % Az

% get azimuth retinotopy
for i=1:length(allen_ap_hva)
azi_ret_hva(i)=allen_map_azi(allen_ap_hva(i), allen_ml_hva(i));
end

for i=1:length(allen_ap_wb)
azi_ret_wb(i)=allen_map_azi(allen_ap_wb(i), allen_ml_wb(i));
end

% get altitude retinotopy
for i=1:length(allen_ap_hva)
alt_ret_hva(i)=allen_map_alt(allen_ap_hva(i), allen_ml_hva(i));
end

for i=1:length(allen_ap_wb)
alt_ret_wb(i)=allen_map_alt(allen_ap_wb(i), allen_ml_wb(i));
end

%% HVA only
% HVA azimuth retinotopy
figure;
k=1;
for i=1:9
    subplot(3,3,k)
    area_data=table2array(data_hva(k+1,3:end));% exclude first injection
    plot(azi_ret_hva,area_data,'o', 'MarkerSize',5,'MarkerEdgeColor',[0.5 0.5 0.5]);
    h = lsline;
    set(h,'color','r')
    [r_hva_azi(i),p_hva_azi(i)]=corr(azi_ret_hva',area_data', 'Rows','complete');
    title({[areas_hva{i} '  r=' num2str(round(r,2))],['p=' num2str(round(p,4))]})
    xlabel('azimuth retinotopy (deg)')
    ylabel('cell count (%)')
    xlim([0 60])
    set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.02 0.05], 'color', 'none', 'box','off')
%     offsetAxes
    k=k+1;
end

% HVA altitude retinotopy
figure;
k=1;
for i=1:9
    subplot(3,3,k)
    area_data=table2array(data_hva(k+1,3:end));% exclude first injection
    plot(alt_ret_hva,area_data,'o', 'MarkerSize',5,'MarkerEdgeColor',[0.5 0.5 0.5]);
    h = lsline;
    set(h,'color','r')
    [r_hva_alt(i),p_hva_alt(i)]=corr(alt_ret_hva',area_data', 'Rows','complete');
    title({[areas_hva{i} '  r=' num2str(round(r,2))],['p=' num2str(round(p,4))]})
    xlabel('altitude retinotopy (deg)')
    ylabel('cell count (%)')
    xlim([-10 15])
    set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.02 0.05], 'color', 'none', 'box','off')
%     offsetAxes
    k=k+1;
end

% Wholebrain azimuth retinotopy
figure;
k=1;
for i=1:15
    subplot(5,3,k)
    area_data=table2array(data_wb_gr(k+1,2:end));
    plot(azi_ret_wb,area_data,'o', 'MarkerSize',5,'MarkerEdgeColor',[0.5 0.5 0.5]);
    h = lsline;
    set(h,'color','r')
    [r_wb_azi(i),p_wb_azi(i)]=corr(azi_ret_wb',area_data','Rows','complete');
    title({[areas_wb_gr{i} '  r=' num2str(round(r,2))],['p=' num2str(round(p,4))]})
    xlabel('azimuth retinotopy (deg)')
    ylabel('cell count (%)')
    box off
    xlim([0 60])
    set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.02 0.05], 'color', 'none', 'box','off')
    k=k+1; 
end

% Wholebrain altitude retinotopy
figure;
k=1;
for i=1:15
    subplot(5,3,k)
    area_data=table2array(data_wb_gr(k+1,2:end));
    plot(alt_ret_wb,area_data,'o', 'MarkerSize',5,'MarkerEdgeColor',[0.5 0.5 0.5]);
    h = lsline;
    set(h,'color','r')
    [r_wb_alt(i),p_wb_alt(i)]=corr(alt_ret_wb',area_data','Rows','complete');
    title({[areas_wb_gr{i} '  r=' num2str(round(r,2))],['p=' num2str(round(p,4))]})
    xlabel('altitude retinotopy (deg)')
    ylabel('cell count (%)')
    box off
    xlim([-10 15])
    set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.02 0.05], 'color', 'none', 'box','off')
    k=k+1;
end

%% Barplot for correlation
% table of r and p
% bar plot of r and p
r_azi=[r_hva_azi,r_wb_azi];
r_alt=[r_hva_alt,r_wb_alt];
p_azi=[p_hva_azi,p_wb_azi];
p_alt=[p_hva_alt,p_wb_alt];
areas_all=[areas_hva; areas_wb_gr];

rp=table(areas_all ,r_azi',p_azi',r_alt',p_alt', 'VariableNames',{'Areas' 'Azi_r' 'Azi_p' 'Alt_r' 'Alt_p'});

figure;
subplot(1,2,1)
    barh(r_azi')
    ax=gca;
    yticks([1:24])    
    ax.YTickLabel = areas_all;
    set(gca,'TickDir','out');
    box off
    set(gca,'Ydir','reverse')
    xlim([-0.8 0.8])
    title('azimuth r')

subplot(1,2,2)
    barh(r_alt')
    ax=gca;
    yticks([1:24])    
    ax.YTickLabel = areas_all;
    set(gca,'TickDir','out');
    box off
    set(gca,'Ydir','reverse')
    xlim([-0.8 0.8])
    title('altitude r')
    

%% bubble plot
figure;
k=1;
for i=1:9
    subplot(3,3,k)
    area_data=table2array(data_hva(k+1 ,3:end));% exclude first injection
    area_data_norm=(area_data-min(area_data))/max(area_data);
    msize=40;
    scatter3(allen_ml_hva,allen_ap_hva,zeros(size(allen_ml_hva)),msize,area_data_norm,'filled','MarkerEdgeColor',[.5 .5 .5]);
    axis equal off 
    box off
    view(2)
    set(gca, 'YDir','reverse')
    xlim([30 570])
    ylim([650 1050])
    hold on
    im=image(im_bound,'XData',0,'YData',0);
    im.AlphaData=max(im_bound,[],3);
    title(areas_hva{i})  
    colormap cool
    k=k+1;
end

figure;
k=1;
for i=1:9
    subplot(3,3,k)
    area_data=table2array(data_hva(k+1,3:end));% exclude first injection
    area_data_norm=(area_data-min(area_data))/max(area_data);
    msize=40;
    scatter3(allen_ml_hva,allen_ap_hva,zeros(size(allen_ml_hva)),area_data*3,area_data_norm,'filled','MarkerEdgeColor',[.5 .5 .5]);
    axis equal off 
    box off
    view(2)
    set(gca, 'YDir','reverse')
    xlim([30 570])
    ylim([650 1050])
    hold on
    im=image(im_bound,'XData',0,'YData',0);
    im.AlphaData=max(im_bound,[],3);
    title(areas_hva{i})  
    colormap cool
    k=k+1;
end

figure;
k=1;
for i=1:15
    subplot(3,5,k)
    area_data=table2array(data_wb_gr(k+1,2:end));% exclude first injection
    area_data_norm=(area_data-min(area_data))/max(area_data);
    msize=40;
    scatter3(allen_ml_wb,allen_ap_wb,zeros(size(allen_ml_wb)), msize, area_data_norm,'filled','MarkerEdgeColor',[.5 .5 .5]);
    axis equal off 
    box off
    view(2)
    set(gca, 'YDir','reverse')
    xlim([30 570])
    ylim([650 1050])
    hold on
    im=image(im_bound,'XData',0,'YData',0);
    im.AlphaData=max(im_bound,[],3);
    title(areas_wb_gr{i})  
    colormap cool
    k=k+1;
end

figure;
k=1;
for i=1:15
    subplot(3,5,k)
    area_data=table2array(data_wb_gr(k+1,2:end));% exclude first injection
    area_data_norm=(area_data-min(area_data))/max(area_data);
    msize=40;
    scatter3(allen_ml_wb,allen_ap_wb,zeros(size(allen_ml_wb)),area_data*3, area_data_norm,'filled','MarkerEdgeColor',[.5 .5 .5]);
    axis equal off 
    box off
    view(2)
    set(gca, 'YDir','reverse')
    xlim([30 570])
    ylim([650 1050])
    hold on
    im=image(im_bound,'XData',0,'YData',0);
    im.AlphaData=max(im_bound,[],3);
    title(areas_wb_gr{i})  
    colormap cool 
    k=k+1;
end

%% HVA only

figure;
k=1;
for i=1:9
    subplot(3,3,k)
    area_data=table2array(data_hva(k+1,2:end));
    plot(ap(k+1:end),area_data(k+1:end),'o');
    h = lsline;
    set(h,'color','r')
    [r,p]=corr(ap(k+1:end)',area_data(k+1:end)');
    title({[areas{i} '  r=' num2str(round(r,2))],['p=' num2str(round(p,2))]})
    xlabel('anatomical AP (um)')
    ylabel('cell count (%)')
    k=k+1;
end

figure;
k=1;
for i=1:9
    subplot(3,3,k)
    area_data=table2array(data_hva(k+1,2:end));
    plot(ml(k+1:end),area_data(k+1:end),'o');
    h = lsline;
    set(h,'color','r')
    [r,p]=corr(ml(k+1:end)',area_data(k+1:end)');
    title({[areas{i} '  r=' num2str(round(r,2))],['p=' num2str(round(p,2))]})
    xlabel('anatomical ML (um)')
    ylabel('cell count (%)')
    k=k+1;
end

%% HVA only

figure;
k=1;
for i=1:9
    subplot(3,3,k)
    area_data=table2array(data_hva(k+1,2:end));
    plot(ap(k+1:end),area_data(k+1:end),'o');
    h = lsline;
    set(h,'color','r')
    [r,p]=corr(ap(k+1:end)',area_data(k+1:end)');
    title({[areas{i} '  r=' num2str(round(r,2))],['p=' num2str(round(p,2))]})
    xlabel('anatomical AP (um)')
    ylabel('cell count (%)')
    k=k+1;
end

figure;
k=1;
for i=1:9
    subplot(3,3,k)
    area_data=table2array(data_hva(k+1,2:end));
    plot(ml(k+1:end),area_data(k+1:end),'o');
    h = lsline;
    set(h,'color','r')
    [r,p]=corr(ml(k+1:end)',area_data(k+1:end)');
    title({[areas{i} '  r=' num2str(round(r,2))],['p=' num2str(round(p,2))]})
    xlabel('anatomical ML (um)')
    ylabel('cell count (%)')
    k=k+1;
end

%%
datafolder='';

load([datafolder filesep 'MLtable']);
areas=cellstr(MLtable.Area);
load([datafolder filesep 'MLtable']);
Inj=table2array(MLtable(:,2:end))';

Inj_new=Inj;
Inj_new(Inj_new==0)=NaN;
figure;
imagesc(blank);
hold on
scatter(ml,ap, 'or'); 
set(gca, 'YDir','reverse')
% plot boundary
% im_bound=imrotate(im_bound, 90);
im=image(im_bound,'XData',0,'YData',0);
im.AlphaData=max(im_bound,[],3);
axis equal
for i=1:length(ap)
ret=blank(ap(i),ml(i));
end
for i=1:length(ap)
ret(i)=blank(ap(i),ml(i));
end
ret=ret';
[~,order_ret]=sort(ret);
ret_new=ret(order_ret);
inj_new=inj(:,order_ret);

% plotting continuous values
figure;
k=1;
for i=1:9
    subplot(3,3,k)
    plot(ret_new,inj_new(i,:),'o');
    h = lsline;
    set(h,'color','r')
    r=corr(ret_new,inj_new(i,:)');
    title([areas{i} '  r=' num2str(round(r,2))])
    xlabel('azimuth retinotopy (deg)')
    ylabel('cell count (%)')
    k=k+1;
end
min(ret)
max(ret)
figure;
k=1;

for i=1:9
    subplot(3,3,k)
    plot(ret_new,inj_new(i,:),'o');
    h = lsline;
    set(h,'color','r')
    r=corr(ret_new,inj_new(i,:)');
    title([areas{i} '  r=' num2str(round(r,2))])
    xlabel('azimuth retinotopy (deg)')
    ylabel('cell count (%)')
    xlim([0 60])
    k=k+1;
end
