% Some scripts for:
% - visualization of all labelled cells per are in 3D volume
% - visualisation of injection sites on retinotopic maps

%% Get inj coords
% get the azimuth/altitude coords for each injection
% load inj coords
load('your path\Inj_coords.mat')
% injections used
inj_used=[3,5,6,7,8,9,10,12,14,15,16,17,18,19,22,27,29,31,33,36,38];
coords_used=Inj_coords(inj_used,:);

% transform the coords to allen CCF coords
addpath(genpath('your path\GitHub\allenCCF'))
bregma=allenCCFbregma();
atlas_resolution=0.01;

allen_ml=-abs(coords_used(:,2))/atlas_resolution+bregma(3);
allen_ap=-coords_used(:,1)/atlas_resolution+bregma(1);

% get retinotopy coords
load('your path\allen_borders.mat')
azi_im=load('your path\meanAzimuthMap.mat');
alt_im=load('your path\meanAltitudeMap.mat');

% azimuth map
allen_map_azi=nan(1320,1140);
azi_im.im(azi_im.im==1000) = nan;
allen_map_azi(501:1100,1:600)=azi_im.im;
figure;
imagesc(allen_map_azi)
hold on
scatter(allen_ml, allen_ap, 'k.')
axis equal
box off
title('azimuth map')

% altitude map
allen_map_alt=nan(1320,1140);
alt_im.im(alt_im.im==1000) = nan;
allen_map_alt(501:1100,1:600)=alt_im.im;
figure;
imagesc(allen_map_alt)
hold on
scatter(allen_ml, allen_ap, 'k.')
axis equal
box off
title('altitude map')

% get altitude retinotopy
for i=1:length(allen_ap)
azi_ret(i)=allen_map_azi(allen_ap(i), allen_ml(i));
end
for i=1:length(allen_ap)
alt_ret(i)=allen_map_alt(allen_ap(i), allen_ml(i));
end

%% Load cell coordinates
% cell coordinates stored in roi tables
save_folder='your path\roi_tables';

cd(save_folder)
filenames=dir;
roi_save_name_suffix=cell(size(filenames,1),1);

for i=1:size(filenames,1)-2
    roi_save_name_suffix{i}=filenames(i+2).name;
end

% load roi_table
for itable=1:length(roi_save_name_suffix)-2
    ROIs(itable) = load(fullfile(save_folder,roi_save_name_suffix{itable}));
end 

%% Filter data by area

% choose an area to look at
area='AUD'; % 'Lateral visual area','VISpm','RSP', etc. 

% get idx of all cells in roi table
% long area names are 'name', shorthand area names are 'acronym'

% % for name
% for inj=1:length(ROIs)
%     area_idx{inj}=contains(ROIs(inj).roi.name, area);
% end

% for acronym
for inj=1:length(ROIs)
    area_idx{inj}=contains(ROIs(inj).roi.acronym, area);
end

% checking the areas that got filtered
A=ROIs(10).roi.acronym(area_idx{10});

%% Assemble cells from all injections
AP=[];ML=[];DV=[];InjNum=[];

for inj=1:length(ROIs)
    ap=ROIs(inj).roi.AP_location(area_idx{inj});
    dv=ROIs(inj).roi.DV_location(area_idx{inj});
    ml=ROIs(inj).roi.ML_location(area_idx{inj});
    inj_num=inj*ones(sum(area_idx{inj}),1); % injection num (not the unique injection ID)
    % concatenate
    AP=[AP;ap];
    ML=[ML;ml];
    DV=[DV;dv];
    InjNum=[InjNum;inj_num];
    tot_area(inj)=length(ap);
end
CellID=1:length(ML);
C=table(CellID',InjNum,ML,AP,DV,'VariableName',{'CellID','InjNum','ML','AP','DV'});

%% Plot cells to V1 azimuth correlation
% make colormap
cm=colormap(cool(60));
%get color idx
for i=1:length(azi_ret)
    azi_ret_col(i,:)=cm(round(azi_ret(i)),:);
end
    
% plot cells
figure;
for inj=1:length(ROIs)%inj=[9,10,12]
    ap=ROIs(inj).roi.AP_location(area_idx{inj});
    dv=ROIs(inj).roi.DV_location(area_idx{inj});
    ml=ROIs(inj).roi.ML_location(area_idx{inj}); 
    h=scatter3(abs(ml), ap, dv, 20,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',azi_ret_col(inj,:));
    alpha = 0.3;
    set(h, 'MarkerFaceAlpha', alpha)%'MarkerEdgeAlpha', alpha,
    hold on
end
xlabel('M-L');ylabel('A-P');zlabel('D-V')
axis(gca,'equal','tight')
set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.01 0.02], 'color', 'none', 'box','off','Zdir','reverse')
colormap(cm)
% colorbar

title([ area ', azimuth color code'])

%% Different views (Azimuth)
figure;
subplot(1,3,1)
    for  inj=1:length(ROIs)%inj=[9,10,12]
        ap=ROIs(inj).roi.AP_location(area_idx{inj});
        dv=ROIs(inj).roi.DV_location(area_idx{inj});
        ml=ROIs(inj).roi.ML_location(area_idx{inj});
        h=scatter3(abs(ml), ap, dv, 20,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',azi_ret_col(inj,:));
        alpha = 0.5;
        set(h, 'MarkerFaceAlpha', alpha)%'MarkerEdgeAlpha', alpha,
        hold on
    end
    xlabel('M-L (mm)');ylabel('A-P (mm)');zlabel('D-V (mm)')
    colormap(cm)
    view(0,0)
    title({ area, 'azimuth color code, coronal view'})
    axis(gca,'equal','tight') 
    set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.01 0.02], 'color', 'none', 'box','off','Zdir','reverse')

    subplot(1,3,2)
    for  inj=1:length(ROIs)%inj=[9,10,12]
        ap=ROIs(inj).roi.AP_location(area_idx{inj});
        dv=ROIs(inj).roi.DV_location(area_idx{inj});
        ml=ROIs(inj).roi.ML_location(area_idx{inj});
        h=scatter3(abs(ml), ap, dv, 20,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',azi_ret_col(inj,:));
        alpha = 0.5;
        set(h, 'MarkerFaceAlpha', alpha)%'MarkerEdgeAlpha', alpha,
        hold on
    end
    xlabel('M-L (mm)');ylabel('A-P (mm)');zlabel('D-V (mm)')
    colormap(cm)
    view(90,0)
    title({ area, 'azimuth color code, saggital view'})
    axis(gca,'equal','tight')
    set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.01 0.02], 'color', 'none', 'box','off','Zdir','reverse')

subplot(1,3,3)
    for  inj=1:length(ROIs)%inj=[9,10,12]
        ap=ROIs(inj).roi.AP_location(area_idx{inj});
        dv=ROIs(inj).roi.DV_location(area_idx{inj});
        ml=ROIs(inj).roi.ML_location(area_idx{inj});
        h=scatter3(abs(ml), ap, dv, 20,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',azi_ret_col(inj,:));
        alpha = 0.5;
        set(h, 'MarkerFaceAlpha', alpha)%'MarkerEdgeAlpha', alpha,
        hold on
    end
    xlabel('M-L (mm)');ylabel('A-P (mm)');zlabel('D-V (mm)')
    colormap(cm)
    view(0,90)
    title({ area, 'azimuth color code, top-down view'})
    axis(gca,'equal','tight')
    set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.01 0.02], 'color', 'none', 'box','off')

%% Plot cells to V1 altitude correlation
% make colormap
cm=colormap(cool(25));% altitude lims [-10 15]

%get color idx
for i=1:length(alt_ret)
    alt_ret_col(i,:)=cm(round(alt_ret(i))+10,:);
end

figure;
for inj=1:length(ROIs)
    ap=ROIs(inj).roi.AP_location(area_idx{inj});
    dv=ROIs(inj).roi.DV_location(area_idx{inj});
    ml=ROIs(inj).roi.ML_location(area_idx{inj});
    scatter3(abs(ml), ap, dv, 20,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',alt_ret_col(inj,:))
    hold on
end
xlabel('M-L');ylabel('A-P');zlabel('D-V')
axis(gca,'equal','tight')
set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.01 0.02], 'color', 'none', 'box','off','Zdir','reverse')
colormap(cm)
% colorbar

title([ area ' altitude color code'])

%% Different views (Altitude)
figure;
subplot(1,3,1)
    for inj=1:length(ROIs)
        ap=ROIs(inj).roi.AP_location(area_idx{inj});
        dv=ROIs(inj).roi.DV_location(area_idx{inj});
        ml=ROIs(inj).roi.ML_location(area_idx{inj});
        h=scatter3(abs(ml), ap, dv, 20,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',alt_ret_col(inj,:));
        alpha = 0.5;
        set(h, 'MarkerFaceAlpha', alpha)
        hold on
    end
    xlabel('M-L');ylabel('A-P');zlabel('D-V')
    colormap(cm)
    view(0,0)
    title({ area, 'altitude color code, coronal view'})
    axis(gca,'equal','tight') 
    set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.01 0.02], 'color', 'none', 'box','off','Zdir','reverse')

    subplot(1,3,2)
    for inj=1:length(ROIs)
        ap=ROIs(inj).roi.AP_location(area_idx{inj});
        dv=ROIs(inj).roi.DV_location(area_idx{inj});
        ml=ROIs(inj).roi.ML_location(area_idx{inj});
        h=scatter3(abs(ml), ap, dv, 20,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',alt_ret_col(inj,:));
        alpha = 0.5;
        set(h, 'MarkerFaceAlpha', alpha)
        hold on
    end
    xlabel('M-L');ylabel('A-P');zlabel('D-V')
    colormap(cm)
    view(90,0)
    title({ area, 'altitude color code, saggital view'})
    axis(gca,'equal','tight')
    set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.01 0.02], 'color', 'none', 'box','off','Zdir','reverse')

subplot(1,3,3)
    for inj=1:length(ROIs)
        ap=ROIs(inj).roi.AP_location(area_idx{inj});
        dv=ROIs(inj).roi.DV_location(area_idx{inj});
        ml=ROIs(inj).roi.ML_location(area_idx{inj});
        h=scatter3(abs(ml), ap, dv, 20,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',alt_ret_col(inj,:));
        alpha = 0.5;
        set(h, 'MarkerFaceAlpha', alpha)
        hold on
    end
    xlabel('M-L');ylabel('A-P');zlabel('D-V')
    colormap(cm)
    view(0,90)
    title({ area, 'altitude color code, top-down view'})
    axis(gca,'equal','tight')
    set(gca, 'fontsize', 10, 'TickDir', 'out', 'TickLength', [0.01 0.02], 'color', 'none', 'box','off')

%% Plot inj color codes
figure;
subplot(1,2,1)
    msize=40;
    h=scatter3(allen_ml,allen_ap,zeros(size(allen_ml)),msize,alt_ret_col,'filled','MarkerEdgeColor',[.5 .5 .5]);
    axis equal off 
    box off
    view(2)
    set(gca, 'YDir','reverse')
    alpha = 0.5;% transparent markers
    set(h, 'MarkerFaceAlpha', alpha)
    xlim([30 570])
    ylim([650 1050])
    hold on
    im=image(im_bound,'XData',0,'YData',0);
    im.AlphaData=max(im_bound,[],3);
    title('V1 altitude color code')  
% add text label
for inj=1:21
    text(allen_ml(inj)+1, allen_ap(inj)+1, num2str(inj), 'Fontsize', 10);
end

subplot(1,2,2)
    msize=40;
    h=scatter3(allen_ml,allen_ap,zeros(size(allen_ml)),msize,azi_ret_col,'filled','MarkerEdgeColor',[.5 .5 .5]);
    axis equal off 
    box off
    view(2)
    set(gca, 'YDir','reverse')
    alpha = 0.5;% transparent markers
    set(h, 'MarkerFaceAlpha', alpha)
    xlim([30 570])
    ylim([650 1050])
    hold on
    im=image(im_bound,'XData',0,'YData',0);
    im.AlphaData=max(im_bound,[],3);
    title('V1 azimuth color code')  
% add text label
for inj=1:21
    text(allen_ml(inj)+1, allen_ap(inj)+1, num2str(inj), 'Fontsize', 10);
end

%% Plot RedGrayBlue map (azi)
figure;
    imagesc(allen_map_azi,[0 60])
    contourcmap('RedGrayBlue',[0:5:60]);
    colorbar
    hold on
    msize=20;
    h=scatter3(allen_ml,allen_ap,zeros(size(allen_ml)),msize,'k','filled');
    axis equal off 
    box off
    view(2)
    set(gca, 'YDir','reverse')
    % transparent markers
%     alpha = 0.5;
%     set(h, 'MarkerFaceAlpha', alpha)
    xlim([30 570])
    ylim([650 1050])
    hold on
    im=image(im_bound,'XData',0,'YData',0);
    im.AlphaData=max(im_bound,[],3);
    title('V1 altitude color code')  

%% Plot RedGrayBlue map (alt)
figure;
    imagesc(allen_map_alt,[-10 15])
    contourcmap('RedGrayBlue',[-10:2.5:15]);
    colorbar
    hold on
    msize=20;
    h=scatter3(allen_ml,allen_ap,zeros(size(allen_ml)),msize,'k','filled');
    axis equal off 
    box off
    view(2)
    set(gca, 'YDir','reverse')
    % transparent markers
%     alpha = 0.5;
%     set(h, 'MarkerFaceAlpha', alpha)
    xlim([30 570])
    ylim([650 1050])
    hold on
    im=image(im_bound,'XData',0,'YData',0);
    im.AlphaData=max(im_bound,[],3);
    title('V1 altitude color code')           

%% Data to include for 3D voxels
% for HVA
x=abs(C.ML);
y=C.AP;
z=C.DV;
cell_ID=C.CellID;

%% Data to include for 3D voxels
% excluding injection closest to PM   
inj_pm=C.InjNum~=9; % 9th inj
cell_ID=C.CellID(inj_pm);
x=abs(C.ML(inj_pm));
y=C.AP(inj_pm);
z=C.DV(inj_pm);

%% Data to include for 3D voxels
% for whole brain    
inj_wb=C.InjNum>5;
cell_ID=C.CellID(inj_wb);
x=abs(C.ML(inj_wb));
y=C.AP(inj_wb);
z=C.DV(inj_wb);

%% Data to include for 3D voxels
% for 3 injections    
inj_3=C.InjNum==9|C.InjNum==10|C.InjNum==12;
cell_ID=C.CellID(inj_3);
x=abs(C.ML(inj_3));
y=C.AP(inj_3);
z=C.DV(inj_3);

%% Binning to 3D voxels
binSize=0.1;%mm
xBins= [floor(min(x)*10)/10: binSize: ceil(max(x)*10)/10];
yBins= [floor(min(y)*10)/10: binSize: ceil(max(y)*10)/10];
zBins= [floor(min(z)*10)/10: binSize: ceil(max(z)*10)/10];
D= zeros(length(xBins),length(yBins),length(zBins));

for icell=1:length(x)
    xi= find((x(icell)>xBins),1,'last');
    yi= find((y(icell)>yBins),1,'last');
    zi= find((z(icell)>zBins),1,'last');
    norm_cell=1/tot_area(C.InjNum(cell_ID(icell))); % per cell weight, normalized by total cells in area per inj    
    norm_ret=(azi_ret(C.InjNum(cell_ID(icell)))-30)/30; %retinotopy normalized to -1 to 1
    D(xi,yi,zi)= D(xi,yi,zi)+norm_cell*norm_ret; % norm cell count * norm retinotopy
end

%% Binning to 3D voxels using nans (azi)
binSize=0.1;%mm
xBins= [floor(min(x)*10)/10: binSize: ceil(max(x)*10)/10];
yBins= [floor(min(y)*10)/10: binSize: ceil(max(y)*10)/10];
zBins= [floor(min(z)*10)/10: binSize: ceil(max(z)*10)/10];
D= nan(length(xBins),length(yBins),length(zBins));

for icell=1:length(x)
    xi= find((x(icell)>xBins),1,'last');
    yi= find((y(icell)>yBins),1,'last');
    zi= find((z(icell)>zBins),1,'last');
    norm_cell=1/tot_area(C.InjNum(cell_ID(icell))); % per cell weight, normalized by total cells in area per inj    
    norm_ret=(azi_ret(C.InjNum(cell_ID(icell)))-30)/30; %retinotopy normalized to -1 to 1 (range 0 to 60deg)
    D(xi,yi,zi)= nansum([D(xi,yi,zi), norm_cell*norm_ret]); % norm cell count * norm retinotopy
end

%% Binning to 3D voxels using nans (alt)
binSize=0.1;%mm
xBins= [floor(min(x)*10)/10: binSize: ceil(max(x)*10)/10];
yBins= [floor(min(y)*10)/10: binSize: ceil(max(y)*10)/10];
zBins= [floor(min(z)*10)/10: binSize: ceil(max(z)*10)/10];
D= nan(length(xBins),length(yBins),length(zBins));

for icell=1:length(x)
    xi= find((x(icell)>xBins),1,'last');
    yi= find((y(icell)>yBins),1,'last'); 
    zi= find((z(icell)>zBins),1,'last');
    norm_cell=1/tot_area(C.InjNum(cell_ID(icell))); % per cell weight, normalized by total cells in area per inj    
    norm_ret=((alt_ret(C.InjNum(cell_ID(icell)))+10)-12.5)/12.5; %retinotopy normalized to -1 to 1 (range -10 to 15deg)
    D(xi,yi,zi)= nansum([D(xi,yi,zi), norm_cell*norm_ret]); % norm cell count * norm retinotopy
end

%% Binning to 3D voxels using nans per inj distribution

figure;
for inj_num=1:21
    
inj_plot=C.InjNum==inj_num;
cell_ID=C.CellID(inj_plot);
x=abs(C.ML(inj_plot));
y=C.AP(inj_plot);
z=C.DV(inj_plot);

binSize=0.1;%mm
xBins= [floor(min(x)*10)/10: binSize: ceil(max(x)*10)/10];
yBins= [floor(min(y)*10)/10: binSize: ceil(max(y)*10)/10];
zBins= [floor(min(z)*10)/10: binSize: ceil(max(z)*10)/10];
D= nan(length(xBins),length(yBins),length(zBins));

for icell=1:length(x)
    xi= find((x(icell)>xBins),1,'last');
    yi= find((y(icell)>yBins),1,'last');
    zi= find((z(icell)>zBins),1,'last');
    D(xi,yi,zi)= nansum([D(xi,yi,zi), 1]); % norm cell count * norm retinotopy
end

D_cor=permute(D,[3,1,2]);%for coronal viewing

clear cm
cm(:,1)=[linspace(0.9,azi_ret_col(inj_num,1),60),azi_ret_col(inj_num,1)*ones(1,41)];
cm(:,2)=[linspace(0.9,azi_ret_col(inj_num,2),60),azi_ret_col(inj_num,2)*ones(1,41)];
cm(:,3)=[linspace(0.9,azi_ret_col(inj_num,3),60),azi_ret_col(inj_num,3)*ones(1,41)];

[~,Ind] = sort(azi_ret);
sp=find(Ind==inj_num);
ax(inj_num)=subplot(5,5,sp);    
    mean_im=(nanmean(D_cor,3));
    h=imagesc(mean_im);%[-0.02 0.02][1 6]
    axis equal tight
    colormap(ax(inj_num),cm)
    set(h,'alphadata',~isnan(mean_im))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
    title([area ' mean coronal inj' num2str(inj_num)])
    colorbar
    hold on
end

%% Plotting coronal view

cm=RedGrayBlue(100,0.8);
% view coronal slices
D_cor=permute(D,[3,1,2]);%for coronal viewing

figure(106);
for i=1:length(yBins)
    subplot(4,ceil(length(yBins)/4),i)
    h=imagesc(D_cor(:,:,i),[-0.03 0.03]);
    axis equal tight
    title(['slice' num2str(i)])
    colormap(cm) 
    set(h,'alphadata',~isnan(D_cor(:,:,i)))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])

end

% Create textbox
annotation(figure(106),'textbox',...
    [0.42 0.96 0.19 0.04],...
    'String',{[area ' azimuth retinotopy coronal sections']},...
    'LineStyle','none',...
    'FontSize',14);

figure;
    mean_im=(nanmean(D_cor,3));
    h=imagesc(mean_im,[-0.02 0.02]);
    axis equal tight
    colormap(cm) 
    set(h,'alphadata',~isnan(mean_im))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
    title([area ' mean coronal'])
    colorbar
    
% Plotting top-down view

% view top-down slices
D_top=permute(D,[2,1,3]);%for top-down viewing

figure(108);
for i=1:length(zBins)
    subplot(4,ceil(length(zBins)/4),i)
    h=imagesc(flipud(D_top(:,:,i)),[-0.03 0.03]);
    axis equal tight
    title(['slice' num2str(i)])
    colormap(cm) 
    set(h,'alphadata',~isnan(flipud(D_top(:,:,i))))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
end

% Create textbox
annotation(figure(108),'textbox',...
    [0.42 0.96 0.19 0.04],...
    'String',{[area ' azimuth retinotopy top-down sections']},...
    'LineStyle','none',...
    'FontSize',14);

figure;
    mean_im=flipud(nanmean(D_top,3));
    h=imagesc(mean_im,[-0.02 0.02]);
    axis equal tight
    title(['slice' num2str(i)])
    colormap(cm) 
    set(h,'alphadata',~isnan(mean_im))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
    title([area ' mean top-down'])
    colorbar
    
% Plotting saggital view

% view top-down slices
D_sag=permute(D,[3,2,1]);%for saggital viewing

figure(109);
for i=1:length(xBins)
    subplot(4,ceil(length(xBins)/4),i)
    h=imagesc(D_sag(:,:,i),[-0.03 0.03]);
    axis equal tight
    title(['slice' num2str(i)])
    colormap(cm) 
    set(h,'alphadata',~isnan(D_sag(:,:,i)))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
end

% Create textbox
annotation(figure(109),'textbox',...
    [0.42 0.96 0.19 0.04],...
    'String',{[area ' azimuth retinotopy sag sections']},...
    'LineStyle','none',...
    'FontSize',14);

figure;
    mean_im=nanmean(D_sag,3);
    h=imagesc(mean_im,[-0.02 0.02]);
    axis equal tight
    title(['slice' num2str(i)])
    colormap(cm) 
    set(h,'alphadata',~isnan(mean_im))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
    title([area ' mean saggital'])
    colorbar
    
%% Plot mean images (alt)
lims=[-0.02 0.02];

h=figure;
subplot(1,3,1)
    mean_im=(nanmean(D_cor,3));
    h=imagesc(mean_im,lims);
    axis equal tight
    colormap(cm) 
    set(h,'alphadata',~isnan(mean_im))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
    xlabel('M-L (mm)');ylabel('D-V (mm)');
    title([area ' mean azimuth coronal'])
    box off
    colorbar
    
subplot(1,3,2)
    mean_im=nanmean(D_sag,3);
    h=imagesc(mean_im,lims);
    axis equal tight 
    title(['slice' num2str(i)])
    colormap(cm) 
    set(h,'alphadata',~isnan(mean_im))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
    xlabel('A-P (mm)');ylabel('D-V (mm)');
    title([area ' mean azimuth saggital'])
    box off
    colorbar
   
subplot(1,3,3)
    mean_im=flipud(nanmean(D_top,3));
    h=imagesc(mean_im,lims);
    axis equal tight
    title(['slice' num2str(i)])
    colormap(cm) 
    set(h,'alphadata',~isnan(mean_im))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
    xlabel('M-L (mm)');ylabel('A-P (mm)');
    title([area ' mean azimuth top-down'])
    box off
    colorbar

set(gcf, 'Position', get(0, 'Screensize'));
    
%% Plot mean images (azi)
lims=[-0.04 0.04];

figure;
subplot(1,3,1)
    mean_im=(nanmean(D_cor,3));
    h=imagesc(mean_im,lims);
    axis equal tight
    colormap(cm) 
    set(h,'alphadata',~isnan(mean_im))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
    xlabel('M-L (mm)');ylabel('D-V (mm)');
    title([area ' mean azimuth coronal'])
    box off
    colorbar
    
subplot(1,3,2)
    mean_im=nanmean(D_sag,3);
    h=imagesc(mean_im,lims);
    axis equal tight 
    title(['slice' num2str(i)])
    colormap(cm) 
    set(h,'alphadata',~isnan(mean_im))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
    xlabel('A-P (mm)');ylabel('D-V (mm)');
    title([area ' mean azimuth saggital'])
    box off
    colorbar
   
subplot(1,3,3)
    mean_im=flipud(nanmean(D_top,3));
    h=imagesc(mean_im,lims);
    axis equal tight
    title(['slice' num2str(i)])
    colormap(cm) 
    set(h,'alphadata',~isnan(mean_im))% Set nan values to transparent: 
    set(gca, 'color', 'white')
    set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
    xlabel('M-L (mm)');ylabel('A-P (mm)');
    title([area ' mean azimuth top-down'])
    box off
    colorbar
    