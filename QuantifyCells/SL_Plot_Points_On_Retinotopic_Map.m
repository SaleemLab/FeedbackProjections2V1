%% ENTER PARAMETERS AND FILE LOCATION

% file location of object points
save_folder = '';

% directory of reference atlas files
annotation_volume_location = 'your github location\allenCCF-master\annotation_volume_10um_by_index.npy';
structure_tree_location = 'your github location\allenCCF-master\structure_tree_safe_2017.csv';

% name of the saved roi tables
roi_save_name_suffix{1} = 'M19114';
roi_save_name_suffix{2} = 'M19115';
roi_save_name_suffix{3} = 'M19116';
roi_save_name_suffix{4} = 'M19117';
roi_save_name_suffix{5} ='M19118';
roi_save_name_suffix{6} ='M19119';
roi_save_name_suffix{7} ='M19121';
roi_save_name_suffix{8} ='M19122';
roi_save_name_suffix{9} ='M19123';

% Assign colors for each injection
InjSiteColors(1,:) = [0,0,0];%Black
InjSiteColors(2,:) = [1,0.6,0];%Orange
InjSiteColors(3,:) = [0.6,0,0.7];%Purple
InjSiteColors(4,:) = [1,0,0];%Red
InjSiteColors(5,:) = [0.87,0.47,0.68];% Pink
InjSiteColors(6,:) = [0.99,0.68,0.57];% Light red
InjSiteColors(7,:) = [0.4,0.6,0.2];%Green
InjSiteColors(8,:) = [0.15,0.2,0.58];%Blue
InjSiteColors(9,:) = [0.3,1,1];%Turquoise

% brain figure black or white
 black_brain = false;

% generate needed values
bregma = allenCCFbregma(); % bregma position in reference data space
atlas_resolution = 0.010; % mm

%% LOAD THE REFERENCE ANNOTATIONS AND ROI TABLE

% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end

% load roi_tables
for itable=1:length(roi_save_name_suffix)
    ROIs(itable) = load(fullfile(save_folder,roi_save_name_suffix{itable}));
end

%% PLOTTING

black_brain = false;
 
fwireframe = plotBrainGrid([], [], [], black_brain); hold on; 
fwireframe.InvertHardcopy = 'off';
    
figure(fwireframe); hold on
for InjSite_num=1:size(ROIs,2)

    % get params
    ap= ROIs(InjSite_num).roi_table.AP_location;
    ml= ROIs(InjSite_num).roi_table.ML_location;
    dv= ROIs(InjSite_num).roi_table.DV_location;
    ann= ROIs(InjSite_num).roi_table.avIndex;

    % use AP, DV, and ML coordinates to get the point's position in the atlas
    curr_objectPoints= NaN(length(ap),3);
    curr_objectPoints(:,1)= -ap/atlas_resolution + bregma(1);
%     curr_objectPoints(:,2)= dv/atlas_resolution + bregma(2);
    curr_objectPoints(:,3)= ml/atlas_resolution + bregma(3);
    
    ROIs(InjSite_num).roi_table.AP_3d_coords=curr_objectPoints(:,1);
%     ROIs(InjSite_num).roi_table.DV_3d_coords=curr_objectPoints(:,2);
    ROIs(InjSite_num).roi_table.ML_3d_coords=curr_objectPoints(:,3);
    
    hp = scatter(curr_objectPoints(:,3), -curr_objectPoints(:,1),5,[InjSiteColors(InjSite_num,:)],'filled');hold on   
% hp = plot3(curr_objectPoints(:,1), curr_objectPoints(:,3), curr_objectPoints(:,2), '.','linewidth',2, 'color',[InjSiteColors(InjSite_num,:) .2],'markers',10);hold on   
     xlim([0 1140]);
     ylim([-1320 0]);
%      zlim([0 800]);
      axis on
      pbaspect([1 1.16 1]);
end

% 
% axis equal 
% axis tight

% plot dummies for legend
 for iarea=1:size(ROIs,2)
    h(iarea)=plot([NaN,NaN],'.','MarkerSize',20,'color',InjSiteColors(iarea,:)); hold on
 end
 
 legend(h,roi_save_name_suffix{1,:},'TextColor',[.0 .0 .0],'Box','off','color','none','Interpreter', 'none');
    

%% Overlay azimuth map and CCF boundary

% Azimuth/Elevation maps
    Azimuth_map='your github location\mean_azimuth_map.npy';
    Azimuth=readNPY(Azimuth_map);
    
% Load Allen borders (tif)
im_dir='your github location\Allen borders';
filename_tif=['CCF_borders_10um' '.tif'];

    InfoImage.tif=imfinfo([im_dir filesep filename_tif]);
    im_tif=zeros(InfoImage.tif.Height,InfoImage.tif.Width,length(InfoImage.tif),'double');

    TifLink = Tiff([im_dir filesep filename_tif], 'r');
    for i=1:length(InfoImage.tif)
       TifLink.setDirectory(i);
       im_tif(:,:,i)=TifLink.read();
    end
    TifLink.close();

% Remove outer area
    Az=Azimuth;
    Az(Az==1000) = nan;
    caxis([-30 90]);
    % Al=Altitude;
    % Al(Al==1000) = nan;
    % caxis([-30 30]);

% Overlay azimuth map
    Az_map=flipdim(Az,1);
    hold on
    Az_m=image(Az_map,'XData',0,'YData',-1100,'AlphaData',0.9);
    hold off

% Overlay CCF boundaries 
    hold on
    im=image(im_tif,'XData',0,'YData',0);
    im.AlphaData=max(im_tif,[],3);
    hold off
    
% Remove outer area
    Az=Azimuth;
    Az(Az==1000) = nan;
    caxis([-30 90]);

% Creat blank matrix 1320x1140
    blank=zeros(size(im_tif));
    blank(501:1100,1:600)=Az;
    contourcmap('jet',[-30:5:90]); % Az
    
    
    
    