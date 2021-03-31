% -----------------------------------------------------------------------------
%      Plot ROIs on the 3D brain model colour-coded by injection site
% -----------------------------------------------------------------------------
% 2019-03 MM
% modified from cortex lab github allenCCF 'Analyze_Clicked_Points.m'

%% ENTER PARAMETERS AND FILE LOCATION
addpath(genpath('your github location\GitHub\allenCCF'));
addpath(genpath('your github location\GitHub\npy-matlab'));

% file location of object points
save_folder = ''; 

% directory of reference atlas files
annotation_volume_location = 'your github location\GitHub\allenCCF\annotation_volume_10um_by_index.npy';
structure_tree_location = 'your github location\Documents\GitHub\allenCCF\structure_tree_safe_2017.csv';

% roi tables to use
roi_save_name_suffix{1} = 'M19117_555.mat';
roi_save_name_suffix{2} = 'M19117_647.mat';
roi_save_name_suffix{3} = 'M19117_488.mat';

% Assign colors for each injection
InjSiteColors(1,:) = [0,1,1];% Cyan
InjSiteColors(2,:) = [1,1,0];% Yellow
InjSiteColors(3,:) = [1,0,1];% Magenta

% brain figure black or white
black_brain = true;

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

% load roi_table
for itable=1:length(roi_save_name_suffix)
    ROIs(itable) = load(fullfile(save_folder,roi_save_name_suffix{itable}));
end

% excluding VISp and white matter for plotting
clear VISp newROIs
for itable=1:length(roi_save_name_suffix)

    VISp_wm{itable}=ismember(ROIs(itable).roi.acronym,{'VISp','VISp1','VISp2/3','VISp4','VISp5','VISp6a', 'VISp6b', 'root', 'scwm', 'or', 'ccg', 'ccb', 'fp', 'cing', 'ec'});
    newROIs(itable).roi=ROIs(itable).roi(~VISp_wm{itable},:);
    
end

%% PLOTTING
% plot 3D brain
fwireframe = plotBrainGrid([], [], [], black_brain); hold on; 
fwireframe.InvertHardcopy = 'off';

for InjSite_num=1:size(newROIs,2)

    % get params
    ap= newROIs(InjSite_num).roi.AP_location;
    dv= newROIs(InjSite_num).roi.DV_location;
    ml= newROIs(InjSite_num).roi.ML_location;
    ann= newROIs(InjSite_num).roi.avIndex;

    % use AP, DV, and ML coordinates to get the point's position in the Allen atlas
    curr_objectPoints= NaN(length(ap),3);
    curr_objectPoints(:,1)= -ap/atlas_resolution + bregma(1);
    curr_objectPoints(:,2)= dv/atlas_resolution + bregma(2);
    curr_objectPoints(:,3)= -ml/atlas_resolution + bregma(3);
    
      
    figure(fwireframe); hold on               
    hp = plot3(curr_objectPoints(:,1), curr_objectPoints(:,3), curr_objectPoints(:,2), '.','linewidth',2, 'color',[InjSiteColors(InjSite_num,:) .2],'markers',4);hold on   
    
    clear ap dv ml ann
end
    
% savefig([save_folder filesep 'InjSite_colour_3d_animal.fig']) 
