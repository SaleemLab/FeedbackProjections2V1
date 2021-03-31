 % -----------------------------------------------------------------------------
%      Plot ROIs on the 3D brain model colour-coded by brain area
% -----------------------------------------------------------------------------
% 2019-02 MM
% modified from cortex lab github allenCCF 'Analyze_Clicked_Points.m'

%% ENTER PARAMETERS AND FILE LOCATION

% add sharp-track and npy-matlab
addpath(genpath('your github location\Documents\GitHub\allenCCF'))
addpath(genpath('your github location\Documents\GitHub\npy-matlab'))

% file location of object points
save_folder = '';

% directory of reference atlas files
annotation_volume_location = 'your github location\Documents\GitHub\allenCCF\annotation_volume_10um_by_index.npy';
structure_tree_location = 'your github location\Documents\GitHub\allenCCF\structure_tree_safe_2017.csv';

% get data file names
cd 'path to your datafiles'
filenames=dir('*Slide*');% whatever common term across of your files
roi_save_name_suffix=cell(size(filenames,1),1);
for i=1:size(filenames,1)
roi_save_name_suffix{i}=filenames(i).name;
end

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

% concatenate tables
roi=[];
for itable=1:length(roi_save_name_suffix)
   roiX= ROIs(itable).roi_table;
   roi=[roi ; roiX];
end

% get params
ap= roi.AP_location;
dv= roi.DV_location;
ml= roi.ML_location;
ann= roi.avIndex;

%% CREATE A CELL ARRAY CONTAINING AREA INDEX GROUPS TO BE PLOTTED IN SAME COLOUR

area{1,1}={'VISp'};
area{1,2}=st.index(ismember(st.acronym,{'VISp','VISp1','VISp2/3','VISp4','VISp5','VISp6a', 'VISp6b'}))';
area{2,1}={'VISpm'};
area{2,2}=st.index(ismember(st.acronym,{'VISpm','VISpm1','VISpm2/3','VISpm4','VISpm5','VISpm6a', 'VISpm6b'}))';
area{3,1}={'VISam'};
area{3,2}=st.index(ismember(st.acronym,{'VISam','VISam1','VISam2/3','VISam4','VISam5','VISam6a', 'VISam6b'}))';
area{4,1}={'VISa'};
area{4,2}=st.index(ismember(st.acronym,{'VISa','VISa1','VISa2/3','VISa4','VISa5','VISa6a', 'VISa6b'}))';
area{5,1}={'VISrl'};
area{5,2}=st.index(ismember(st.acronym,{'VISrl','VISrl1','VISrl2/3','VISrl4','VISrl5','VISrl6a', 'VISrl6b'}))';
area{6,1}={'VISal'};
area{6,2}=st.index(ismember(st.acronym,{'VISal','VISal1','VISal2/3','VISal4','VISal5','VISal6a', 'VISal6b'}))';
area{7,1}={'VISl'};
area{7,2}=st.index(ismember(st.acronym,{'VISl','VISl1','VISl2/3','VISl4','VISl5','VISl6a', 'VISl6b'}))';
area{8,1}={'VISli'};
area{8,2}=st.index(ismember(st.acronym,{'VISli','VISli1','VISli2/3','VISli4','VISli5','VISli6a', 'VISli6b'}))';
area{9,1}={'VISpor'};
area{9,2}=st.index(ismember(st.acronym,{'VISpor','VISpor1','VISpor2/3','VISpor4','VISpor5','VISpor6a', 'VISpor6b'}))';
area{10,1}={'VISpl'};
area{10,2}=st.index(ismember(st.acronym,{'VISpl','VISpl1','VISpl2/3','VISpl4','VISpl5','VISpl6a', 'VISpl6b'}))';
area{11,1}={'RSPagl'};
area{11,2}=st.index(ismember(st.acronym,{'RSPagl','RSPagl1','RSPagl2/3','RSPagl4','RSPagl5','RSPagl6a', 'RSPagl6b'}))';
area{12,1}={'RSPd'};
area{12,2}=st.index(ismember(st.acronym,{'RSPd','RSPd1','RSPd2/3','RSPd4','RSPd5','RSPd6a', 'RSPd6b'}))';
area{13,1}={'RSPv'};
area{13,2}=st.index(ismember(st.acronym,{'RSPv','RSPv1','RSPv2/3','RSPv4','RSPv5','RSPv6a', 'RSPv6b'}))';
area{14,1}={'TEa'};
area{14,2}=st.index(ismember(st.acronym,{'TEa','TEa1','TEa2/3','TEa4','TEa5','TEa6a', 'TEa6b'}))';
area{15,1}={'AUDpo'};
area{15,2}=st.index(ismember(st.acronym,{'AUDpo','AUDpo1','AUDpo2/3','AUDpo4','AUDpo5','AUDpo6a', 'AUDpo6b'}))';
area{16,1}={'ENTl'};
area{16,2}=st.index(ismember(st.acronym,{'ENTl','ENTl1','ENTl2','ENTl2/3','ENTl2a','ENTl2b', 'ENTl3','ENTl4','ENTl4/5','ENTl5', 'ENTl5/6','ENTl6a', 'ENTl6b'}))';
area{17,1}={'ENTm'};
area{17,2}=st.index(ismember(st.acronym,{'ENTm','ENTm1','ENTm2','ENTm2/3','ENTm2a','ENTm2b', 'ENTm3','ENTm4','ENTm4/5','ENTm5', 'ENTm5/6','ENTm6a', 'ENTm6b'}))';
area{18,1}={'PERI'};
area{18,2}=st.index(ismember(st.acronym,{'PERI','PERI1','PERI2/3','PERI5','PERI6a','PERI6b'}))';
area{19,1}={'ECT'};
area{19,2}=st.index(ismember(st.acronym,{'ECT','ECT1','ECT2/3','ECT5','ECT6a','ECT6b'}))';
area{20,1}={'ACA'};
area{20,2}=st.index(ismember(st.acronym,{'ACA','ACA1','ACA2/3','ACA5','ACA6a','ACA6b'}))';
area{21,1}={'ACAd'};
area{21,2}=st.index(ismember(st.acronym,{'ACAd','ACAd1','ACAd2/3','ACAd5','ACAd6a','ACAd6b'}))';
area{22,1}={'ACAv'};
area{22,2}=st.index(ismember(st.acronym,{'ACAv','ACAv1','ACAv2/3','ACAv5','ACAv6a','ACAv6b'}))';
area{23,1}={'MOs'};
area{23,2}=st.index(ismember(st.acronym,{'MOs','MOs1','MOs2/3','MOs5','MOs6a','MOs6b'}))';
area{24,1}={'CLA'};
area{24,2}=st.index(ismember(st.acronym,{'CLA'}))';
area{25,1}={'LA'};
area{25,2}=st.index(ismember(st.acronym,{'LA'}))';
area{26,1}={'LGd'};
area{26,2}=st.index(ismember(st.acronym,{'LGd','LGd-sh','LGd-co','LGd-ip'}))';
area{27,1}={'LP'};
area{27,2}=st.index(ismember(st.acronym,{'LP'}))';
area{28,1}={'LD'};
area{28,2}=st.index(ismember(st.acronym,{'LD'}))';

%---------INSERT MORE AREAS OF INTEREST HERE IF NEEDED-----------
AreaColors = [1 1 1; .7 0 0;  1 .6 0; .99 .68 .57; 0 .27 .1; .4 .6 .2; .7 .95 .3; 1 1 0.698;.15 .2 .58;.3 1 1; 1 .35 .65; .6 0 .7; .33 .15, .56; .7 .7 1; .02 .35 .55; .65 .4 .25; 1 .75 0; .4 .15 .02; .93 .44 .08; .42 .68 .84; .42 .68 .84; .42 .68 .84; .42 .32 .64; .35 .71 .67; .004 .4 .37; .87 .47 .68; .95 .71 .85; .99 .89 .94]; % add colors as needed
% AreaColors = [.7 0 0;  1 .6 0; .99 .68 .57; 0 .27 .1; .4 .6 .2; .7 .95 .3; 1 1 0.698;.15 .2 .58;.3 1 1; 1 .35 .65; .6 0 .7; .33 .15, .56; .7 .7 1; .02 .35 .55; .65 .4 .25; 1 .75 0; .4 .15 .02; .93 .44 .08; .42 .68 .84; .42 .68 .84; .42 .68 .84; .42 .32 .64; .35 .71 .67; .004 .4 .37; .87 .47 .68; .95 .71 .85; .99 .89 .94];
% order of colors: {'white','gold','orange','light pink','forest green','fern','green apple','light yellow','turquoise','navy blue','bubble gum','overcast sky','rawhide','purple','red', 'dark purple', 'dark blue' };

%% PLOTTING

% use AP, DV, and ML coordinates to get the point's position in the atlas
    curr_objectPoints= NaN(length(ap),3);
    curr_objectPoints(:,1)= -ap/atlas_resolution + bregma(1);
    curr_objectPoints(:,2)= dv/atlas_resolution + bregma(2);
    curr_objectPoints(:,3)= ml/atlas_resolution + bregma(3);

% plot 3D brain
fwireframe = plotBrainGrid([], [], [], black_brain); hold on; 
fwireframe.InvertHardcopy = 'off';
   
    for point = 1:size(curr_objectPoints,1)
    
       figure(fwireframe); hold on

           if ~isempty(find([area{:,2}]==ann{point}, 1))
               
               for iarea=1:size(area,1)
                   if  ~isempty(find(area{iarea,2}==ann{point}, 1))
                       area_num=iarea;
                   else
                       ;
                   end
               end
            hp = plot3(curr_objectPoints(point,1), curr_objectPoints(point,3), curr_objectPoints(point,2), '.','linewidth',2, 'color',[AreaColors(area_num,:) .2],'Markers',2);hold on
           else
            hp = plot3(curr_objectPoints(point,1), curr_objectPoints(point,3), curr_objectPoints(point,2), '.','linewidth',2, 'color',[0 0 0],'Markers',2);hold on
           end

    end

% plot dummies for legend
 for iarea=1:size(area,1)
    h(iarea)=plot([NaN,NaN],'.','MarkerSize',24,'color',AreaColors(iarea,:)); hold on
 end
    h(size(area,1)+1)=plot([NaN,NaN],'.','MarkerSize',24,'color',[.5 .5 .5]);

    legend(h,[[area{:,1}], {'other'}],'TextColor',[.7 .7 .7],'Box','off','color','none');
    
    title(['Animal ' num2str(animal)],'Color','w'); 
   
    
 savefig([save_folder filesep 'area_colour_3d_animal' num2str(animal) '.fig']) 
 
%% MOVIE MAKING

idx=1;

for i=1:36
   camorbit(-2.5,0,'data',[0 0 1])
   pause(0.05)
   drawnow
   frame = getframe(1);
   im{idx} = frame2im(frame);
   idx=idx+1;
end

for i=1:36
   camorbit(2.5,0,'data',[0 0 1])
   pause(0.05)
   drawnow
   frame = getframe(1);
   im{idx} = frame2im(frame);
   idx=idx+1;
end

for i=1:36
   camorbit(-2.5,0,'data',[1 0 0])
   pause(0.05)
   drawnow
   frame = getframe(1);
   im{idx} = frame2im(frame);
   idx=idx+1;
end

for i=1:36
   camorbit(2.5,0,'data',[1 0 0])
   pause(0.05)
   drawnow
   frame = getframe(1);
   im{idx} = frame2im(frame);
   idx=idx+1;
end

filename = 'brainAnimated_fast.gif'; % Specify the output file name
for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,[save_folder filesep filename],'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,[save_folder filesep filename],'gif','WriteMode','append','DelayTime',0.1);
    end
end

