% Load all ROI tables first
area_acr={'AUD','MO','SS','ACA','CLA','TEa','RSPagl','RSPd','RSPv','ECT','PRE','POST','PAR','ENTl','ENTm','VISpm'};%VISpm','VISam','VISrl','VISal','VISli','VISpor','VISpl'};%'VISa' 'VISl'
layers={'1','2/3','4','5','6'};

%% Filter data by area
%get idx of all cells in roi table
inj_used=[9,10,12,14,15,16,17,18,19,22,27,29,31,33,36,38];% only 6:end for whole brain

counts=nan(length(area_acr)*length(layers),length(inj_used)+1);    
L=array2table(counts,'VariableNames',{'lamina','inj9','inj10','inj12','inj14','inj15','inj16','inj17','inj18','inj19','inj22','inj27','inj29','inj31','inj33','inj36','inj38'});
L.lamina=cell(size(L,1),1);

%%
k=1;
for i=1:length(area_acr)
    % pick area
    area=area_acr{i};
    if strcmp(area,'POST')||strcmp(area,'PRE')||strcmp(area,'PAR') % these areas only have 3 layers
        layers={'1','2','3'};
    else
        layers={'1','2/3','4','5','6'};
    end
    
    for j=1:length(layers)
        %pick layer
        layer=layers{j};
        % label with area-layer
        L.lamina(k)={[area,layer]};
     
        % get idx of cells in this layer and area
        cells=0;
        for inj=1:length(inj_used)
            inj_num=inj+5; % first 5 injections not used for whole brain data
            % filter by area and layer 
            area_idx{inj}=contains(ROIs(inj_num).roi.acronym, area)& contains(ROIs(inj_num).roi.acronym, layer); 
            % cell count
            L{k,inj+1}= cells + nansum(area_idx{inj});
        end

        % checking the cells that got filtered
        A=ROIs(5+5).roi.acronym(area_idx{5});
        
        k=k+1;
    end
end
writetable(L,[filename '.xlsx'])

    
    
    

