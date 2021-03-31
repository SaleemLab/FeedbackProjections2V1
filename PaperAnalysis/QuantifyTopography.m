function T=QuantifyTopography(T, ROIs, azi_ret, alt_ret, thresh, ret)

if nargin<5
    thresh=0.004;
end

if nargin<6
    ret={'azimuth','altitude'};
end

for r=1:length(ret)
    retinotopy=ret{r};
    
for iarea=1:size(T,1)
    area=T.area{iarea}

%% Filter data by area
area_type=T.area_type{iarea};

switch area_type
    
case 'acr'
    % for acronym
    for inj=1:length(ROIs)
        area_idx{inj}=contains(ROIs(inj).roi.acronym, area);
    end
    
case 'name'
    % for name
    for inj=1:length(ROIs)
        area_idx{inj}=contains(ROIs(inj).roi.name, area);
    end
    
end

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

%% Data to include for 3D voxels
area_loc=T.area_loc{iarea};

switch area_loc
    
case 'hva'
    x=abs(C.ML);
    y=C.AP;
    z=C.DV;
    cell_ID=C.CellID;
    
case 'wb'
    inj_wb=C.InjNum>5;
    x=abs(C.ML(inj_wb));
    y=C.AP(inj_wb);
    z=C.DV(inj_wb);
    cell_ID=C.CellID(inj_wb);
    
end

%% Binning to 3D voxels using nans
clear D D_dist
binSize=0.1;%mm
xBins= [floor(min(x)*10)/10: binSize: ceil(max(x)*10)/10];
yBins= [floor(min(y)*10)/10: binSize: ceil(max(y)*10)/10];
zBins= [floor(min(z)*10)/10: binSize: ceil(max(z)*10)/10];
D= nan(length(xBins),length(yBins),length(zBins));

switch retinotopy
    case 'azimuth'
    for icell=1:length(x)
        xi= find((x(icell)>xBins),1,'last');
        yi= find((y(icell)>yBins),1,'last');
        zi= find((z(icell)>zBins),1,'last');
        if ~isempty(xi) && ~isempty(yi) && ~isempty(zi)
            norm_cell=1/tot_area(C.InjNum(cell_ID(icell))); % per cell weight, normalized by total cells in area per inj    
            norm_ret=(azi_ret(C.InjNum(cell_ID(icell)))-30)/30; %retinotopy normalized to -1 to 1 (range 0 to 60deg)
            D(xi,yi,zi)= nansum([D(xi,yi,zi), norm_cell*norm_ret]); % norm cell count * norm retinotopy
        else
            ;
        end
    end
    
    case 'altitude'
    for icell=1:length(x)
        xi= find((x(icell)>xBins),1,'last');
        yi= find((y(icell)>yBins),1,'last'); 
        zi= find((z(icell)>zBins),1,'last');
        if ~isempty(xi) && ~isempty(yi) && ~isempty(zi)
            norm_cell=1/tot_area(C.InjNum(cell_ID(icell))); % per cell weight, normalized by total cells in area per inj    
            norm_ret=((alt_ret(C.InjNum(cell_ID(icell)))+10)-12.5)/12.5; %retinotopy normalized to -1 to 1 (range -10 to 15deg)
            D(xi,yi,zi)= nansum([D(xi,yi,zi), norm_cell*norm_ret]); % norm cell count * norm retinotopy
        else
            ;
        end
    end

end

%% Collapse to 2D top-down 
D_top=permute(D,[2,1,3]);%for top-down viewing
D_mean=nanmean(D_top,3);

%% Shuffle distribution
clear D_sh D_sh_dist D_sh_mean D_sh_perm
binSize=0.1;%mm
xBins= [floor(min(x)*10)/10: binSize: ceil(max(x)*10)/10];
yBins= [floor(min(y)*10)/10: binSize: ceil(max(y)*10)/10];
zBins= [floor(min(z)*10)/10: binSize: ceil(max(z)*10)/10];

switch retinotopy
    
    case 'azimuth'
        
    for perm=1:100
        switch area_loc
            case 'hva'
            randInd=randperm(length(azi_ret)); % shuffle the injection IDs
            azi_ret_perm=azi_ret(randInd);
            
            case 'wb'
            randInd=randperm(16); % shuffle the injection IDs
            azi_ret_perm=azi_ret([1:5,randInd+5]); % only injection#>5 shuffled
        end
        
        D_sh= nan(length(xBins),length(yBins),length(zBins)); % initialize with nans

        for icell=1:length(x)
            xi= find((x(icell)>xBins),1,'last');
            yi= find((y(icell)>yBins),1,'last');
            zi= find((z(icell)>zBins),1,'last');
            if ~isempty(xi) && ~isempty(yi) && ~isempty(zi)
                norm_cell=1/tot_area(C.InjNum(cell_ID(icell))); % per cell weight, normalized by total cells in area per inj    
                norm_ret=(azi_ret_perm(C.InjNum(cell_ID(icell)))-30)/30; %retinotopy normalized to -1 to 1 (range 0 to 60deg)
                D_sh(xi,yi,zi)= nansum([D_sh(xi,yi,zi), norm_cell*norm_ret]); % norm cell count * norm retinotopy
            else
            ;
            end         
        end
        
        D_sh_perm(:,:,:,perm)=D_sh;
        D_sh_top=permute(D_sh,[2,1,3]);%for top-down viewing
        D_sh_mean(:,:,perm)=nanmean(D_sh_top,3);
        
    end
   
        
    case 'altitude'
        
    for perm=1:100        
        switch area_loc
            case 'hva'
            randInd=randperm(length(azi_ret)); % shuffle the injection IDs
            alt_ret_perm=alt_ret(randInd);
            
            case 'wb'
            randInd=randperm(16); % shuffle the injection IDs
            alt_ret_perm=alt_ret([1:5,randInd+5]); % only injection#>5 shuffled
        end
        
        D_sh= nan(length(xBins),length(yBins),length(zBins));% initialize with nans
        
        for icell=1:length(x)
            xi= find((x(icell)>xBins),1,'last');
            yi= find((y(icell)>yBins),1,'last');
            zi= find((z(icell)>zBins),1,'last');
            if ~isempty(xi) && ~isempty(yi) && ~isempty(zi)
                norm_cell=1/tot_area(C.InjNum(cell_ID(icell))); % per cell weight, normalized by total cells in area per inj    
                norm_ret=((alt_ret_perm(C.InjNum(cell_ID(icell)))+10)-12.5)/12.5; %retinotopy normalized to -1 to 1 (range -10 to 15deg)
                D_sh(xi,yi,zi)= nansum([D_sh(xi,yi,zi), norm_cell*norm_ret]); % norm cell count * norm retinotopy
            else
            ;
            end
        end
        
        D_sh_perm(:,:,:,perm)=D_sh;
        D_sh_top=permute(D_sh,[2,1,3]);%for top-down viewing
        D_sh_mean(:,:,perm)=nanmean(D_sh_top,3);
    
    end

end

%% Count biased pixels (based on z-score)
clear mean_vals std_vals z_score

% in 3D
for i = 1:size(D,1)
  for j = 1:size(D,2)
    for k = 1:size(D,3)
       
    mean_vals(i,j,k) = nanmean(reshape(D_sh_perm(i,j,k,:),1,[]));
    std_vals(i,j,k) = nanstd(reshape(D_sh_perm(i,j,k,:),1,[]));
    z_score(i,j,k) = (D(i,j,k) - mean_vals(i,j,k))./std_vals(i,j,k);
    
    end
  end
end

pix_b_01=100*sum(abs(z_score(:))>1.65)./sum(~isnan(D(:)));
pix_b_05=100*sum(abs(z_score(:))>1.96)./sum(~isnan(D(:)));

pix_all=sum(~isnan(D(:)));

T.pix_count{iarea}=pix_all;

cm=RedGrayBlue(100,0.8);

z_score_top=permute(z_score,[2,1,3]);%for top-down viewing
z_score_cor=permute(z_score,[3,1,2]);%for coronal viewing
z_score_sag=permute(z_score,[3,2,1]);%for saggital viewing

z_score_top_mean=nanmean(z_score_top,3);
z_score_top_mean_exp=sign(z_score_top_mean).*exp(abs(z_score_top_mean));

z_score_cor_mean=nanmean(z_score_cor,3);
z_score_cor_mean_exp=sign(z_score_cor_mean).*exp(abs(z_score_cor_mean));

z_score_sag_mean=nanmean(z_score_sag,3);
z_score_sag_mean_exp=sign(z_score_sag_mean).*exp(abs(z_score_sag_mean));

figure(101);
    subplot(1,3,1)
    h=imagesc(rot90(z_score_top_mean,2),[-1.5 1.5]);%[-2.5 2.5]
    axis equal; colorbar; box off; axis off
    title([area ' top-down z-score map ' retinotopy])
    colormap(cm) 
    set(h,'alphadata',~isnan(rot90(z_score_top_mean_exp,2)))% Set nan values to transparent: 
    set(gca, 'color', 'white')%     set(gca, 'color', 'none')

subplot(1,3,2)
    h=imagesc(fliplr(z_score_cor_mean),[-1.5 1.5]);%[-2.5 2.5]
    axis equal; colorbar; box off; axis off
    title([area ' coronal z-score map ' retinotopy])
    colormap(cm) 
    set(h,'alphadata',~isnan(fliplr(z_score_cor_mean_exp)))% Set nan values to transparent: 
    set(gca, 'color', 'white')%     set(gca, 'color', 'none')

subplot(1,3,3)
    h=imagesc(fliplr(z_score_sag_mean),[-1.5 1.5]);%[-2.5 2.5]
    axis equal; colorbar; box off; axis off
    title([area ' saggital z-score map ' retinotopy])
    colormap(cm) 
    set(h,'alphadata',~isnan(fliplr(z_score_sag_mean_exp)))% Set nan values to transparent: 
    set(gca, 'color', 'white')%     set(gca, 'color', 'none')

drawnow

% pause

    switch retinotopy
        case 'azimuth'
        T.azi_bias_01{iarea}=pix_b_01;    
        T.azi_bias_05{iarea}=pix_b_05;    

        case 'altitude'
        T.alt_bias_01{iarea}=pix_b_01;    
        T.alt_bias_05{iarea}=pix_b_05;    
        
    end
    

end
end