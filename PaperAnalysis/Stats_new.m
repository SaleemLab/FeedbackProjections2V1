%% Load data
datafolder='';

load([datafolder filesep 'AllenData_FB']); 
load([datafolder filesep 'MLtable']); 
load([datafolder filesep 'APtable']); 

Allengroups=cellstr(AllenData.Area);
groups=cellstr(MLtable.Area);

AllenData=table2array(AllenData(:,2:end));
Inj=table2array(MLtable(:,2:end))';
Inj_new=Inj;
Inj_new(Inj_new==0)=NaN; % replace zeros with NaNs

M_all=[MLtable.Md4,MLtable.Md5_2,MLtable.Ms6,MLtable.Md6,MLtable.Md7,MLtable.Md7_2,MLtable.Md12,MLtable.Md13]';
L_all=[MLtable.Ls4,MLtable.Ls6,MLtable.Ls7,MLtable.Ld5,MLtable.Ld5_2,MLtable.Ld7,MLtable.Ld8,MLtable.Ld9,MLtable.Ld11,MLtable.Ld12,MLtable.Ld13]';


%% Normality test
figure;
for i=1:size(ML_all,2)
    x=ML_all(:,i);
    [h,p] = kstest(x);%Kolmogorov-Smirnov test
    ML_ks_h(i)=h;
    ML_ks_p(i)=p;

    subplot(4,5,i)
    cdfplot(x)
    hold on
    x_values = linspace(min(x),max(x));
    plot(x_values,normcdf(x_values,0,1),'r-')
    legend('Empirical CDF','Standard Normal CDF','Location','best')
end

figure;
plot(ML_ks_p)

%% Multiple Comparisons
figure;
[p,tbl,stats] = kruskalwallis(Inj,groups,'on'); 
figure;
c = multcompare(stats, 'Display','on');   
 
mult_mat=NaN(22,22);
for i=1:22
    start=find(c(:,1)==i,1);
    mult_mat(i,i+1:22)=c(start:start+21-i,6);
end
cm = RedWhiteBlue;

% make grid
    [X,Y] = meshgrid(0.5:length(groups)+0.5,0.5:length(groups)+0.5);
    
% plot with imagesc
    sl=log10(.05); % significance level=.05
    figure;
    sig=mult_mat<0.05;
    imagesc(sig); title('multcompare log p-value'); colormap(gca,flipud(gray)); colorbar; axis square; axis tight
    set(gca,'ytick',[1:22],'yticklabel',groups)
    set(gca,'xtick',[1:22],'xticklabel',groups); xtickangle(90)
    set(gca,'TickDir','out','box','off','color','none'); hold on; 
    line(X,Y, 'Color',[.7 .7 .7]); hold on
    line(Y,X, 'Color',[.7 .7 .7]); hold off
    
    figure;
    mult_mat_new=rot90(mult_mat,2);
    sig_new=mult_mat_new<0.05;
    imagesc(sig_new); title('multcompare log p-value'); colormap(gca,flipud(gray)); colorbar; axis square; axis tight
    set(gca,'ytick',[1:22],'yticklabel',flipud(groups))
    set(gca,'xtick',[1:22],'xticklabel',flipud(groups)); xtickangle(90)
    set(gca,'TickDir','out','box','off','color','none'); hold on; 
    line(X,Y, 'Color',[.7 .7 .7]); hold on
    line(Y,X, 'Color',[.7 .7 .7]); hold off

%% 1. Compare M to L, D to S
% Q: per area, is there a preference for projecting to M or L? try ttest2, ranksum function (some other rank tests?) comparisons to be made = 'Ms vs Ls', 'Md vs Ld', 'M_all vs L_all'

% Generate an empty table for test results
    Area=MLtable.Area;
    MLs=NaN(size(MLtable,1),1);
    MLd=NaN(size(MLtable,1),1);
    MLall=NaN(size(MLtable,1),1);
    
% Comparisons to be made
    comp = {'Ms vs Ls','Md vs Ld','Mall vs Lall'};
    comb =[2,4; 1,3; 5,6]; % translating to combinations in third dimention of MLmat matrix


%% Wilcoxon ranksum

ranksum_p=table(Area,MLs,MLd,MLall);
ranksum_h=table(Area,MLs,MLd,MLall);
    
% Run all comparisons
    for icomb=1:size(comb,1)
        for iarea=1:size(MLtable,1)
            X=comb(icomb,1);
            Y=comb(icomb,2);
            [p,h] = ranksum(MLmat(iarea,:,X),MLmat(iarea,:,Y)); 
            ranksum_p{iarea,icomb+1}=p; % first column is area names
            ranksum_h{iarea,icomb+1}=h; % first column is area names
        end
    end
    
% plot
    ranksum_p_mat=table2array(ranksum_p(:,2:end));
    ranksum_h_mat=table2array(ranksum_h(:,2:end));
    
    figure;
    subplot(1,2,1)
    imagesc(log10(ranksum_p_mat)); title('ranksum log p-value'); colormap(gca,flipud(cm)); colorbar
    set(gca,'ytick',[1:29],'yticklabel',Area)
    set(gca,'xtick',[1:6],'xticklabel',comp); xtickangle(45)
    
    subplot(1,2,2)
    imagesc(ranksum_h_mat); title('ranksum significance'); colormap(gca,gray); colorbar
    set(gca,'ytick',[1:29],'yticklabel',Area)
    set(gca,'xtick',[1:6],'xticklabel',comp); xtickangle(45)
    
%% Unpaired ttest

ttest2_p=table(Area,MLs,MLd,MLall);
ttest2_h=table(Area,MLs,MLd,MLall);
     
     % Run all comparisons
    for icomb=1:size(comb,1)
        for iarea=1:size(MLtable,1)
            X=comb(icomb,1);
            Y=comb(icomb,2);
            [h,p] = ttest2(MLmat(iarea,:,X),MLmat(iarea,:,Y)); 
            ttest2_p{iarea,icomb+1}=p; % first column is area names
            ttest2_h{iarea,icomb+1}=h; % first column is area names
        end
    end
    
    % plot
    ttest2_p_mat=table2array(ttest2_p(:,2:end));
    ttest2_h_mat=table2array(ttest2_h(:,2:end));
    
    figure;
    subplot(1,2,1)
    imagesc(log10(ttest2_p_mat)); title('ttestt2 log p-value'); colormap(gca,flipud(cm)); colorbar
    set(gca,'ytick',[1:29],'yticklabel',Area)
    set(gca,'xtick',[1:6],'xticklabel',comp); xtickangle(45)
    
    subplot(1,2,2)
    imagesc(ttest2_h_mat); title('ttest2 significance'); colormap(gca,gray); colorbar
    set(gca,'ytick',[1:29],'yticklabel',Area)
    set(gca,'xtick',[1:6],'xticklabel',comp); xtickangle(45)

%% 2. Per injection, is there a bias in the areas it receives inputs from? 

% Paired t-test
% Create empty matrix for results (ML)
ML_area_pairs_p=NaN(length(Area),length(Area),ML_inj_cat);
ML_area_pairs_h=NaN(length(Area),length(Area),ML_inj_cat);

% Pairwise comparison between areas per injection category
for iinj=1:ML_inj_cat
    for iarea1=1:size(MLmat,1)
        for iarea2=1:iarea1-1
        [h,p,ci,stats] = ttest(MLmat(iarea1,:,iinj),MLmat(iarea2,:,iinj)); 
        if stats.tstat>0
            ML_area_pairs_p(iarea1,iarea2,iinj)=log10(p); % = negative = red
        else
            ML_area_pairs_p(iarea1,iarea2,iinj)=-log10(p); % = positive = blue
        end
            ML_area_pairs_h(iarea1,iarea2,iinj)=h;
        end
    end
    
% make NaN entries go to 0
    N = isnan(ML_area_pairs_p);
    ML_area_pairs_p(N)=0;
    
% make grid
    [X,Y] = meshgrid(0.5:length(Area)+0.5,0.5:length(Area)+0.5);
    
% plot with imagesc
    sl=log10(.01); % significance level=.05
    figure;
    subplot(1,2,1)
    imagesc(ML_area_pairs_p(:,:,iinj),[sl -sl]); title([ML_inj_names(iinj) ' ttest log p-value']); colormap(gca,flipud(cm)); colorbar; axis square; axis tight
    set(gca,'ytick',[1:26],'yticklabel',Area,'YColor','red')
    set(gca,'xtick',[1:26],'xticklabel',Area,'XColor','blue'); xtickangle(90)
    set(gca,'TickDir','out','box','off','color','none'); hold on; 
    line(X,Y, 'Color',[.7 .7 .7]); hold on
    line(Y,X, 'Color',[.7 .7 .7]); hold off
    
    subplot(1,2,2)
    imagesc(ML_area_pairs_h(:,:,iinj)); title([ML_inj_names(iinj) ' ttest significance']); colormap(gca,flipud(gray)); colorbar; axis square; axis tight
    set(gca,'ytick',[1:26],'yticklabel',Area)
    set(gca,'xtick',[1:26],'xticklabel',Area); xtickangle(90)
    set(gca,'TickDir','out','box','off','color','none'); hold on; 
    line(X,Y, 'Color',[.7 .7 .7]); hold on
    line(Y,X, 'Color',[.7 .7 .7]); hold off
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.96, 0.96]);
    
end

























