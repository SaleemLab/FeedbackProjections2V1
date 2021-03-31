% For clustering cell counts per area from CTB injection in V1
% 2019-12 MM

%% Load data
% HVA data
dir='';
% cell count normalized
load([dir filesep 'all_inj_rel_hva_new.mat']);
% names
load([dir filesep 'all_inj_names_new.mat']);
load([dir filesep 'hva_names.mat']);
% injecion site coords transformed to Allen coords
load([dir filesep 'injMLcoords_transf.mat']) 
load([dir filesep 'injAPcoords_transf.mat']) 
% allen CCF boundary image
load([dir filesep 'allen_borders.mat']);

data=data'; % data needs to be sample-by-variable format

%% Plot injection sites
figure;
imagesc(im_bound); 
colormap(flipud(gray))
axis equal tight
hold on
scatter(ml,ap, '.r');
for i=1:size(data,1)
    text(ml(i),ap(i),inj{i}, 'interpreter','none')
end

%% Plot injection sites (left hemi)
figure;
imagesc(im_bound); 
colormap(flipud(gray))
axis equal tight
hold on
scatter(ml-570,ap, '.r');
for i=1:size(data,1)
    text(ml(i)-570,ap(i),inj{i}, 'interpreter','none')
end

%% PCA
[coeff,score,latent,tsquared,explained] = pca(data);

figure;
pareto(explained)
title('explained variance of Principal Components')

figure;
scatter(score(:,1),score(:,2))
hold on;
scatter(score(4,1),score(4,2),'r');
scatter(score(9,1),score(9,2),'g');
scatter(score(sort_idx(2:8),1),score(sort_idx(2:8),2),'m');
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

for i=1:size(data,1)
text(score(i,1),score(i,2),inj{i})
end

figure;
scatter3(score(:,1),score(:,2),score(:,3));
hold on;
scatter3(score(4,1),score(4,2),score(4,3),'r');
scatter3(score(9,1),score(9,2),score(9,3),'g');
scatter3(score(sort_idx(2:8),1),score(sort_idx(2:8),2),score(sort_idx(2:8),3),'m');

biplot(coeff(:,1:3),'scores',score(:,1:3),'varlabels',areas);

%% k-means
% % using all data, 2 clusters
% [cidx2,cmeans2] = kmeans(data,2,'dist','sqeuclidean');
% [silh2,h] = silhouette(data,cidx2,'sqeuclidean');

% using first 4PCs, 2 clusters
[~,sort_idx] = sort(ml);
M_point=1;
L_point=length(ml)-2;
seed = [data(sort_idx(M_point),:); data(sort_idx(L_point),:)]; % 2nd most L and 2nd most M
[cidx2,cmeans2] = kmeans(data(:,:),2,'dist','sqeuclidean','Display','iter', 'Start', seed);

 
% Plot 2 clusters
figure;
imagesc(im_bound); 
colormap(flipud(gray))
axis equal tight
hold on

for i=1:length(cidx2)
    switch cidx2(i)
        case 1
            col='m';
        case 2
            col='c';
        case 3
            col='g';
        case 4
            col='b';
        otherwise
            col=[0.5 0.5 0.5];
    end
    
scatter(ml(i),ap(i), 30, 'o', 'MarkerFaceColor', col,  'MarkerEdgeColor', 'none');
end
scatter(ml(sort_idx(M_point)),ap(sort_idx(M_point)), 30, 'o', 'MarkerEdgeColor', 'k', 'LineWidth',1); hold on
scatter(ml(sort_idx(L_point)),ap(sort_idx(L_point)), 30, 'o', 'MarkerEdgeColor', 'k', 'LineWidth',1);
title('2 cluster k-means on full data with seeds')
box off
axis off

%% replot scatter with k-means clusters

figure;

scatter(score(cidx2==1,1),score(cidx2==1,2),40, 'o', 'MarkerFaceColor', 'm',  'MarkerEdgeColor', 'w', 'LineWidth',1); hold on
scatter(score(cidx2==2,1),score(cidx2==2,2),40, 'o', 'MarkerFaceColor', 'c',  'MarkerEdgeColor', 'w', 'LineWidth',1);hold on
scatter(score(sort_idx(M_point),1),score(sort_idx(M_point),2), 40, 'o', 'MarkerEdgeColor', 'k', 'LineWidth',1); hold on
scatter(score(sort_idx(L_point),1),score(sort_idx(L_point),2), 40, 'o', 'MarkerEdgeColor', 'k', 'LineWidth',1);

xlabel('1st Principal Component')
ylabel('2nd Principal Component')

% for i=1:size(data,1)
% text(score(i,1),score(i,2),inj{i})
% end

%% Plot 3 clusters
% using first 4PCs, 3 clusters
[cidx2,cmeans2] = kmeans(data(:,:),2,'dist','sqeuclidean','Display','iter', 'Start', seed);
[cidx3,cmeans3] = kmeans(score(:,1:4),3,'dist','sqeuclidean');
[silh3,h] = silhouette(score(:,1:4),cidx3,'sqeuclidean');

figure;
imagesc(im_bound); 
colormap(flipud(gray))
axis equal tight
hold on

for i=1:length(cidx3)
    if cidx3(i)==1
        col='r';
    elseif cidx3(i)==2
        col='b';
    else 
        col='g';
    end
    
scatter(ml(i),ap(i), 30, 'o', 'MarkerFaceColor', col,  'MarkerEdgeColor', 'none');

end

title('3 cluster k-means with 4PCs')


%% Plot 3 clusters with same seed
% using first 4PCs, 3 clusters
[~,sort_idx] = sort(ml);
M_point=1;
L_point=length(ml)-2;
mid_point=13;
seed3 = [data(sort_idx(M_point),:); data(sort_idx(L_point),:); data(sort_idx(mid_point),:)]; % 2nd most L and 2nd most M

[cidx3,cmeans3] = kmeans(data(:,:),3,'dist','sqeuclidean','Display','iter', 'Start', seed3);
[silh3,h] = silhouette(score(:,1:4),cidx3,'sqeuclidean');

figure;
imagesc(im_bound); 
colormap(flipud(gray))
axis equal tight
hold on

for i=1:length(cidx3)
    if cidx3(i)==1
        col='m';
    elseif cidx3(i)==2
        col='c';
    else 
        col='r';
    end
    
scatter(ml(i),ap(i), 30, 'o', 'MarkerFaceColor', col,  'MarkerEdgeColor', 'none');

end
scatter(ml(sort_idx(M_point)),ap(sort_idx(M_point)), 30, 'o', 'MarkerEdgeColor', 'k', 'LineWidth',1); hold on
scatter(ml(sort_idx(L_point)),ap(sort_idx(L_point)), 30, 'o', 'MarkerEdgeColor', 'k', 'LineWidth',1);hold on
scatter(ml(sort_idx(mid_point)),ap(sort_idx(mid_point)), 30, 'o', 'MarkerEdgeColor', 'k', 'LineWidth',1);
title('3 cluster k-means with full data with seed')

%% replot scatter with k-means clusters

figure;

scatter(score(cidx3==1,1),score(cidx3==1,2),40, 'o', 'MarkerFaceColor', 'm',  'MarkerEdgeColor', 'w', 'LineWidth',1); hold on
scatter(score(cidx3==2,1),score(cidx3==2,2),40, 'o', 'MarkerFaceColor', 'c',  'MarkerEdgeColor', 'w', 'LineWidth',1);hold on
scatter(score(cidx3==3,1),score(cidx3==3,2),40, 'o', 'MarkerFaceColor', 'r',  'MarkerEdgeColor', 'w', 'LineWidth',1);hold on

scatter(score(sort_idx(M_point),1),score(sort_idx(M_point),2), 40, 'o', 'MarkerEdgeColor', 'k', 'LineWidth',1); hold on
scatter(score(sort_idx(L_point),1),score(sort_idx(L_point),2), 40, 'o', 'MarkerEdgeColor', 'k', 'LineWidth',1); hold on
scatter(score(sort_idx(mid_point),1),score(sort_idx(mid_point),2), 40, 'o', 'MarkerEdgeColor', 'k', 'LineWidth',1);

xlabel('1st Principal Component')
ylabel('2nd Principal Component')

% for i=1:size(data,1)
% text(score(i,1),score(i,2),inj{i})
% end

