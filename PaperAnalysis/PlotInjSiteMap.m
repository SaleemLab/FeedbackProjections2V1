%% inj site plot
cd 'your path\GitHub\allenCCF' 
bregma=allenCCFbregma();
atlas_resolution=0.01;
% load boundary
InfoImage.low=imfinfo('your path\CCF_borders_10um.tif');
im=zeros(InfoImage.low.Height,InfoImage.low.Width,length(InfoImage.low),'double');
TifLink = Tiff('your path\CCF_borders_10um.tif', 'r');
for i=1:length(InfoImage.low)
TifLink.setDirectory(i);
im(:,:,i)=TifLink.read();
end
TifLink.close();

% translate coords to allen map coords
ml=abs(inj_coords(2,:))/atlas_resolution+bregma(3);
ap=-inj_coords(1,:)/atlas_resolution+bregma(1);

figure;
imagesc(im); colormap(flipud(gray))
axis equal tight
hold on
scatter(ml,ap, '.r');
for i=1:21 
text(ml(i),ap(i),inj{i}, 'interpreter','none')
end

%%
figure;
imagesc(im); colormap(flipud(gray))
axis equal tight
hold on

for i=1:length(cidx2)
    if cidx2(i)==1
        col='r';
    else
        col='b';
    end
    
scatter(ml(i),ap(i), 30, 'o', 'MarkerFaceColor', col,  'MarkerEdgeColor', 'none');

end

% title('2 cluster k-means with 4PCs')
title('2 cluster k-means with all data')

%%
figure;
imagesc(im); colormap(flipud(gray))
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
% text(ml(i),ap(i),inj{i}, 'interpreter','none', 'Color', col)
end

% title('3 cluster k-means with 4PCs')
title('3 cluster k-means with all data')

%%
figure;
imagesc(im); colormap(flipud(gray))
axis equal tight
hold on

for i=1:22
    switch T1(i)
    case 1
        col='r';
    case 2
        col='b';
    case 3
        col='m';
    case 4
        col='g';
    end
scatter(ml(i),ap(i), 30, 'o', 'MarkerFaceColor', col,  'MarkerEdgeColor', 'none');
% text(ml(i),ap(i),inj{i}, 'interpreter','none', 'Color', col)
end

%%
figure;
imagesc(im); colormap(flipud(gray))
axis equal tight
hold on

group=[2,11,10,20,9,1];

for i=1:22
    if sum(ismember(group,i))==0
        col='b';
    else
        col='r';
    end
scatter(ml(i),ap(i), 30, 'o', 'MarkerFaceColor', col,  'MarkerEdgeColor', 'none');
% text(ml(i),ap(i),inj{i}, 'interpreter','none', 'Color', col)
end
title('pc1 vs pc2')

