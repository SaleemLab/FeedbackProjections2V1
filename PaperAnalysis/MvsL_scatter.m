datafolder='';

load([datafolder filesep 'MLtable']); 

Allengroups=cellstr(AllenData.Area);
groups=cellstr(MLtable.Area);

M_all=[MLtable.Ms6,MLtable.Md6,MLtable.Md7,MLtable.Md7_2,MLtable.Md12,MLtable.Md13,MLtable.Cd11]';
L_all=[MLtable.Ls6,MLtable.Ls7,MLtable.Ld7,MLtable.Ld8,MLtable.Ld9,MLtable.Ld11,MLtable.Ld12,MLtable.Ld13,MLtable.Cs7]';

% scatter plot of M vs L
M_med=nanmedian(M_all,1);
L_med=nanmedian(L_all,1);
err_M=std(M_all,1,'omitnan')/sqrt(size(M_all,1));
err_L=std(L_all,1,'omitnan')/sqrt(size(L_all,1));

figure;
scatter(L_med,M_med,20,[0,0,1],'filled'); hold on
eb(1) = errorbar(L_med,M_med,err_L, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(L_med,M_med,err_M, 'vertical', 'LineStyle', 'none');
set(eb, 'color',[.75,.75,.75] , 'LineWidth', 0.5)
axis equal square
u=refline(1,0);
u.Color=[.7,.7,.7];
box off

for i=1:length(groups)
text(L_med(i)+0.3,M_med(i),groups(i), 'interpreter','none','fontsize',8)
end

xlim([0 22])
ylim([0 22])

ylabel('Medial injections (median)','fontsize',12)
xlabel('Lateral injections (median)','fontsize',12)

%inset
figure;
scatter(L_med,M_med,20,[0,0,1],'filled'); hold on
eb(1) = errorbar(L_med,M_med,err_L, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(L_med,M_med,err_M, 'vertical', 'LineStyle', 'none');
set(eb, 'color',[.75,.75,.75] , 'LineWidth', 0.5)
axis equal square
u=refline(1,0);
u.Color=[.7,.7,.7];
box off

for i=1:length(groups)
text(L_med(i)+0.3,M_med(i),groups(i), 'interpreter','none','fontsize',8)
end

xlim([0 5])
ylim([0 5])

ylabel('Medial injections (median)','fontsize',12)
xlabel('Lateral injections (median)','fontsize',12)


%% Wilcoxon ranksum 

        for iarea=1:size(MLtable,1)
            [p,h] = ranksum(M_all(:,iarea),L_all(:,iarea)); 
            ranksum_p(iarea)=p; % first column is area names
            ranksum_h(iarea)=h; % first column is area names
        end
    
% plot
cm=RedWhiteBlue;
    figure;
    subplot(1,2,1)
    imagesc(log10(ranksum_p')); title('ranksum log p-value'); colormap(gca,flipud(cm)); colorbar
    set(gca,'ytick',[1:29],'yticklabel',groups)
    set(gca,'xtick',[1:6],'xticklabel','M vs L'); xtickangle(45)
    
    subplot(1,2,2)
    imagesc(ranksum_h'); title('ranksum significance'); colormap(gca,gray); colorbar
    set(gca,'ytick',[1:29],'yticklabel',groups)
    set(gca,'xtick',[1:6],'xticklabel','M vs L'); xtickangle(45)

%% Unpaired ttest   

        for iarea=1:size(MLtable,1)
            [h,p] = ttest2(M_all(:,iarea),L_all(:,iarea)); 
            ttest2_p(iarea)=p; % first column is area names
            ttest2_h(iarea)=h; % first column is area names
        end
    
    % plot
    figure;
    subplot(1,2,1)
    imagesc(log10(ttest2_p')); title('ttestt2 log p-value'); colormap(gca,flipud(cm)); colorbar
    set(gca,'ytick',[1:29],'yticklabel',groups)
    set(gca,'xtick',[1:6],'xticklabel','M vs L'); xtickangle(45)
    
    subplot(1,2,2)
    imagesc(ttest2_h'); title('ttest2 significance'); colormap(gca,gray); colorbar
    set(gca,'ytick',[1:29],'yticklabel',groups)
    set(gca,'xtick',[1:6],'xticklabel','M vs L'); xtickangle(45)
    
%     % save table and figure
%     writetable(ttest2_p,[datafolder filesep 'ttest2_p.xlsx']);
%     savefig([datafolder filesep 'ttest2_result.fig']);    
