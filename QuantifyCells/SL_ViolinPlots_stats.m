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

%% Violin plots

figure;
    h1=violin(AllenData','bw',2,'jitter',0.75,'xlabel',Allengroups');
    set(h1,'FaceColor',[0.62,0.79,0.88],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5) %'LineWidth',1.5
    ax=gca;
    set(gca, 'box','off','TickLength',[0.01 0.01],'TickDir','out','fontsize',10,'color','none','xtick',[1:30]);
    xticklabels(Allengroups)
    ylim([0 10])
    xlabel('Areas')
    ylabel('Relative volume (mm^3)')
    title('Allen Brain Institute Data (V1 input connectivity)','fontsize',20)
    
figure;
    h2=violin(Inj,'bw',2,'jitter',0.75,'xlabel',groups'); 
    set(h2,'FaceColor',[0.62,0.79,0.88],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5) %'LineWidth',1.5
    set(gca, 'box','off','TickLength',[0.01 0.01],'TickDir','out','fontsize',10,'color','none','xtick',[1:30]);
    xticklabels(groups)
    xlabel('Areas','fontsize',15)
    ylabel('Normalized cell count (%)','fontsize',15)
    title('CTB Injection Data (Retrograde)','fontsize',20)

% Log version
logInj=log(Inj_new);
figure;
    h2=violin(logInj,'bw',2,'jitter',0.75,'xlabel',groups'); 
    set(h2,'FaceColor',[0.62,0.79,0.88],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5) %'LineWidth',1.5
    set(gca, 'box','off','TickLength',[0.01 0.01],'TickDir','out','fontsize',10,'color','none','xtick',[1:30]);
    xticklabels(groups)
    ylim([-5 20])
    xlabel('Areas','fontsize',15)
    ylabel('Log normalized cell count (%)','fontsize',15)
    title('CTB Injection Data (Retrograde)','fontsize',20)
    
%% Plot Medial vs Lateral
% All areas
figure;
    h3=violin(M_all,'bw',2,'jitter',0.75,'x',[0.5 2.5 4.5 6.5 8.5 10.5 12.5 14.5 16.5 18.5 20.5 22.5 24.5 26.5 28.5 30.5 32.5 34.5 36.5 38.5 40.5 42.5]);
    set(h3,'FaceColor',[1.0,0.549,0],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
    h4=violin(L_all,'bw',2,'jitter',0.75,'x',[1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43]);
    set(h4,'FaceColor',[0.196,0.804,0.196],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
    set(gca, 'box','off','TickLength',[0.01 0.01],'TickDir','out','fontsize',10,'color','none','xtick',[0.75:2:42.75]);
    xticklabels(groups')
    ylim([-2 45])
    xlabel('Areas','fontsize',15)
    ylabel('Relative cell count (%)','fontsize',15)
    title('CTB injection, all areas','fontsize',20)
    
% Excluding visual areas 
figure;
    h3=violin(M_all(:,1:13),'bw',2,'jitter',0.75,'x',[0.5 2.5 4.5 6.5 8.5 10.5 12.5 14.5 16.5 18.5 20.5 22.5 24.5]);
    set(h3,'FaceColor',[1.0,0.549,0],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
    h4=violin(L_all(:,1:13),'bw',2,'jitter',0.75,'x',[1 3 5 7 9 11 13 15 17 19 21 23 25]);
    set(h4,'FaceColor',[0.196,0.804,0.196],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
    set(gca, 'box','off','TickLength',[0.01 0.01],'TickDir','out','fontsize',10,'color','none','xtick',[0.75:2:24.75]);
    xticklabels(groups(1:13)')
    ylim([-2 45])
    xlabel('Areas','fontsize',15)
    ylabel('Relative cell count (%)','fontsize',15)
    title('CTB injection, non-visual areas','fontsize',20)
 
% HVAs
figure;
    h3=violin(M_all(:,14:end),'bw',2,'jitter',0.75,'x',[2 4 6 8 10 12 14 16 18]);
    set(h3,'FaceColor',[1.0,0.549,0],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
    h4=violin(L_all(:,14:end),'bw',2,'jitter',0.75,'x',[2.5 4.5 6.5 8.5 10.5 12.5 14.5 16.5 18.5]);
    set(h4,'FaceColor',[0.196,0.804,0.196],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
    set(gca, 'box','off','TickLength',[0.01 0.01],'TickDir','out','fontsize',10,'color','none','xtick',[2.25:2:18.25]);
    xticklabels(groups(14:end)')
    ylim([-2 45])
    xlim([1 20])
    xlabel('Areas','fontsize',15)
    ylabel('Relative cell count (%)','fontsize',15)
    title('CTB injection, HVAs','fontsize',20)    
    
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

xlim([0 20])
ylim([0 20])

ylabel('Medial injections (median)','fontsize',12)
xlabel('Lateral injections (median)','fontsize',12)

%% Plot Anterior vs Posterior

% All areas
figure;
subplot(2,1,1)
h3=violin(A_all,'bw',2,'jitter',0.75,'x',[0.5 2.5 4.5 6.5 8.5 10.5 12.5 14.5 16.5 18.5 20.5 22.5 24.5 26.5 28.5 30.5 32.5 34.5 36.5 38.5 40.5 42.5]);
set(h3,'FaceColor',[0.306,0.702,0.827],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
h4=violin(P_all,'bw',2,'jitter',0.75,'x',[1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43]);
set(h4,'FaceColor',[0.969,0.408,0.631],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
set(gca, 'box','off','TickLength',[0.01 0.01],'TickDir','out','fontsize',10,'color','none','xtick',[0.75:2:42.75]);
xticklabels(groups')
ylim([-2 45])
xlabel('Areas','fontsize',15)
ylabel('Relative cell count (%)','fontsize',15)
title('CTB injection, all areas','fontsize',20)
    
% Excluding visual areas 
figure;
    h3=violin(A_all(:,1:13),'bw',2,'jitter',0.75,'x',[0.5 2.5 4.5 6.5 8.5 10.5 12.5 14.5 16.5 18.5 20.5 22.5 24.5]);
    set(h3,'FaceColor',[0.306,0.702,0.827],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
    h4=violin(P_all(:,1:13),'bw',2,'jitter',0.75,'x',[1 3 5 7 9 11 13 15 17 19 21 23 25]);
    set(h4,'FaceColor',[0.969,0.408,0.631],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
    set(gca, 'box','off','TickLength',[0.01 0.01],'TickDir','out','fontsize',10,'color','none','xtick',[0.75:2:24.75]);
    xticklabels(groups(1:13)')
    ylim([-2 45])
    xlabel('Areas','fontsize',15)
    ylabel('Relative cell count (%)','fontsize',15)
    title('CTB injection, non-visual areas','fontsize',20)

% HVAs
figure;
    subplot(2,1,1)
    h3=violin(A_all(:,14:end),'bw',2,'jitter',0.75,'x',[2 4 6 8 10 12 14 16 18]);
    set(h3,'FaceColor',[0.306,0.702,0.827],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
    h4=violin(P_all(:,14:end),'bw',2,'jitter',0.75,'x',[2.5 4.5 6.5 8.5 10.5 12.5 14.5 16.5 18.5]);
    set(h4,'FaceColor',[0.969,0.408,0.631],'FaceAlpha',0.3,'EdgeColor',[0,0,0],'LineWidth',1.5)
    set(gca, 'box','off','TickLength',[0.01 0.01],'TickDir','out','fontsize',10,'color','none','xtick',[2.25:2:18.25]);
    xticklabels(groups(14:end)')
    ylim([-2 45])
    xlim([1 20])
    xlabel('Areas','fontsize',15)
    ylabel('Relative cell count (%)','fontsize',15)
    title('CTB injection, HVAs','fontsize',20)
 
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
    sig=mult_mat<0.05;
    figure;
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
    
%% T-test and Wilcoxon ranksum

% Collecting partial data from master table
Md = [MLtable.Md5,MLtable.Md6,MLtable.Md6_2,MLtable.Md7,MLtable.Md7_2,MLtable.Md12,MLtable.Md13];
Ms = [MLtable.Ms5,MLtable.Ms6];
Ld = [MLtable.Ld5,MLtable.Ld7,MLtable.Ld8,MLtable.Ld9,MLtable.Ld11,MLtable.Ld12,MLtable.Ld13];
Ls = [MLtable.Ls5,MLtable.Ls6,MLtable.Ls7];
M_all = [Md,Ms];
L_all = [Ld,Ls];
C_all=[MLtable.Cs7,MLtable.Cd11];

inj_cat =7; % category of injections
inj_num = 7;% max number of injections per injection category
inj_names = {'Md','Ms','Ld','Ls','M all','L all','C_all'};

% Collate to a single matrix
MLmat = NaN(size(MLtable,1),inj_num,inj_cat);
MLmat(:,1:size(Md,2),1)= Md;
MLmat(:,1:size(Ms,2),2)= Ms;
MLmat(:,1:size(Ld,2),3)= Ld;
MLmat(:,1:size(Ls,2),4)= Ls;
MLmat(:,1:size(M_all,2),5)= M_all;
MLmat(:,1:size(L_all,2),6)= L_all;
MLmat(:,1:size(C_all,2),7)= C_all;

% Compare M to L, D to S
% Q: per area, is there a preference for projecting to M or L? comparisons to be made = 'Ms vs Ls', 'Md vs 
% Ld', 'M_all vs L_all'

% Generate an empty table for test results
    Area=MLtable.Area;
    MLs=NaN(size(MLtable,1),1);
    MLd=NaN(size(MLtable,1),1);
    MLall=NaN(size(MLtable,1),1);
    
% Comparisons to be made
    comp = {'Ms vs Ls','Md vs Ld','Mall vs Lall'};
    comb =[2,4; 1,3; 5,6]; % translating to combinations in third dimention of MLmat matrix

% Wilcoxon ranksum
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
    
% save table and figure
%     writetable(ranksum_p,[datafolder filesep 'ranksum_p.xlsx']);
%     savefig([datafolder filesep 'ranksum_result.fig']);

% Unpaired ttest
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
    
%     % save table and figure
%     writetable(ttest2_p,[datafolder filesep 'ttest2_p.xlsx']);
%     savefig([datafolder filesep 'ttest2_result.fig']);
