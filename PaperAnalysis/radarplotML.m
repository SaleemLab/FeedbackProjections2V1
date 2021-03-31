% Load data 
% rows are samples, column are areas
% polar plot
figure;
rho = [0.540409097,0.572503291,0.643619123,0.334528076,0.388908389,0.330258957,0.511217283,0.361178046,0.342997626;...
    0.459590903,0.427496709,0.356380877,0.665471924,0.611091611,0.669741043,0.488782717,0.638821954,0.657002374];

h=polarplot(rho,'-or');
hP=h.Parent;                % retrieve polar axes handle (parent of polar plot)
hP.RTickLabel={'PM','AM','A','RL','AL','L','LI','POR','PL'};

figure;
order=[3,2,1,9,8,7,6,5,4];
radarplot(data_norm(:,order),{'A','AM','PM','PL','POR','LI','L','AL','RL'},[],[],[],4);
legend('Medial injections','Lateral injections')
title('normalized cell count (%)')

data_norm_sum=sum(data_norm,1);

for i=1:2
    for j=1:9
        data_norm_mod(i,j)=data_norm(i,j)/data_norm_sum(j);
    end
end

figure;
order=[3,2,1,9,8,7,6,5,4];
radarplot(data_norm_mod(:,order),{'A','AM','PM','PL','POR','LI','L','AL','RL'},[],[],[],4);
legend('Medial injections','Lateral injections')
title('normalized cell count M vs L')

%% RadarPlot
data=data';
%index of M and L categories
indM=[1:7,19,21];
indL=[8:18,20];
meanM=mean(data(indM,:),1); 
meanL=mean(data(indL,:),1); 

% for i=1:size(data,2)
%     modM(i)=log(meanM(i)/meanL(i));
%     modL(i)=log(meanL(i)/meanM(i));
% end

for i=1:size(data,2)
    modM(i)=meanM(i)/(meanL(i)+meanM(i));
    modL(i)=meanL(i)/(meanL(i)+meanM(i));
end

mat=[modM(order);modL(order)];

figure;
order=[3,2,1,9,8,7,6,5,4];
radarplot(mat,{'A','AM','PM','PL','POR','LI','L','AL','RL'},[],[],[],4);
legend('Medial injections','Lateral injections')
title('normalized cell count M vs L')  












































