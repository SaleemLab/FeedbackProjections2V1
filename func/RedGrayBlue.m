function cmap=RedGrayBlue(m,g)
% m is the length of the colormap
% g is gray level
if nargin<1
    m=100;
end
if nargin<2
    g=0.8;
end
map(:,1)=[linspace(0,g,2*m+1),linspace(g,1,2*m+1)]';
map(:,2)=[linspace(0,g,2*m+1),linspace(g,0,2*m+1)]';
map(:,3)=[linspace(1,g,2*m+1),linspace(g,0,2*m+1)]';
cmap(:,1)=map(1:2:end,1);
cmap(:,2)=map(1:2:end,2);
cmap(:,3)=map(1:2:end,3);

end
