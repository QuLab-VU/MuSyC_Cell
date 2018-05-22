%This code is to generate example dose response surfaces in matlab
%Set default plot parameters
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultLineLinewidth',4)
set(0, 'DefaultAxesFontWeight','bold')
set(0, 'DefaultAxesLineWidth',4)
set(0, 'DefaultFigureColor','w')
set(0, 'DefaultTextFontSize',12)
set(0, 'DefaultTextFontWeight','bold')

%Set parameters for surfaces
d1 = logspace(-5,0,1000);
d2 = logspace(-5,0,1000);
h1 = 1;
h2 = 1;
C1 = .005;
C2 = .001;
E0 = .05;
E1 = .005;
E2 = .01;
beta = [0];
alpha = [1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Null surface plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure('color','w')
figure(h)
clf
p = 1;
k = 1;
for i = 1:length(d1)
    for j = 1:length(d2)
        E3 = -beta(p)*min(E1,E2)+min(E1,E2);
        Ed(i,j) = (C1^h1*C2^h2*E0 + d1(i).^h1*C2^h2*E1 + d2(j).^h2*C1^h1*E2 + alpha(k)*d1(i)^h1*d2(j)^h2*E3)/(C1^h1*C2^h2 + d1(i)^h1*C2^h2 + d2(j)^h2*C1^h1 + alpha(k)*d1(i)^h1*d2(j)^h2);
    end
end
%Downloaded from http://www.kennethmoreland.com/color-maps/
map = diverging_map(0:.001:1,[.23,.299,.745],[.706,.016,.15]);
[X,Y] = meshgrid(d1,d2);
C = contour(X,Y,Ed,20);
clf
hold on
surf(X,Y,Ed,'lines','None')
colormap(map)
set(gca,'xscale','log','yscale','log','linewidth',2.,'xticklabels',[],'yticklabels',[],'zticklabels',[],'xtick',logspace(-5,-1,5),'ytick',logspace(-5,-1,5))
grid off
axis([min(d1),max(d1),min(d2),max(d2),-.005,max(Ed(:))])
view([67,46])
xlabel('log([Drug 1])')
ylabel('log([Drug 2])')
zlabel('DIP')
hold on
scatter3(d1,min(d2)*ones(1000,1),Ed(1,:),30,[.7,.7,.7],'filled')
scatter3(min(d1)*ones(1000,1),d2,Ed(:,1),30,'k','filled')
scatter3(min(d1),min(d2),max(Ed(:)),150,'r','filled')
scatter3(max(d1),min(d2),Ed(1,end),150,'r','filled')
scatter3(min(d1),max(d2),Ed(end,1),150,'r','filled')
scatter3(max(d1),max(d2),min(Ed(:)),150,'r','filled')
[~,m_idx] = min((d2-C2).^2);
scatter3(C2,min(d2),Ed(1,m_idx),150,'r','filled')
[~,m_idx] = min((d1-C1).^2);
scatter3(min(d2),C1,Ed(m_idx,1),150,'r','filled')

s = contourdata(C); 
for i = 1:length(s)
    cont = [];
    for j = 1:length(s(i).xdata)
            cont(j) = (C1^h1*C2^h2*E0 + s(i).xdata(j).^h1*C2^h2*E1 + s(i).ydata(j).^h2*C1^h1*E2 + alpha(k)*s(i).xdata(j)^h1*s(i).ydata(j)^h2*E3)/(C1^h1*C2^h2 + s(i).xdata(j)^h1*C2^h2 + s(i).ydata(j)^h2*C1^h1 + alpha(k)*s(i).xdata(j)^h1*s(i).ydata(j)^h2);                    
    end
    %plot3(s(i).xdata,s(i).ydata,cont,'color',map(i,:),'linewidth',1)
    plot3(s(i).ydata,s(i).xdata,cont,'k','linewidth',1)
end
set(gca,'ZGrid','on')
scatter3(ones(1000,1)*min(d1),d2,zeros(1000,1),10,'k')
scatter3(d1,max(d2)*ones(1000,1),zeros(1000,1),10,'k')
scatter3(ones(1000,1)*max(d1),d2,zeros(1000,1),10,'k')
scatter3(d1,min(d2)*ones(1000,1),zeros(1000,1),10,'k')
xlabel([])
ylabel([])
zlabel([])
saveas(h,'Plots/Figure1B.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Single Dose response curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure('color','w')
hold on
set(gca,'xscale','log','linewidth',2.,'xticklabels',[],'yticklabels',[],'xtick',logspace(-5,-1,5))
scatter(d1,Ed(1,:),30,[.7,.7,.7],'filled')
scatter(d2,Ed(:,1),30,'k','filled')
scatter(min(d1),max(Ed(:)),150,'r','filled')
scatter(max(d1),Ed(1,end),150,'r','filled')
scatter(max(d2),Ed(end,1),150,'r','filled')
[~,m_idx] = min((d2-C2).^2);
scatter(C2,Ed(1,m_idx),150,'r','filled')
[~,m_idx] = min((d1-C1).^2);
scatter(C1,Ed(m_idx,1),150,'r','filled')
set(gca,'YGrid','on')
scatter(d1,zeros(1000,1),10,'k','filled')

saveas(h,'Plots/Figure1A.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Different quadrants of synergy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find the bounds for plots
beta = [-.02,.01];
alpha = [.05,10];
mx_Ed = -1
m_Ed = 1
for k = 1:length(alpha)
    for p = 1:length(beta)
        for i = 1:length(d1)
            for j = 1:length(d2)
                E3 = -beta(p)+min(E1,E2);
                Ed(i,j) = (C1^h1*C2^h2*E0 + d1(i).^h1*C2^h2*E1 + d2(j).^h2*C1^h1*E2 + alpha(k)*d1(i)^h1*d2(j)^h2*E3)/(C1^h1*C2^h2 + d1(i)^h1*C2^h2 + d2(j)^h2*C1^h1 + alpha(k)*d1(i)^h1*d2(j)^h2);
            end
        end
        m_Ed = min([m_Ed,min(Ed(:))]);
        mx_Ed = max([mx_Ed,max(Ed(:))]);
    end
end

for k = 1:length(alpha)
    for p = 1:length(beta)
        for i = 1:length(d1)
            for j = 1:length(d2)
                E3 = -beta(p)+min(E1,E2);
                Ed(i,j) = (C1^h1*C2^h2*E0 + d1(i).^h1*C2^h2*E1 + d2(j).^h2*C1^h1*E2 + alpha(k)*d1(i)^h1*d2(j)^h2*E3)/(C1^h1*C2^h2 + d1(i)^h1*C2^h2 + d2(j)^h2*C1^h1 + alpha(k)*d1(i)^h1*d2(j)^h2);
            end
        end        
        h = figure('color','w');
        C = contour(X,Y,Ed,20);
        clf
        s = contourdata(C);        %subplot(length(alpha),length(beta),(k-1)*length(alpha)+p)
        [X,Y] = meshgrid(d1,d2);
        surf(X,Y,Ed,'lines','None')
        set(gca,'xscale','log','yscale','log','linewidth',2.,'xticklabels',[],'yticklabels',[],'zticklabels',[],'xtick',logspace(-5,-1,5),'ytick',logspace(-5,-1,5))
        map = diverging_map(0:.001:1,[.23,.299,.745],[.706,.016,.15]);
        colormap(map)
        caxis([m_Ed,mx_Ed])
        grid off
        axis([min(d1),max(d1),min(d2),max(d2),m_Ed,mx_Ed])
        view([67,46])
        hold on
        map = diverging_map(0:1/length(s):1,[.23,.299,.745],[.706,.016,.15]);
        scatter3(d1,min(d2)*ones(1000,1),Ed(1,:),30,[.7,.7,.7],'filled')
        scatter3(min(d1)*ones(1000,1),d2,Ed(:,1),30,'k','filled')
        scatter3(min(d1),min(d2),max(Ed(:)),150,'r','filled')
        scatter3(max(d1),min(d2),Ed(1,end),150,'r','filled')
        scatter3(min(d1),max(d2),Ed(end,1),150,'r','filled')
        scatter3(max(d1),max(d2),Ed(end,end),150,'r','filled')
        [~,m_idx] = min((d2-C2).^2);
        scatter3(C2,min(d2),Ed(1,m_idx),150,'r','filled')
        [~,m_idx] = min((d1-C1).^2);
        scatter3(min(d2),C1,Ed(m_idx,1),150,'r','filled')

        s = contourdata(C); 
        for i = 1:length(s)
            cont = [];
            for j = 1:length(s(i).xdata)
                    cont(j) = (C1^h1*C2^h2*E0 + s(i).xdata(j).^h1*C2^h2*E1 + s(i).ydata(j).^h2*C1^h1*E2 + alpha(k)*s(i).xdata(j)^h1*s(i).ydata(j)^h2*E3)/(C1^h1*C2^h2 + s(i).xdata(j)^h1*C2^h2 + s(i).ydata(j)^h2*C1^h1 + alpha(k)*s(i).xdata(j)^h1*s(i).ydata(j)^h2);                    
            end
            %plot3(s(i).xdata,s(i).ydata,cont,'color',map(i,:),'linewidth',1)
            plot3(s(i).ydata,s(i).xdata,cont,'k','linewidth',1)
        end
        set(gca,'ZGrid','on')
        scatter3(ones(1000,1)*min(d1),d2,zeros(1000,1),10,'k')
        scatter3(d1,max(d2)*ones(1000,1),zeros(1000,1),10,'k')
        scatter3(ones(1000,1)*max(d1),d2,zeros(1000,1),10,'k')
        scatter3(d1,min(d2)*ones(1000,1),zeros(1000,1),10,'k')
        xlabel([])
        ylabel([])
        zlabel([])
        %title(sprintf('alpha:%i,beta:%i',alpha(k),beta(p)))
        saveas(h,sprintf('Plots/alpha_%i_beta_%i.png',alpha(k),beta(p)))
        
    end
end



