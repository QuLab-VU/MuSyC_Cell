%This code is to generate example dose response surfaces in matlab
%Set default plot parameters
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultLineLinewidth',4)
set(0, 'DefaultAxesFontWeight','bold')
set(0, 'DefaultAxesLineWidth',4)
set(0, 'DefaultFigureColor','w')
set(0, 'DefaultTextFontSize',18)
set(0, 'DefaultTextFontWeight','bold')

%Set parameters for surfaces
d1 = logspace(-5,0,100);
d2 = logspace(-5,0,100);
h1 = 1;
h2 = 1;
C1 = .005;
C2 = .001;
E0 = .05;
E1 = .01;
E2 = .03;
beta = 0.001;
alpha1 = 1;
alpha2 = 1;
gamma1 = 1;
gamma2 = 1;
r1 = 100*E0;
r2 = 100*E0;
Ed = zeros(length(d1),length(d2));

v = [1,logspace(0,1.5,50),10^(1.5)*ones(1,15),logspace(1.5,-1,100),1e-1*ones(1,15),logspace(-1,0,50)];
cnt = 1;
%Vary alpha
for alpha1 = v
    alpha2 = alpha1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Null surface plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = figure('color','w','pos',[10,10,1024,512]);
    figure(h)
    clf
    for i = 1:length(d1)
        for j = 1:length(d2)
            E3 = -beta*(E0-min(E1,E2))+min(E1,E2);
            Ed(i,j) = (C1^(2*gamma1*h1)*E2*d2(j)^h2*r1 + C2^(2*gamma2*h2)*E1*d1(i)^h1*r2 + C1^(2*gamma1*h1)*C2^(gamma2*h2)*E0*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*E0*r2 + E3*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1*(alpha1*d2(j))^(gamma1*h2) + E3*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*E2*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*E1*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r1*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*E2*d1(i)^h1*r1*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*E3*d2(j)^h2*r1*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*E1*d2(j)^h2*r2*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*E3*d1(i)^h1*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*C2^(gamma2*h2)*E1*d1(i)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E2*d2(j)^h2*r2)/(C1^(2*gamma1*h1)*C2^(gamma2*h2)*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*r2 + C1^(2*gamma1*h1)*d2(j)^h2*r1 + C2^(2*gamma2*h2)*d1(i)^h1*r2 + C1^(gamma1*h1)*C2^(gamma2*h2)*r1*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*d1(i)^h1*r1*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*d2(j)^h2*r1*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*d1(i)^h1*r2*(alpha1*d2(j))^(gamma1*h2) + C2^(gamma2*h2)*d2(j)^h2*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*d1(i)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*d2(j)^h2*r2 + alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1*(alpha1*d2(j))^(gamma1*h2) + alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1);
        end
    end
    %Downloaded from http://www.kennethmoreland.com/color-maps/
    map = diverging_map(0:.001:1,[.23,.299,.745],[.706,.016,.15]);
    [X,Y] = meshgrid(d1,d2);
    C = contour(X,Y,Ed,20);
    clf
    subtightplot(1,2,1,[0,0,0])
    hold on
    surf(X,Y,Ed,'lines','None')
    colormap(map)
    set(gca,'xscale','log','yscale','log','linewidth',2.,'xticklabels',[],'yticklabels',[],'zticklabels',[],'xtick',logspace(-5,-1,5),'ytick',logspace(-5,-1,5))
    grid off
    axis([min(d1),max(d1),min(d2),max(d2),-.03,max(Ed(:))])
    view([67,23])
    xlabel('log([Drug 1])')
    ylabel('log([Drug 2])')
    zlabel('DIP')
    hold on
    scatter3(d1,min(d2)*ones(100,1),Ed(1,:),30,[.7,.7,.7],'filled')
    scatter3(min(d1)*ones(100,1),d2,Ed(:,1),30,'k','filled')
    scatter3(min(d1),min(d2),max(Ed(:)),150,'r','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(max(d1),min(d2),Ed(1,end),150,'m','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(min(d1),max(d2),Ed(end,1),150,'m','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(max(d1),max(d2),Ed(end,end),150,'b','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min((d2-C2).^2);
    scatter3(C2,min(d2),Ed(1,m_idx),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min((d1-C1).^2);
    scatter3(min(d2),C1,Ed(m_idx,1),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    
    [~,m_idx] = min(((max(Ed(:,end))-max(Ed(:,end)-min(Ed(:,end)))/2)-Ed(:,end)).^2);
    scatter3(max(d1),d2(m_idx),Ed(m_idx,end),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min(((max(Ed(end,:))-max(Ed(end,:)-min(Ed(end,:)))/2)-Ed(end,:)).^2);
    scatter3(d1(m_idx),max(d2),Ed(end,m_idx),150,'c','filled','MarkerEdgeColor','k','linewidth',2)

    s = contourdata(C); 
    for i = 1:length(s)
        cont = zeros(length(s(i).xdata),1);
        for j = 1:length(s(i).xdata)
                cont(j) = (C1^(2*gamma1*h1)*E2*s(i).ydata(j)^h2*r1 + C2^(2*gamma2*h2)*E1*s(i).xdata(j)^h1*r2 + C1^(2*gamma1*h1)*C2^(gamma2*h2)*E0*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*E0*r2 + E3*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + E3*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*E2*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*E1*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*E2*s(i).xdata(j)^h1*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*E3*s(i).ydata(j)^h2*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*E1*s(i).ydata(j)^h2*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*E3*s(i).xdata(j)^h1*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*C2^(gamma2*h2)*E1*s(i).xdata(j)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E2*s(i).ydata(j)^h2*r2)/(C1^(2*gamma1*h1)*C2^(gamma2*h2)*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*r2 + C1^(2*gamma1*h1)*s(i).ydata(j)^h2*r1 + C2^(2*gamma2*h2)*s(i).xdata(j)^h1*r2 + C1^(gamma1*h1)*C2^(gamma2*h2)*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*s(i).xdata(j)^h1*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*s(i).ydata(j)^h2*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*s(i).xdata(j)^h1*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C2^(gamma2*h2)*s(i).ydata(j)^h2*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*s(i).xdata(j)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*s(i).ydata(j)^h2*r2 + alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1);
        end
        %plot3(s(i).xdata,s(i).ydata,cont,'color',map(i,:),'linewidth',1)
        plot3(s(i).ydata,s(i).xdata,cont,'k','linewidth',1)
    end
    set(gca,'ZGrid','on')
    scatter3(ones(100,1)*min(d1),d2,zeros(100,1),10,'k')
    scatter3(d1,max(d2)*ones(100,1),zeros(100,1),10,'k')
    scatter3(ones(100,1)*max(d1),d2,zeros(100,1),10,'k')
    scatter3(d1,min(d2)*ones(100,1),zeros(100,1),10,'k')
    xlabel([])
    ylabel([])
    zlabel([])
    if alpha1<1
        title('Antagonistic Potency','fontsize',12)
    else
        title('Synergistic Potency','fontsize',12)
    end
    
    subtightplot(1,2,2,[0,0,0])
    hold on
    plot([0,0],[-.5,1],'k')
    plot([-1,2],[0,0],'k')
    scatter(log10(alpha1),beta,200,'r','filled')
    axis off
    text(.6,.08,'<--log(\alpha)-->')
    ht = text(.15,.65,'<--\beta-->');
    text(-.7,-.1,'(0,0)')
    set(ht,'Rotation',-90)
    
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if cnt==1
      imwrite(imind,cm,'combo_surf_mov','gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,'combo_surf_mov','gif','WriteMode','append'); 
    end 
    cnt = cnt + 1;
    close(h)
end
      

alpha1 = 1;
alpha2 = 1;
v = [0.001,linspace(0.001,1,50),1*ones(1,15),linspace(1,-.5,100),-.5*ones(1,15),linspace(-.5,0.001,50)];
for beta = v
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Null surface plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = figure('color','w','pos',[10,10,1024,512]);
    figure(h)
    clf
    for i = 1:length(d1)
        for j = 1:length(d2)
            E3 = -beta*(E0-min(E1,E2))+min(E1,E2);
            Ed(i,j) = (C1^(2*gamma1*h1)*E2*d2(j)^h2*r1 + C2^(2*gamma2*h2)*E1*d1(i)^h1*r2 + C1^(2*gamma1*h1)*C2^(gamma2*h2)*E0*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*E0*r2 + E3*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1*(alpha1*d2(j))^(gamma1*h2) + E3*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*E2*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*E1*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r1*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*E2*d1(i)^h1*r1*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*E3*d2(j)^h2*r1*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*E1*d2(j)^h2*r2*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*E3*d1(i)^h1*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*C2^(gamma2*h2)*E1*d1(i)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E2*d2(j)^h2*r2)/(C1^(2*gamma1*h1)*C2^(gamma2*h2)*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*r2 + C1^(2*gamma1*h1)*d2(j)^h2*r1 + C2^(2*gamma2*h2)*d1(i)^h1*r2 + C1^(gamma1*h1)*C2^(gamma2*h2)*r1*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*d1(i)^h1*r1*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*d2(j)^h2*r1*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*d1(i)^h1*r2*(alpha1*d2(j))^(gamma1*h2) + C2^(gamma2*h2)*d2(j)^h2*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*d1(i)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*d2(j)^h2*r2 + alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1*(alpha1*d2(j))^(gamma1*h2) + alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1); 
        end
    end
    %Downloaded from http://www.kennethmoreland.com/color-maps/
    map = diverging_map(0:.001:1,[.23,.299,.745],[.706,.016,.15]);
    [X,Y] = meshgrid(d1,d2);
    C = contour(X,Y,Ed,20);
    clf
    subtightplot(1,2,1,[0,0,0])
    hold on
    surf(X,Y,Ed,'lines','None')
    colormap(map)
    set(gca,'xscale','log','yscale','log','linewidth',2.,'xticklabels',[],'yticklabels',[],'zticklabels',[],'xtick',logspace(-5,-1,5),'ytick',logspace(-5,-1,5))
    grid off
    axis([min(d1),max(d1),min(d2),max(d2),-.03,max(Ed(:))])
    view([67,23])
    xlabel('log([Drug 1])')
    ylabel('log([Drug 2])')
    zlabel('DIP')
    hold on
    scatter3(d1,min(d2)*ones(100,1),Ed(1,:),30,[.7,.7,.7],'filled')
    scatter3(min(d1)*ones(100,1),d2,Ed(:,1),30,'k','filled')
    scatter3(min(d1),min(d2),max(Ed(:)),150,'r','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(max(d1),min(d2),Ed(1,end),150,'m','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(min(d1),max(d2),Ed(end,1),150,'m','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(max(d1),max(d2),Ed(end,end),150,'b','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min((d2-C2).^2);
    scatter3(C2,min(d2),Ed(1,m_idx),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min((d1-C1).^2);
    scatter3(min(d2),C1,Ed(m_idx,1),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    
    [~,m_idx] = min(((max(Ed(:,end))-max(Ed(:,end)-min(Ed(:,end)))/2)-Ed(:,end)).^2);
    scatter3(max(d1),d2(m_idx),Ed(m_idx,end),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min(((max(Ed(end,:))-max(Ed(end,:)-min(Ed(end,:)))/2)-Ed(end,:)).^2);
    scatter3(d1(m_idx),max(d2),Ed(end,m_idx),150,'c','filled','MarkerEdgeColor','k','linewidth',2)

    s = contourdata(C); 
    for i = 1:length(s)
        cont = zeros(length(s(i).xdata),1);
        for j = 1:length(s(i).xdata)
                cont(j) = (C1^(2*gamma1*h1)*E2*s(i).ydata(j)^h2*r1 + C2^(2*gamma2*h2)*E1*s(i).xdata(j)^h1*r2 + C1^(2*gamma1*h1)*C2^(gamma2*h2)*E0*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*E0*r2 + E3*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + E3*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*E2*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*E1*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*E2*s(i).xdata(j)^h1*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*E3*s(i).ydata(j)^h2*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*E1*s(i).ydata(j)^h2*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*E3*s(i).xdata(j)^h1*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*C2^(gamma2*h2)*E1*s(i).xdata(j)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E2*s(i).ydata(j)^h2*r2)/(C1^(2*gamma1*h1)*C2^(gamma2*h2)*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*r2 + C1^(2*gamma1*h1)*s(i).ydata(j)^h2*r1 + C2^(2*gamma2*h2)*s(i).xdata(j)^h1*r2 + C1^(gamma1*h1)*C2^(gamma2*h2)*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*s(i).xdata(j)^h1*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*s(i).ydata(j)^h2*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*s(i).xdata(j)^h1*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C2^(gamma2*h2)*s(i).ydata(j)^h2*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*s(i).xdata(j)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*s(i).ydata(j)^h2*r2 + alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1); 
        end
        %plot3(s(i).xdata,s(i).ydata,cont,'color',map(i,:),'linewidth',1)
        plot3(s(i).ydata,s(i).xdata,cont,'k','linewidth',1)
    end
    set(gca,'ZGrid','on')
    scatter3(ones(100,1)*min(d1),d2,zeros(100,1),10,'k')
    scatter3(d1,max(d2)*ones(100,1),zeros(100,1),10,'k')
    scatter3(ones(100,1)*max(d1),d2,zeros(100,1),10,'k')
    scatter3(d1,min(d2)*ones(100,1),zeros(100,1),10,'k')
    xlabel([])
    ylabel([])
    zlabel([])
    if beta<0
        title('Antagonistic Efficacy','fontsize',12)
    else
        title('Synergistic Efficacy','fontsize',12)
    end
    
    subtightplot(1,2,2,[0,0,0])
    hold on
    plot([0,0],[-.5,1],'k')
    plot([-1,2],[0,0],'k')
    scatter(log10(alpha1),beta,200,'r','filled')
    axis off
    text(.6,.08,'<--log(\alpha)-->')
    ht = text(.15,.65,'<--\beta-->');
    text(-.7,-.1,'(0,0)')
    set(ht,'Rotation',-90)
    
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if cnt==1
      imwrite(imind,cm,'combo_surf_mov','gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,'combo_surf_mov','gif','WriteMode','append'); 
    end 
    cnt = cnt + 1;
    close(h)
end

va = [1,logspace(0,1.5,50),10^(1.5)*ones(1,15),logspace(1.5,-1,100),1e-1*ones(1,15),logspace(-1,0,50)];
vb = [0.001,linspace(0.001,1,50),1*ones(1,15),linspace(1,-.5,100),-.5*ones(1,15),linspace(-.5,0.001,50)];
for vv = 1:length(v)
    beta  = vb(vv);
    alpha1 = va(vv);
    alpha2 = va(vv);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Null surface plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = figure('color','w','pos',[10,10,1024,512]);
    figure(h)
    clf
    for i = 1:length(d1)
        for j = 1:length(d2)
            E3 = -beta*(E0-min(E1,E2))+min(E1,E2);
            Ed(i,j) = (C1^(2*gamma1*h1)*E2*d2(j)^h2*r1 + C2^(2*gamma2*h2)*E1*d1(i)^h1*r2 + C1^(2*gamma1*h1)*C2^(gamma2*h2)*E0*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*E0*r2 + E3*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1*(alpha1*d2(j))^(gamma1*h2) + E3*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*E2*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*E1*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r1*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*E2*d1(i)^h1*r1*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*E3*d2(j)^h2*r1*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*E1*d2(j)^h2*r2*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*E3*d1(i)^h1*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*C2^(gamma2*h2)*E1*d1(i)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E2*d2(j)^h2*r2)/(C1^(2*gamma1*h1)*C2^(gamma2*h2)*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*r2 + C1^(2*gamma1*h1)*d2(j)^h2*r1 + C2^(2*gamma2*h2)*d1(i)^h1*r2 + C1^(gamma1*h1)*C2^(gamma2*h2)*r1*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*d1(i)^h1*r1*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*d2(j)^h2*r1*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*d1(i)^h1*r2*(alpha1*d2(j))^(gamma1*h2) + C2^(gamma2*h2)*d2(j)^h2*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*d1(i)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*d2(j)^h2*r2 + alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1*(alpha1*d2(j))^(gamma1*h2) + alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1); 
        end
    end
    %Downloaded from http://www.kennethmoreland.com/color-maps/
    map = diverging_map(0:.001:1,[.23,.299,.745],[.706,.016,.15]);
    [X,Y] = meshgrid(d1,d2);
    C = contour(X,Y,Ed,20);
    clf
    subtightplot(1,2,1,[0,0,0])
    hold on
    surf(X,Y,Ed,'lines','None')
    colormap(map)
    set(gca,'xscale','log','yscale','log','linewidth',2.,'xticklabels',[],'yticklabels',[],'zticklabels',[],'xtick',logspace(-5,-1,5),'ytick',logspace(-5,-1,5))
    grid off
    axis([min(d1),max(d1),min(d2),max(d2),-.03,max(Ed(:))])
    view([67,23])
    xlabel('log([Drug 1])')
    ylabel('log([Drug 2])')
    zlabel('DIP')
    hold on
    scatter3(d1,min(d2)*ones(100,1),Ed(1,:),30,[.7,.7,.7],'filled')
    scatter3(min(d1)*ones(100,1),d2,Ed(:,1),30,'k','filled')
    scatter3(min(d1),min(d2),max(Ed(:)),150,'r','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(max(d1),min(d2),Ed(1,end),150,'m','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(min(d1),max(d2),Ed(end,1),150,'m','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(max(d1),max(d2),Ed(end,end),150,'b','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min((d2-C2).^2);
    scatter3(C2,min(d2),Ed(1,m_idx),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min((d1-C1).^2);
    scatter3(min(d2),C1,Ed(m_idx,1),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    
    [~,m_idx] = min(((max(Ed(:,end))-max(Ed(:,end)-min(Ed(:,end)))/2)-Ed(:,end)).^2);
    scatter3(max(d1),d2(m_idx),Ed(m_idx,end),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min(((max(Ed(end,:))-max(Ed(end,:)-min(Ed(end,:)))/2)-Ed(end,:)).^2);
    scatter3(d1(m_idx),max(d2),Ed(end,m_idx),150,'c','filled','MarkerEdgeColor','k','linewidth',2)

    s = contourdata(C); 
    for i = 1:length(s)
        cont = zeros(length(s(i).xdata),1);
        for j = 1:length(s(i).xdata)
                cont(j) = (C1^(2*gamma1*h1)*E2*s(i).ydata(j)^h2*r1 + C2^(2*gamma2*h2)*E1*s(i).xdata(j)^h1*r2 + C1^(2*gamma1*h1)*C2^(gamma2*h2)*E0*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*E0*r2 + E3*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + E3*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*E2*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*E1*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*E2*s(i).xdata(j)^h1*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*E3*s(i).ydata(j)^h2*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*E1*s(i).ydata(j)^h2*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*E3*s(i).xdata(j)^h1*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*C2^(gamma2*h2)*E1*s(i).xdata(j)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E2*s(i).ydata(j)^h2*r2)/(C1^(2*gamma1*h1)*C2^(gamma2*h2)*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*r2 + C1^(2*gamma1*h1)*s(i).ydata(j)^h2*r1 + C2^(2*gamma2*h2)*s(i).xdata(j)^h1*r2 + C1^(gamma1*h1)*C2^(gamma2*h2)*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*s(i).xdata(j)^h1*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*s(i).ydata(j)^h2*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*s(i).xdata(j)^h1*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C2^(gamma2*h2)*s(i).ydata(j)^h2*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*s(i).xdata(j)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*s(i).ydata(j)^h2*r2 + alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1); 
        end
        %plot3(s(i).xdata,s(i).ydata,cont,'color',map(i,:),'linewidth',1)
        plot3(s(i).ydata,s(i).xdata,cont,'k','linewidth',1)
    end
    set(gca,'ZGrid','on')
    scatter3(ones(100,1)*min(d1),d2,zeros(100,1),10,'k')
    scatter3(d1,max(d2)*ones(100,1),zeros(100,1),10,'k')
    scatter3(ones(100,1)*max(d1),d2,zeros(100,1),10,'k')
    scatter3(d1,min(d2)*ones(100,1),zeros(100,1),10,'k')
    xlabel([])
    ylabel([])
    zlabel([])
    if beta<=0 && alpha<=1
        title('Antagonistic Efficacy, Antagonistic Potency','fontsize',12)
    elseif beta>=0 && alpha<=1
        title('Synergistic Efficacy, Antagonistic Potency','fontsize',12)
    elseif beta<=0 && alpha>=1
        title('Antagonistic Efficacy, Synergistic Potency','fontsize',12)
    elseif beta>=0 && alpha>=1
        title('Synergistic Efficacy, Synergistic Potency','fontsize',12)
    end
    
    subtightplot(1,2,2,[0,0,0])
    hold on
    plot([0,0],[-.5,1],'k')
    plot([-1,2],[0,0],'k')
    scatter(log10(alpha1),beta,200,'r','filled')
    axis off
    text(.6,.08,'<--log(\alpha)-->')
    ht = text(.15,.65,'<--\beta-->');
    text(-.7,-.1,'(0,0)')
    set(ht,'Rotation',-90)
    
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if cnt==1
      imwrite(imind,cm,'combo_surf_mov','gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,'combo_surf_mov','gif','WriteMode','append'); 
    end 
    cnt = cnt + 1;
    close(h)
end


va = [1,logspace(0,1.5,50),10^(1.5)*ones(1,15),logspace(1.5,-1,100),1e-1*ones(1,15),logspace(-1,0,50)];
vb = [0.001,linspace(0.001,1,50),1*ones(1,15),linspace(1,-.5,100),-.5*ones(1,15),linspace(-.5,0.001,50)];
for vv = 1:length(v)
    beta  = vb(end-vv+1);
    alpha1 = va(vv);
    alpha2 = va(vv);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Null surface plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = figure('color','w','pos',[10,10,1024,512]);
    figure(h)
    clf
    for i = 1:length(d1)
        for j = 1:length(d2)
            E3 = -beta*(E0-min(E1,E2))+min(E1,E2);
            Ed(i,j) = (C1^(2*gamma1*h1)*E2*d2(j)^h2*r1 + C2^(2*gamma2*h2)*E1*d1(i)^h1*r2 + C1^(2*gamma1*h1)*C2^(gamma2*h2)*E0*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*E0*r2 + E3*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1*(alpha1*d2(j))^(gamma1*h2) + E3*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*E2*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*E1*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r1*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*E2*d1(i)^h1*r1*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*E3*d2(j)^h2*r1*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*E1*d2(j)^h2*r2*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*E3*d1(i)^h1*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*C2^(gamma2*h2)*E1*d1(i)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E2*d2(j)^h2*r2)/(C1^(2*gamma1*h1)*C2^(gamma2*h2)*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*r2 + C1^(2*gamma1*h1)*d2(j)^h2*r1 + C2^(2*gamma2*h2)*d1(i)^h1*r2 + C1^(gamma1*h1)*C2^(gamma2*h2)*r1*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*r2*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*d1(i)^h1*r1*(alpha1*d2(j))^(gamma1*h2) + C1^(gamma1*h1)*d2(j)^h2*r1*(alpha2*d1(i))^(gamma2*h1) + C2^(gamma2*h2)*d1(i)^h1*r2*(alpha1*d2(j))^(gamma1*h2) + C2^(gamma2*h2)*d2(j)^h2*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*d1(i)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*d2(j)^h2*r2 + alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1*(alpha1*d2(j))^(gamma1*h2) + alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2*(alpha2*d1(i))^(gamma2*h1) + C1^(gamma1*h1)*alpha1^(gamma1*h2)*d2(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*alpha2^(gamma2*h1)*d1(i)^(h1*(gamma2 + 1))*r1); 
        end
    end
    %Downloaded from http://www.kennethmoreland.com/color-maps/
    map = diverging_map(0:.001:1,[.23,.299,.745],[.706,.016,.15]);
    [X,Y] = meshgrid(d1,d2);
    C = contour(X,Y,Ed,20);
    clf
    subtightplot(1,2,1,[0,0,0])
    hold on
    surf(X,Y,Ed,'lines','None')
    colormap(map)
    set(gca,'xscale','log','yscale','log','linewidth',2.,'xticklabels',[],'yticklabels',[],'zticklabels',[],'xtick',logspace(-5,-1,5),'ytick',logspace(-5,-1,5))
    grid off
    axis([min(d1),max(d1),min(d2),max(d2),-.03,max(Ed(:))])
    view([67,23])
    xlabel('log([Drug 1])')
    ylabel('log([Drug 2])')
    zlabel('DIP')
    hold on
    scatter3(d1,min(d2)*ones(100,1),Ed(1,:),30,[.7,.7,.7],'filled')
    scatter3(min(d1)*ones(100,1),d2,Ed(:,1),30,'k','filled')
    scatter3(min(d1),min(d2),max(Ed(:)),150,'r','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(max(d1),min(d2),Ed(1,end),150,'m','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(min(d1),max(d2),Ed(end,1),150,'m','filled','MarkerEdgeColor','k','linewidth',2)
    scatter3(max(d1),max(d2),Ed(end,end),150,'b','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min((d2-C2).^2);
    scatter3(C2,min(d2),Ed(1,m_idx),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min((d1-C1).^2);
    scatter3(min(d2),C1,Ed(m_idx,1),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    
    [~,m_idx] = min(((max(Ed(:,end))-max(Ed(:,end)-min(Ed(:,end)))/2)-Ed(:,end)).^2);
    scatter3(max(d1),d2(m_idx),Ed(m_idx,end),150,'c','filled','MarkerEdgeColor','k','linewidth',2)
    [~,m_idx] = min(((max(Ed(end,:))-max(Ed(end,:)-min(Ed(end,:)))/2)-Ed(end,:)).^2);
    scatter3(d1(m_idx),max(d2),Ed(end,m_idx),150,'c','filled','MarkerEdgeColor','k','linewidth',2)

    s = contourdata(C); 
    for i = 1:length(s)
        cont = zeros(length(s(i).xdata),1);
        for j = 1:length(s(i).xdata)
                cont(j) = (C1^(2*gamma1*h1)*E2*s(i).ydata(j)^h2*r1 + C2^(2*gamma2*h2)*E1*s(i).xdata(j)^h1*r2 + C1^(2*gamma1*h1)*C2^(gamma2*h2)*E0*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*E0*r2 + E3*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + E3*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*E2*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*E1*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*E0*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*E2*s(i).xdata(j)^h1*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*E3*s(i).ydata(j)^h2*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*E1*s(i).ydata(j)^h2*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*E3*s(i).xdata(j)^h1*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*C2^(gamma2*h2)*E1*s(i).xdata(j)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*E2*s(i).ydata(j)^h2*r2)/(C1^(2*gamma1*h1)*C2^(gamma2*h2)*r1 + C1^(gamma1*h1)*C2^(2*gamma2*h2)*r2 + C1^(2*gamma1*h1)*s(i).ydata(j)^h2*r1 + C2^(2*gamma2*h2)*s(i).xdata(j)^h1*r2 + C1^(gamma1*h1)*C2^(gamma2*h2)*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*s(i).xdata(j)^h1*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + C1^(gamma1*h1)*s(i).ydata(j)^h2*r1*(alpha2*s(i).xdata(j))^(gamma2*h1) + C2^(gamma2*h2)*s(i).xdata(j)^h1*r2*(alpha1*s(i).ydata(j))^(gamma1*h2) + C2^(gamma2*h2)*s(i).ydata(j)^h2*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*C2^(gamma2*h2)*s(i).xdata(j)^h1*r1 + C1^(gamma1*h1)*C2^(gamma2*h2)*s(i).ydata(j)^h2*r2 + alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1*(alpha1*s(i).ydata(j))^(gamma1*h2) + alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2*(alpha2*s(i).xdata(j))^(gamma2*h1) + C1^(gamma1*h1)*alpha1^(gamma1*h2)*s(i).ydata(j)^(h2*(gamma1 + 1))*r2 + C2^(gamma2*h2)*alpha2^(gamma2*h1)*s(i).xdata(j)^(h1*(gamma2 + 1))*r1); 
        end
        %plot3(s(i).xdata,s(i).ydata,cont,'color',map(i,:),'linewidth',1)
        plot3(s(i).ydata,s(i).xdata,cont,'k','linewidth',1)
    end
    set(gca,'ZGrid','on')
    scatter3(ones(100,1)*min(d1),d2,zeros(100,1),10,'k')
    scatter3(d1,max(d2)*ones(100,1),zeros(100,1),10,'k')
    scatter3(ones(100,1)*max(d1),d2,zeros(100,1),10,'k')
    scatter3(d1,min(d2)*ones(100,1),zeros(100,1),10,'k')
    xlabel([])
    ylabel([])
    zlabel([])
    if beta<=0 && alpha<=1
        title('Antagonistic Efficacy, Antagonistic Potency','fontsize',12)
    elseif beta>=0 && alpha<=1
        title('Synergistic Efficacy, Antagonistic Potency','fontsize',12)
    elseif beta<=0 && alpha>=1
        title('Antagonistic Efficacy, Synergistic Potency','fontsize',12)
    elseif beta>=0 && alpha>=1
        title('Synergistic Efficacy, Synergistic Potency','fontsize',12)
    end
    
    subtightplot(1,2,2,[0,0,0])
    hold on
    plot([0,0],[-.5,1],'k')
    plot([-1,2],[0,0],'k')
    scatter(log10(alpha1),beta,200,'r','filled')
    axis off
    text(.6,.08,'<--log(\alpha)-->')
    ht = text(.15,.65,'<--\beta-->');
    text(-.7,-.1,'(0,0)')
    set(ht,'Rotation',-90)
    
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if cnt==1
      imwrite(imind,cm,'combo_surf_mov','gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,'combo_surf_mov','gif','WriteMode','append'); 
    end 
    cnt = cnt + 1;
    close(h)
end



