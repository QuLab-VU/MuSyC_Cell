function varargout = synGUI(varargin)
% SYNGUI MATLAB code for synGUI.fig
%      SYNGUI, by itself, creates a new SYNGUI or raises the existing
%      singleton*.
%
%      H = SYNGUI returns the handle to a new SYNGUI or the handle to
%      the existing singleton*.
%
%      SYNGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYNGUI.M with the given input arguments.
%
%      SYNGUI('Property','Value',...) creates a new SYNGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before synGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to synGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help synGUI

% Last Modified by GUIDE v2.5 25-Oct-2017 16:19:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @synGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @synGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before synGUI is made visible.
function synGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to synGUI (see VARARGIN)

% Choose default command line output for synGUI
handles.axes1;
handles.axes2;
handles.axes3;
handles.output = hObject;


handles = InitializeHandles(handles);

%Set default plot parameters
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultLineLinewidth',3)
set(0, 'DefaultAxesFontWeight','bold')
set(0, 'DefaultAxesLineWidth',3)
set(0, 'DefaultFigureColor','w')
set(0, 'DefaultTextFontSize',12)
set(0, 'DefaultTextFontWeight','bold')

%Show model
handles = GetDefaults(handles);
handles = modelPic(handles);
handles = InitializeHandles(handles);
handles = PlotFun(handles);

%save('handles.mat','handles')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes synGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function [handles] = GetDefaults(handles)
    h1 = get(handles.slider6,'Value');
    h2 = get(handles.slider7,'Value');
    C1 = 10^get(handles.slider4,'Value');
    C2 = 10^get(handles.slider5,'Value');
    E0 = get(handles.slider1,'Value');
    E1 = get(handles.slider2,'Value');
    E2 = get(handles.slider3,'Value');
    alpha1 = get(handles.slider9,'Value');
    alpha2 = get(handles.slider10,'Value');
    beta = get(handles.slider11,'Value');
    trans_ratio = get(handles.slider12,'Value');
    scale_rate = get(handles.slider13,'Value');
    handles.defaults = [h1,h2,C1,C2,E0,E1,E2,alpha1,alpha2,beta,trans_ratio,scale_rate];

function [handles] = InitializeHandles(handles)
    h1 = get(handles.slider6,'Value');
    h2 = get(handles.slider7,'Value');
    C1 = 10^get(handles.slider4,'Value');
    C2 = 10^get(handles.slider5,'Value');
    E0 = get(handles.slider1,'Max')+.005;
    E1 = get(handles.slider2,'Min');
    E2 = get(handles.slider3,'Min');
    alpha1 = get(handles.slider9,'Value');
    alpha2 = get(handles.slider10,'Value');
    beta = get(handles.slider11,'Max');
    trans_ratio = get(handles.slider12,'Value');
    scale_rate = get(handles.slider13,'Value');    
    E3 = (-beta*(E0-min(E1,E2))+min(E1,E2))-.005;
    d1 = logspace(-5,0,1000);
    d2 = logspace(-5,0,1000);
    r1 = scale_rate*E0;
    r2 = r1/10^trans_ratio;
    mx_Ed = -1;
    m_Ed = 1;
    %Find surface
    for i = 1:length(d1)
        for j = 1:length(d2)
            Ed(i,j) = (E3*r1*(alpha2*d2(j))^h2*(alpha1*d1(i)^2)^h1 + E3*r2*(alpha1*d1(i))^h1*(alpha2*d2(j)^2)^h2 + C2^h2*E0*r1*(C1*alpha1*d1(i))^h1 + C1^h1*E0*r2*(C2*alpha2*d2(j))^h2 + C2^h2*E1*r1*(alpha1*d1(i)^2)^h1 + C1^h1*E2*r2*(alpha2*d2(j)^2)^h2 + E3*d2(j)^h2*r1*(C1*alpha1*d1(i))^h1 + E3*d1(i)^h1*r2*(C2*alpha2*d2(j))^h2 + C1^(2*h1)*C2^h2*E0*r1 + C1^h1*C2^(2*h2)*E0*r2 + C1^(2*h1)*E2*d2(j)^h2*r1 + C2^(2*h2)*E1*d1(i)^h1*r2 + E1*r2*(C2*d2(j))^h2*(alpha1*d1(i))^h1 + E2*r1*(C1*d1(i))^h1*(alpha2*d2(j))^h2 + C2^h2*E1*r1*(C1*d1(i))^h1 + C1^h1*E2*r2*(C2*d2(j))^h2)/(r1*(alpha2*d2(j))^h2*(alpha1*d1(i)^2)^h1 + r2*(alpha1*d1(i))^h1*(alpha2*d2(j)^2)^h2 + C2^h2*r1*(C1*alpha1*d1(i))^h1 + C1^h1*r2*(C2*alpha2*d2(j))^h2 + C2^h2*r1*(alpha1*d1(i)^2)^h1 + C1^h1*r2*(alpha2*d2(j)^2)^h2 + d2(j)^h2*r1*(C1*alpha1*d1(i))^h1 + d1(i)^h1*r2*(C2*alpha2*d2(j))^h2 + C1^(2*h1)*C2^h2*r1 + C1^h1*C2^(2*h2)*r2 + C1^(2*h1)*d2(j)^h2*r1 + C2^(2*h2)*d1(i)^h1*r2 + r1*(C1*d1(i))^h1*(alpha2*d2(j))^h2 + r2*(C2*d2(j))^h2*(alpha1*d1(i))^h1 + C2^h2*r1*(C1*d1(i))^h1 + C1^h1*r2*(C2*d2(j))^h2);
        end
    end
    handles.m_Ed = min(Ed(:));
    handles.mx_Ed = max(Ed(:));


function [handles] = modelPic(handles)
    axes(handles.axes4)
    cla
    hold on
    E0 = get(handles.slider1,'Value');
    E1 = get(handles.slider2,'Value');
    E2 = get(handles.slider3,'Value');
    beta = get(handles.slider11,'Value');
    E3 = (-beta*(E0-min(E1,E2))+min(E1,E2));
    plot([.5,.5],[.5,3.5],'k')
    plot([.5,.5],[6.5,9.5],'k')
    plot([0,2.5],[2,2],'k')
    plot([0,2.5],[8,8],'k')
    plot([7.5,7.5],[.5,3.5],'k')
    plot([7.5,7.5],[6.5,9.5],'k')
    plot([7,9.5],[2,2],'k')
    plot([7,9.5],[8,8],'k')
    p1 = [1.25 3.5];        
    p2 = [1.25 6.5];                         
    dp = p2-p1;                       
    quiver(p1(1),p1(2),dp(1),dp(2),0,'k')
    p2 = [1.25 3.5];        
    p1 = [1.25 6.5];                         
    dp = p2-p1;                       
    quiver(p1(1),p1(2),dp(1),dp(2),0,'k')
    p1 = [8.25 3.5];        
    p2 = [8.25 6.5];                         
    dp = p2-p1;                       
    quiver(p1(1),p1(2),dp(1),dp(2),0,'k')
    p2 = [8.25 3.5];        
    p1 = [8.25 6.5];                         
    dp = p2-p1;                       
    quiver(p1(1),p1(2),dp(1),dp(2),0,'k')
    p1 = [3 2];        
    p2 = [6.5 2];                         
    dp = p2-p1;                       
    quiver(p1(1),p1(2),dp(1),dp(2),0,'k')
    p2 = [3 2];        
    p1 = [6.5 2];                         
    dp = p2-p1;                       
    quiver(p1(1),p1(2),dp(1),dp(2),0,'k')
    p1 = [3 7];        
    p2 = [6.5 7];                         
    dp = p2-p1;                       
    quiver(p1(1),p1(2),dp(1),dp(2),0,'k')
    p2 = [3 7];        
    p1 = [6.5 7];                         
    dp = p2-p1;                       
    quiver(p1(1),p1(2),dp(1),dp(2),0,'k')
    plot([.5,2.5],[2,2+15*E1*2],'r')
    plot([.5,2.5],[8,8+15*E0*2],'b')
    plot([7.5,9.5],[2,2+15*E3*2],'m')
    plot([7.5,9.5],[8,8+15*E2*2],'k')
    text(1.25,9,'E0')
    text(8.25,9,'E2')
    text(1.25,3,'E1')
    text(8.25,3,'E3')
    xlim([0,10])
    ylim([0,10])
    set(gca,'xcolor','none','ycolor','none')
        
function [handles] = PlotFun(handles)
    fid = fopen('proverbs.csv','r');
    for i=1:randi(525,1)
        str = fgets(fid);
    end
    popbox = msgbox('Rendering Graph','Hold please...');
    h1 = get(handles.slider6,'Value');
    h2 = get(handles.slider7,'Value');
    C1 = 10^get(handles.slider4,'Value');
    C2 = 10^get(handles.slider5,'Value');
    E0 = get(handles.slider1,'Value');
    E1 = get(handles.slider2,'Value');
    E2 = get(handles.slider3,'Value');
    alpha1 = get(handles.slider9,'Value');
    alpha2 = get(handles.slider10,'Value');
    beta = get(handles.slider11,'Value');
    trans_ratio = get(handles.slider12,'Value');
    scale_rate = get(handles.slider13,'Value'); 
    
    E3 = (-beta*(E0-min(E1,E2))+min(E1,E2));

    d1 = logspace(-5,0,1000);
    d2 = logspace(-5,0,1000);
    [X,Y] = meshgrid(d1,d2);

    r1 = scale_rate*E0;
    r2 = r1/10^trans_ratio;
    
    %Find surface
    for i = 1:length(d1)
        for j = 1:length(d2)
            Ed(i,j) = (E3*r1*(alpha2*d2(j))^h2*(alpha1*d1(i)^2)^h1 + E3*r2*(alpha1*d1(i))^h1*(alpha2*d2(j)^2)^h2 + C2^h2*E0*r1*(C1*alpha1*d1(i))^h1 + C1^h1*E0*r2*(C2*alpha2*d2(j))^h2 + C2^h2*E1*r1*(alpha1*d1(i)^2)^h1 + C1^h1*E2*r2*(alpha2*d2(j)^2)^h2 + E3*d2(j)^h2*r1*(C1*alpha1*d1(i))^h1 + E3*d1(i)^h1*r2*(C2*alpha2*d2(j))^h2 + C1^(2*h1)*C2^h2*E0*r1 + C1^h1*C2^(2*h2)*E0*r2 + C1^(2*h1)*E2*d2(j)^h2*r1 + C2^(2*h2)*E1*d1(i)^h1*r2 + E1*r2*(C2*d2(j))^h2*(alpha1*d1(i))^h1 + E2*r1*(C1*d1(i))^h1*(alpha2*d2(j))^h2 + C2^h2*E1*r1*(C1*d1(i))^h1 + C1^h1*E2*r2*(C2*d2(j))^h2)/(r1*(alpha2*d2(j))^h2*(alpha1*d1(i)^2)^h1 + r2*(alpha1*d1(i))^h1*(alpha2*d2(j)^2)^h2 + C2^h2*r1*(C1*alpha1*d1(i))^h1 + C1^h1*r2*(C2*alpha2*d2(j))^h2 + C2^h2*r1*(alpha1*d1(i)^2)^h1 + C1^h1*r2*(alpha2*d2(j)^2)^h2 + d2(j)^h2*r1*(C1*alpha1*d1(i))^h1 + d1(i)^h1*r2*(C2*alpha2*d2(j))^h2 + C1^(2*h1)*C2^h2*r1 + C1^h1*C2^(2*h2)*r2 + C1^(2*h1)*d2(j)^h2*r1 + C2^(2*h2)*d1(i)^h1*r2 + r1*(C1*d1(i))^h1*(alpha2*d2(j))^h2 + r2*(C2*d2(j))^h2*(alpha1*d1(i))^h1 + C2^h2*r1*(C1*d1(i))^h1 + C1^h1*r2*(C2*d2(j))^h2);
        end
    end
    
    %Single drug plots
    axes(handles.axes2)
    cla
    hold on
    set(gca,'xscale','log','linewidth',2.,'xticklabels',[],'yticklabels',[0],'xtick',logspace(-5,-0,6),'ytick',[0])
    plot(d1,Ed(1,:),'color',[.7,.7,.7])
    plot(d2,Ed(:,1),'color','k')
    legend('Drug 2','Drug 1','northeast')
    
    scatter(min(d1),max(Ed(:)),150,'r','filled')
    scatter(max(d1),Ed(1,end),150,'r','filled')
    scatter(max(d2),Ed(end,1),150,'r','filled')
    [~,m_idx] = min((d2-C2).^2);
    scatter(C2,Ed(1,m_idx),150,'r','filled')
    text(C2,Ed(1,m_idx)-.004,'C2')
    [~,m_idx] = min((d1-C1).^2);
    scatter(C1,Ed(m_idx,1),150,'r','filled')
    text(C1,Ed(m_idx,1)-.004,'C1')
    set(gca,'YGrid','on')
    scatter(d1,zeros(1000,1),10,'k','filled')
    text(min(d1),E0+.002,'E0')
    text(max(d1),E1+.004,'E1')
    text(max(d1),E2+.004,'E2')
    ylabel('DIP Rate')
    xlabel('log(Drug)')
    
    %alpha beta plot
    axes(handles.axes3)
    cla
    hold on
    scatter(min(log10([alpha1,alpha2])),beta,150,'r','filled')
    xlim([-1,1])
    ylim([-.25,.5])
    plot([0,0],[-.25,1.25],'k')
    plot([-1,1],[0,0],'k')
    text(.005,.002,'(0,0)')
    xlabel('min(log(alpha1,alpha2))')
    ylabel('beta')
%     str = sprintf('Antagonistic Potency \n <-------------------');
%     text(-.5,.025,str)
%     str = sprintf('Synergistic Potency \n ------------------->');
%     text(.5,.025,str)
%     
%     str = sprintf('Antagonistic Efficacy \n <-------------------');
%     h = text(.5,-.005,str);
%     set(h,'Rotation',90);
%     str = sprintf('Synergistic Efficacy \n ------------------->');
%     h = text(.5,.005,str);
%     set(h,'Rotation',90);

    axes(handles.axes1)
    C = contour(X,Y,Ed,20);
    cla
    s = contourdata(C);        %subplot(length(alpha),length(beta),(k-1)*length(alpha)+p)
    [X,Y] = meshgrid(d1,d2);
    surf(X,Y,Ed,'lines','None')
    set(gca,'xscale','log','yscale','log','linewidth',2.,'xticklabels',[],'yticklabels',[],'zticklabels',[],'xtick',logspace(-5,-1,5),'ytick',logspace(-5,-1,5))
    %Downloaded from http://www.kennethmoreland.com/color-maps/
    map = diverging_map(0:.001:1,[.23,.299,.745],[.706,.016,.15]);
    colormap(map)
    caxis([handles.m_Ed,handles.mx_Ed])
    grid off
    axis([min(d1),max(d1),min(d2),max(d2),handles.m_Ed,handles.mx_Ed])
    view([67,46])
    hold on
    map = diverging_map(0:1/length(s):1,[.23,.299,.745],[.706,.016,.15]);
    scatter3(d1,min(d2)*ones(1000,1),Ed(1,:),30,[.7,.7,.7],'filled')
    scatter3(min(d1)*ones(1000,1),d2,Ed(:,1),30,'k','filled')
    scatter3(min(d1),min(d2),Ed(1,1),150,'r','filled')
    scatter3(max(d1),min(d2),Ed(1,end),150,'m','filled')
    scatter3(min(d1),max(d2),Ed(end,1),150,'m','filled')
    scatter3(max(d1),max(d2),Ed(end,end),150,'b','filled')
    [~,m_idx] = min((d2-C2).^2);
    scatter3(C2,min(d2),Ed(1,m_idx),150,'c','filled')
    [~,m_idx] = min((d1-C1).^2);
    scatter3(min(d2),C1,Ed(m_idx,1),150,'c','filled')
    
    [~,m_idx] = min((Ed(end,:)-(max(Ed(end,:))-abs(Ed(end,1)-Ed(end,end))/2)).^2);
    scatter3(d1(m_idx),max(d2),Ed(end,m_idx),150,'c','filled')    
    [~,m_idx] = min((Ed(:,end)-(max(Ed(:,end))-abs(Ed(1,end)-Ed(end,end))/2)).^2);
    scatter3(max(d1),d2(m_idx),Ed(m_idx,end),150,'c','filled')

    
    
    s = contourdata(C); 
    for i = 1:length(s)
        cont = [];
        for j = 1:length(s(i).xdata)
                cont(j) = (E3*r1*(alpha2*s(i).xdata(j))^h2*(alpha1*s(i).ydata(j)^2)^h1 + E3*r2*(alpha1*s(i).ydata(j))^h1*(alpha2*s(i).xdata(j)^2)^h2 + C2^h2*E0*r1*(C1*alpha1*s(i).ydata(j))^h1 + C1^h1*E0*r2*(C2*alpha2*s(i).xdata(j))^h2 + C2^h2*E1*r1*(alpha1*s(i).ydata(j)^2)^h1 + C1^h1*E2*r2*(alpha2*s(i).xdata(j)^2)^h2 + E3*s(i).xdata(j)^h2*r1*(C1*alpha1*s(i).ydata(j))^h1 + E3*s(i).ydata(j)^h1*r2*(C2*alpha2*s(i).xdata(j))^h2 + C1^(2*h1)*C2^h2*E0*r1 + C1^h1*C2^(2*h2)*E0*r2 + C1^(2*h1)*E2*s(i).xdata(j)^h2*r1 + C2^(2*h2)*E1*s(i).ydata(j)^h1*r2 + E1*r2*(C2*s(i).xdata(j))^h2*(alpha1*s(i).ydata(j))^h1 + E2*r1*(C1*s(i).ydata(j))^h1*(alpha2*s(i).xdata(j))^h2 + C2^h2*E1*r1*(C1*s(i).ydata(j))^h1 + C1^h1*E2*r2*(C2*s(i).xdata(j))^h2)/(r1*(alpha2*s(i).xdata(j))^h2*(alpha1*s(i).ydata(j)^2)^h1 + r2*(alpha1*s(i).ydata(j))^h1*(alpha2*s(i).xdata(j)^2)^h2 + C2^h2*r1*(C1*alpha1*s(i).ydata(j))^h1 + C1^h1*r2*(C2*alpha2*s(i).xdata(j))^h2 + C2^h2*r1*(alpha1*s(i).ydata(j)^2)^h1 + C1^h1*r2*(alpha2*s(i).xdata(j)^2)^h2 + s(i).xdata(j)^h2*r1*(C1*alpha1*s(i).ydata(j))^h1 + s(i).ydata(j)^h1*r2*(C2*alpha2*s(i).xdata(j))^h2 + C1^(2*h1)*C2^h2*r1 + C1^h1*C2^(2*h2)*r2 + C1^(2*h1)*s(i).xdata(j)^h2*r1 + C2^(2*h2)*s(i).ydata(j)^h1*r2 + r1*(C1*s(i).ydata(j))^h1*(alpha2*s(i).xdata(j))^h2 + r2*(C2*s(i).xdata(j))^h2*(alpha1*s(i).ydata(j))^h1 + C2^h2*r1*(C1*s(i).ydata(j))^h1 + C1^h1*r2*(C2*s(i).xdata(j))^h2);
        end
        %plot3(s(i).xdata,s(i).ydata,cont,'color',map(i,:),'linewidth',1)
        plot3(s(i).xdata,s(i).ydata,cont,'k','linewidth',1)
    end
    set(gca,'ZGrid','on')
    scatter3(ones(1000,1)*min(d1),d2,zeros(1000,1),10,'k')
    scatter3(d1,max(d2)*ones(1000,1),zeros(1000,1),10,'k')
    scatter3(ones(1000,1)*max(d1),d2,zeros(1000,1),10,'k')
    scatter3(d1,min(d2)*ones(1000,1),zeros(1000,1),10,'k')
    
    xlabel('log(Drug2)')
    ylabel('log(Drug1)')
    zlabel('DIP rate')
    
    text(max(d1)-.3,max(d2)-.3,Ed(end,end)+.004,'E3')
    set(gca,'ztick',sort(unique([get(gca,'ztick'),0])))
    tix = get(gca,'ztick');
    tix_lab = cell(length(tix));
    tix_lab{tix==0}='0';
    set(gca,'zticklabels',tix_lab)
    rotate3d on
    close(popbox)



% --- Outputs from this function are returned to the command line.
function varargout = synGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider9_Callback(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider10_Callback(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider11_Callback(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.defaults
x = handles.defaults;
h1 = x(1);
h2 = x(2);
C1 = x(3);
C2 = x(4);
E0 = x(5);
E1 = x(6);
E2 = x(7);
alpha1 = x(8);
alpha2 = x(9);
beta = x(10);
trans_ratio = x(11);
scale_rate = x(12);
set(handles.slider6,'Value',h1);
set(handles.slider7,'Value',h2);
set(handles.slider4,'Value',log10(C1));
set(handles.slider5,'Value',log10(C2));
set(handles.slider1,'Value',E0);
set(handles.slider2,'Value',E1);
set(handles.slider3,'Value',E2);
set(handles.slider9,'Value',alpha1);
set(handles.slider10,'Value',alpha2);
set(handles.slider11,'Value',beta);
set(handles.slider12,'Value',trans_ratio);
set(handles.slider13,'Value',scale_rate);
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);


% --- Executes on slider movement.
function slider12_Callback(hObject, eventdata, handles)
% hObject    handle to slider12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider13_Callback(hObject, eventdata, handles)
% hObject    handle to slider13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = PlotFun(handles);
handles = modelPic(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
