%% Supplementary Material to "Beyond the Boundary - Chapter 5"
%
%%% Introduction
% This is a demonstration script accompanying the fifth chapter of my dissertation (Beyond the Boundary). Its purpose is to
% recreate Figure 5.1 and 5.2. This script should thus be used as an example to build your own scripts only. It is not a function or toolbox, not
% optimized for speed or functionality! Evaluate the code per Cell (using MATLAB's Cell Mode) and inspect the workspace to
% see what is going on.
%
% Create Figures 5.1 - 5.2

%% Data / Toolboxes / Scripts used
%
% * Fred's toolbox: https://github.com/FredHasselman/toolboxML
% * Data is available in the GithHub repository: <https://github.com/FredHasselman/BTB-Supplemental-Material-CHAPTER5>

%% Author / Version / License
%
% Repository: <https://github.com/FredHasselman/BTB-Supplemental-Material-CHAPTER4/ *BTB Chapter 4 on GitHub*>
%
% Created by: <http://www.fredhasselman.com/ *Fred Hasselman*> / January 2011
% Affiliations: <http://www.ru.nl/bsi _Behavioural Science Institute_> - <http://www.ru.nl _Radboud University Nijmegen_>

%% PREP
%
% Uncomment next line to clear and close everything... detergent grade!
% omo

%Change ... to the path on your machine where you stored the files
% path=pwd;
% cd(path)
% load('btb_ch4_stimfeatures_ori.mat');



%% Figure 5.1
omo;


K=[0; 0.2; 0.4; 0.6; -.5];
X=[-2.5:.001:2.5];
V=zeros(length(K),length(X));
denom=40;

for k = 1:length(K);
 
 V(k,:) = (K(k).*X - (X.^2)/2 + (X.^4)/4)/denom;
 
 [well(k).m well(k).I]=pickpeaks(V(k,:),100,2,'throughs');
 
end;

stim = {'/bAk/','/dAk/'};
wght = {'normal','bold'};

% stim = {'Novel','Repeated'};
% wght = {'normal','bold'};
indX = [.6 -2]; indY = [-.06 -.06];
up = .007;


right=[3];
titl = {'n = 1','n = 2','n = 3','n = 4','n = 5'};
nSTIMS=5;
Nc = 3;
CS = 0;
K = -.2;
E= -1;
% 
% [nCS]=numel(CS);
% [nK] =numel(K);
% [nE] =numel(E);
% [~,nNc]=size(Nc);
% 
% X = [-2:.001:2];
% 
% e=1;
% 
% s=tic;
% for n = 1:nSTIMS
%  
%  if (n-Nc)<0
%   Phi=0;
%  else
%   Phi=1;
%  end
%  
%  k(n) = K + n/Nc + (E(e)*Phi);
%  V(n).pot = k(n).*X - (X.^2)./2 + (X.^4)./4;
%  
%  [well(n).m well(n).I]=pickpeaks(V(n).pot,100,2,'throughs');
%  
% end
% 
% stim = {'Novel','Repeated'};
% wght = {'bold','normal'};
% indX = [.6 -2]; indY = [-.1 -.1];
% up = .08;
% 
% right=[5];
% titl = {'n = 1','n = 2','n_c = 3','n = 4','n = 5'};

f0=figure;
maximize(f0);

for n=1:Nc
 subplot(1,5,n)
 plot(X,V(n).pot,'k','LineWidth',2); hold on;
  ylim([-1 1]); xlim([-2.5 2.5]);
 axis square
 
 if n==Nc
  lr=1;
 lbl=1;
 else
  lr=2;
  
 lbl=2;
 end
 
 
%  if (X(well(n).I(lr)) > 0)
%   lbl=1;
%  else
%   lbl=2;
%  end
%  

% text(x(well(n).I(lr))-.6,well(n).m(lr)-.01,stim,'FontName','courier','FontWeight',wght);
 
text(indX(lbl),indY(lbl),stim{lbl},'FontName','courier',...
 'FontWeight',wght{lbl});

 oa=oaxes([0,0,0],'LineWidth',1,'XTickLabel',[],'YTickLabel',[],'Title',titl{n},...
  'TickLength',[2 2],'ArrowAspectRatio',1,'YLabel',{'V(x)',''},'XLabel',{'','x'});
 
 %oa=oaxes([0,0,0],'Title',titl{n},'YLabel',{'V(x)',''},'XLabel',{'','x'});
 
 
 plot(X(well(n).I(1)),well(n).m(1)+up,'ok','MarkerFaceColor',[.5 .5 .5],...
  'MarkerSize',12);

 
 hold on;
 
end
 plot(X(well(n).I(2)),well(n).m(2)+up,'ok','MarkerFaceColor',[.8 .8 .8],...
  'MarkerSize',12);

% Uncomment if you want to save a figure
% grab('btb_ch5_Figure1',0);



%% Figure 5.2 model with varying k
omo;

step=7;
denom=40;
lu=linspace(0,1,step);
L=[lu lu];
Le=lu(end);
E=1.2;
k(1)=-1;
X=[-2.5:.001:2.5];
nc = round((length(L)/2));
V=zeros(length(L),length(X));


V(1,:) = (k(1).*X - (X.^2)/2 + (X.^4)/4)/denom;
[well(1).m well(1).I]=pickpeaks(V(1,:),100,2,'throughs');

for n = 2:length(L);
 
 if (n-nc) < 0
  Phi=0;
 else
  Phi=1;
 end
 
 k(n) = k(1) + L(n) + E/2 + E*Phi*(L(n)-Le);
 
 V(n,:) = (k(n).*X - (X.^2)/2 + (X.^4)/4)/denom;
 
 [well(n).m well(n).I]=pickpeaks(V(n,:),100,2,'throughs');
 
 
end;

right=[4;5;12];

titl = {'$$1\Rightarrow$$','$$2\Rightarrow$$','$$3\Rightarrow$$','$$4\Rightarrow$$',...
 '$$5\Rightarrow$$','$$6\Rightarrow$$','$$7\Rightarrow$$','$$1\Leftarrow$$',...
 '$$2\Leftarrow$$','$$3\Leftarrow$$','$$4\Leftarrow$$','$$5\Leftarrow$$',...
 '$$6\Leftarrow$$','$$7\Leftarrow$$'};

f0=figure;
maximize(f0);

stim = {'/bAk/','/dAk/'};
wght = {'normal','bold'};
indX = [.6 -2]; indY = [-.07 -.07];
up = .011;

for n=1:length(L)
 subplot(3,nc,n)
 plot(X,V(n,:),'k','LineWidth',2); hold on;
 ylim([-.08 .08]); xlim([-2.5 2.5]);
 axis square;
 
 if ismember(n,right)
  lr=2;
 else
  lr=1;
 end
 
 if (X(well(n).I(lr)) > 0)
  lbl=1;
 else
  lbl=2;
 end
 
text(indX(lbl),indY(lbl),stim{lbl},'FontName','courier',...
 'FontWeight',wght{lbl});
 

 oa=oaxes([0,0,0],'LineWidth',1,'XTickLabel',[],'YTickLabel',[],...
  'TickLength',[2 2],'ArrowAspectRatio',1,'YLabel',{'V(x)',''},'XLabel',{'','x'});
 
 text('Interpreter','latex','String',titl{n},'Position',[-.3 .09]);
 
 plot(X(well(n).I(lr)),well(n).m(lr)+up,'ok','MarkerFaceColor',[.5 .5 .5],...
  'MarkerSize',12);
 
 hold on;
 
end

jump = [6;11]; from = jump + [-1;1]; lr = [2; 2];
x1 = [.76 .73; .505 .535];
y1 = [.855 .825; .555 .525];

for i = 1:length(jump)
 
 subplot(3,nc,jump(i))
 
 plot(X(well(from(i)).I(lr(i)))+.05,V(jump(i),well(from(i)).I(lr(i)))+up,'o',...
  'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.7 .7 .7],'MarkerSize',12);
 hold on;
 
 text(indX(i),indY(i),stim{i},'FontName','courier',...
 'FontWeight',wght{i});

 annotation('arrow',x1(i,:),y1(i,:),'HeadStyle','plain','HeadWidth',6,'HeadLength',6,...
  'LineWidth',1.1);
 
end

% Uncomment if you want to save a figure
% grab('btb_ch5_Figure2',0);

%% Figure 5.3
omo;

coupD = 0; coupP = .3; coupN = 0.7;
k = 0;

[x,y] = meshgrid([-4:1:4]);
r2 = x.^2 + y.^2; r4 = x.^4 + y.^4;
z1N = k*(x - y) - 1.3*r2 + .13*r4 + coupN*(x.*y);
z1D = k*(x - y) - 1.3*r2 + .1*r4 + coupD*(x.*y);
z1P = k*(x - y) - 1.3*r2 + .1*r4 + coupP*(x.*y);

[a,b] = meshgrid([-4:.2:4]);
p2 = a.^2 + b.^2; p4 = a.^4 + b.^4;
z2N = k*(a - b) - 1.3*p2 + .13*p4 + coupN*(a.*b);
z2D = k*(a - b) - 1.3*p2 + .1*p4 + coupD*(a.*b);
z2P = k*(a - b) - 1.3*p2 + .1*p4 + coupP*(a.*b);

subplot(1,3,1)
[dXD,dYD]=gradient(-z1D); hold on;
quiver(x,y,dXD,dYD,'Color',[0.2 0.2 0.2]); hold on;
contour(a,b,z2D); hold on;
surface(a,b,z2D+40);
view(-45,26); grid off; axis off;
axis square;
title('PRELINGUITSIC CHILD (No coupling)');

subplot(1,3,2)
[dXP,dYP]=gradient(-z1P); hold on;
quiver(x,y,dXP,dYP,'Color',[0.2 0.2 0.2]); hold on;
contour(a,b,z2P); hold on;
surface(a,b,z2P+44);
view(-45,26); grid off; axis off;
axis square;
title('DYSLEXIC READER (Weak coupling)');

subplot(1,3,3)
[dXN,dYN]=gradient(-z1N); hold on;
quiver(x,y,dXN,dYN,'Color',[0.2 0.2 0.2]); hold on;
contour(a,b,z2N); hold on;
surface(a,b,z2N+60);
view(-45,26); grid off; axis off;
axis square;
title('AVERAGE READER (Strong coupling)');

colormap(gray)


% Uncomment if you want to save a figure
% grab('btb_ch5_Figure3',0);


%% Figure 5.4
omo;

f0=figure;
maximize(f0);

CS  = 0;
k   = [-1 -.5 0 .5 1];
ttl = {'k = -1','k < -k_{c}','k = 0','k > k_{c}','k = 1'};

cm = [linspace(.2,1,64)' linspace(.2,1,64)' linspace(.2,1,64)'];
colormap(cm);

for j=1:length(k)
 
 [x,y] = meshgrid([-2.4:.2:2.4]);
 z = k(j)*(x-y) - .5*(x.^2 + y.^2) + .1*(x.^4 + y.^4) + CS*(x.*y);
 
 subplot(3,length(k),j)
 
 surface(x,y,z+40,'EdgeColor','none');
 view(-45,26); grid off; axis tight; axis off; axis square;
 
 title(ttl{j});
 
 
 subplot(3,length(k),j+length(k))
 
 [x,y] = meshgrid([-2.4:.1:2.4]);
 z = k(j)*(x-y) - .5*(x.^2 + y.^2) + .1*(x.^4 + y.^4) + CS*(x.*y);
 cors = x-y;
 line(x-y,z,'Color',[0.7 0.7 0.7]);
 ylim([-5 5]);
 hold on
 axis on
 
 cnt = 0;
 
 for i=1:length(x)
  z45(i) = z((length(x)+1)-i,i);
  xy(i) = cors((length(x)+1)-i,i);
 end
 
 lz = z45;
 up=[.5 .33 0.36 .38 .5];
 lr=[1 1 2 1 1];
 
 stim = {'/bAk/','/dAk/','?'};
 wght = {'normal','bold','bold'};
 
 
 plot(xy,lz,'-k','Linewidth',3);
 grid off; axis off;
 axis square;
 
 [well(j).m well(j).I]=pickpeaks(lz,10,2,'throughs');
 
 if j==4
  plot(xy(well(j).I(2))+.3,well(j).m(2)+up(j),'ok','MarkerEdgeColor',[.5 .5 .5],...
   'MarkerFaceColor',[.7 .7 .7],'MarkerSize',10); hold on;
  
  plot(0,-1,'ok','MarkerFaceColor',[.5 .5 .5],'MarkerSize',10); hold on;
  
  text(2.2,-5.5,stim{1},'FontName','courier','FontWeight',wght{1});hold on;
  
  annotation('arrow',[.71 .69],[.51 .48],'HeadStyle','plain','HeadWidth',...
   6,'HeadLength',6,'LineWidth',1.1);hold on;
  
 else
  plot(xy(well(j).I(lr(j))),well(j).m(lr(j))+up(j),'ok','MarkerFaceColor',[.5 .5 .5],...
   'MarkerSize',10); hold on;
 end
 
 lbl  = [1 1 1 3 2];
 indX = [2.8 2.6 2.4 -0.1 -4.8];
 indY = [-5.5 -5.5 -5.5 -5.5 -5.5];
 
 text(indX(j),indY(j),stim{lbl(j)},'FontName','courier','FontWeight',wght{lbl(j)});
 
 subplot(3,length(k),j+length(k)*2)
 
 lms=2.0;
 
 if j==round(length(k)/2)
  stp=.8;
  arw=.8;
  lin=15;
 else
  stp=.5;
  arw=.9;
  lin=25;
 end
 
 [xc,yc] = meshgrid([-lms:.01:lms]);
 zc = k(j)*(xc-yc) - .5*(xc.^2 + yc.^2) + .1*(xc.^4 + yc.^4) + CS*(xc.*yc);
 contour(xc,yc,zc,lin); hold on;
 
 indX = [ 1.85  1.8  1.6 1.2 -1.85];
 indY = [-1.85 -1.8 -1.6 1.8  1.85];
 
 if j==4
  plot(1,-1,'ok','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[.5 .5 .5],...
   'MarkerSize',10); hold on;
 end
 
 plot(indX(j),indY(j),'ok','MarkerFaceColor',[.5 .5 .5],'MarkerSize',10); hold on;
 
 [xq,yq] = meshgrid([-lms:stp:lms]);
 zq = k(j)*(xq-yq) - .5*(xq.^2 + yq.^2) + .1*(xq.^4 + yq.^4) + CS*(xq.*yq);
 [dX,dY]=gradient(-zq); hold on;
 quiver(xq,yq,dX,dY,arw,'Color','k','LineWidth',1); hold on;
 
 set(gca,'XTick',[],'YTick',[]);
 box on
 axis tight
 axis square
 
end


% Uncomment if you want to save a figure
% grab('btb_ch5_Figure4',0);

%% Figure 5.5 is a .tex file

%% Figure 5.6a and 5.6b
load('btb_ch5_ROCdata.mat');

fig1 = figure(1);
maximize(fig1);
a(1) = subplot(2,4,[1 2]);
a(2) = subplot(2,4,3);
a(3) = subplot(2,4,4);
a(4) = subplot(2,4,[5 6]);
a(5) = subplot(2,4,7);
a(6) = subplot(2,4,8);

fig2 = figure(2);
maximize(fig2);
a(7) = subplot(2,4,[1 2]);
a(8) = subplot(2,4,3);
a(9) = subplot(2,4,4);
a(10) = subplot(2,4,[5 6]);
a(11) = subplot(2,4,7);
a(12) = subplot(2,4,8);

set([fig1, fig2], 'NextPlot','add');
set(a, 'NextPlot','add');

axPRED=[1 4 7 10];
axROC =[2 5 8 11];
axPR  =[3 6 9 12];

str={'ok','vk','^k','dk'};
lns ={'-','-.'};
mclr=[0 0 0;.5 .5 .5];

ttlPRED={'Predicted probability - None','Predicted probability - Slowed Down',...
  'Predicted probability - Amplified','Predicted probability - Both'};
ttlROC={'ROC - None','ROC - Slowed Down','ROC - Amplified','ROC - Both'};
ttlPR={'PR - None','PR - Slowed Down','PR - Amplified','PR - Both'};

d=[-.1 .1];

for man=1:4
  
  axes(a(axPRED(man)))
  
  % Logit predictions
  for grp=1:2
    h(grp).h=ploterr([1:7]+d(grp),ROCdata(grp,man).fixedPRED(:,1),[],...
      {ROCdata(grp,man).fixedPRED(:,2) ROCdata(grp,man).fixedPRED(:,3)},...
      [lns{grp},str{man}],'hhy',0.5);hold on;
     set(h(grp).h(1,1),'MarkerFaceColor',mclr(grp,:),'MarkerSize',12) 
    
  end

  legend([h(1).h(1,1) h(2).h(1,1)],'Average reader','Dyslexic reader','Location','NorthEast');
   
  set(gca,'Xlim',[0.5 7.5],'Ylim',[0 1]);
  xlabel('Stimulus Pair'); ylabel('Predicted Probability for Perceiving DIFFERENT');
  title(ttlPRED{man});
  line(xlim,[0.5 0.5],'Color',[0.7 0.7 0.7],'LineStyle','--'); hold on;
  
  clear h
  
  %ROC
  axes(a(axROC(man)))
  
  plot([-.05 1.05],[-.05 1.05],'--','Color',[.7 .7 .7]);hold on;
  plot([-.05 .5],[1.05 .5],'-','Color',[.7 .7 .7]); hold on;
  
  for grp=1:2
    
    hr=ploterrb(ROCdata(grp,man).X(:,1),ROCdata(grp,man).Y(:,1),...
      {ROCdata(grp,man).X(:,2) ROCdata(grp,man).X(:,3)},...
      {ROCdata(grp,man).Y(:,2) ROCdata(grp,man).Y(:,3)},str{man}); hold on;
    set(hr,'MarkerSize',14,'MarkerFaceColor',mclr(grp,:),'Color','k');

    plot(ROCdata(grp,man).Xa(:,1),ROCdata(grp,man).Ya(:,1),...
      lns{grp},'Color','k'); hold on;
    
    
    Thresholds=ROCdata(grp,man).obsCIrnd(1:7,1);
     for sp = 1:numel(Thresholds)
      ind=find(Thresholds(sp)==ROCdata(grp,man).T,1,'first');
      text(ROCdata(grp,man).X(ind,1)-0.01,ROCdata(grp,man).Y(ind,1),...
      num2str(sp),'Color','w','FontSize',12);
    end
 
  end
  
  axis square;
  set(gca,'Xlim',[-0.05 1.05],'Ylim',[-0.05 1.05],'YTick', 0:.1:1,...
    'XTick', 0:.1:1,'Box','on');
  xlabel('False Positive Rate (1-Specificity)');ylabel('True Positive Rate (Sensitivity)');
  title(ttlROC{man});
  legend('Random','Unbiased','Average Reader','Dyslexic Reader','Location','SouthEast');
  
  
  %PR
  axes(a(axPR(man)))
  
  plot([-0.05 1.05],fliplr([-0.05 1.05]),'--','Color',[.7 .7 .7]);hold on;
  plot([.5 1.05],[.5 1.05],'-','Color',[.7 .7 .7]);
  
  for grp=1:2
    
    hpr=ploterrb(ROCdata(grp,man).Xpr(:,1),ROCdata(grp,man).Ypr(:,1),...
      {ROCdata(grp,man).Xpr(:,2) ROCdata(grp,man).Xpr(:,3)},...
      {ROCdata(grp,man).Ypr(:,2) ROCdata(grp,man).Ypr(:,3)},str{man});
    set(hpr,'MarkerSize',14,'MarkerFaceColor',mclr(grp,:),'Color','k');
    
    plot(ROCdata(grp,man).Xapr(:,1),ROCdata(grp,man).Yapr(:,1),...
      lns{grp},'Color','k'); hold on;
    
     Thresholds=ROCdata(grp,man).obsCIrnd(1:7,1);
   
    for sp = 1:numel(Thresholds)
      ind=find(Thresholds(sp)==ROCdata(grp,man).Tpr,1,'first');
      text(ROCdata(grp,man).Xpr(ind,1)-0.01,ROCdata(grp,man).Ypr(ind,1),...
      num2str(sp),'Color','w','FontSize',12);
    end
    
  end
  
  axis square
  set(gca,'Xlim',[-0.05 1.05],'Ylim',[-0.05 1.05],'YTick', 0:.1:1,...
    'XTick', 0:.1:1,'Box','on');
  xlabel('True Positive Rate (Recall)');ylabel('Positive Predictive Value (Precision)');
  title(ttlPR{man});
  legend('Random','Unbiased','Average Reader','Dyslexic Reader','Location','SouthWest');
  
  
end

% figure(fig1)
% grab('btb_ch5_Figure6a',0)
% 
% figure(fig2)
% grab('btb_ch5_Figure6b',0)

%% Table 5.2

%Tables of AUC
for man=1:4 
  for grp=1:2
    tabROC{man,grp} = [num2str(ROCdata(grp,man).AUC(1),'%1.2f'),...
      ' (',num2str(ROCdata(grp,man).AUC(2),'%1.2f'),',',num2str(ROCdata(grp,man).AUC(3),'%1.2f'),')'];
    tabPR{man,grp} = [num2str(ROCdata(grp,man).AUCpr(1),'%1.2f'),...
      ' (',num2str(ROCdata(grp,man).AUCpr(2),'%1.2f'),',',num2str(ROCdata(grp,man).AUCpr(3),'%1.2f'),')'];
  end
end

tabROC
tabPR

%% Figure 5.15
omo;

nSTIMS=5;
Nc = 5;
CS = 0;
K = -.2;
E= -1;

[nCS]=numel(CS);
[nK] =numel(K);
[nE] =numel(E);
[~,nNc]=size(Nc);

X = [-2:.001:2];

e=1;

s=tic;
for n = 1:nSTIMS
 
 if (n-Nc)<0
  Phi=0;
 else
  Phi=1;
 end
 
 k(n) = K + n/Nc + (E(e)*Phi);
 V(n).pot = k(n).*X - (X.^2)./2 + (X.^4)./4;
 
 [well(n).m well(n).I]=pickpeaks(V(n).pot,100,2,'throughs');
 
end

stim = {'Novel','Repeated'};
wght = {'bold','normal'};
indX = [.6 -2]; indY = [-.1 -.1];
up = .08;

right=[5];
titl = {'n = 1','n = 2','n = 3','n = 4','n_c = 5'};

f0=figure;
maximize(f0);

for n=1:Nc
 subplot(1,5,n)
 plot(X,V(n).pot,'k','LineWidth',2); hold on;
  ylim([-1 1]); xlim([-2.5 2.5]);
 axis square
 
 if n==Nc
  lr=1;
 lbl=1;
 else
  lr=2;
  
 lbl=2;
 end
 
text(indX(lbl),indY(lbl),stim{lbl},'FontName','courier',...
 'FontWeight',wght{lbl});

 oa=oaxes([0,0,0],'LineWidth',1,'XTickLabel',[],'YTickLabel',[],'Title',titl{n},...
  'TickLength',[2 2],'ArrowAspectRatio',1,'YLabel',{'V(x)',''},'XLabel',{'','x'});
  
 plot(X(well(n).I(1)),well(n).m(1)+up,'ok','MarkerFaceColor',[.5 .5 .5],...
  'MarkerSize',12);

 hold on;
 
end
 plot(X(well(n).I(2)),well(n).m(2)+up,'ok','MarkerFaceColor',[.8 .8 .8],...
  'MarkerSize',12);

% Uncomment if you want to save a figure
% grab('btb_ch5_Figure15',0);

