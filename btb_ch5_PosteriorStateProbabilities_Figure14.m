%% Supplementary Material to "Beyond the Boundary - Chapter 5"
%
%%% Introduction
% This is a demonstration script accompanying the fifth chapter of my dissertation (Beyond the Boundary). Its purpose is to
% recreate Figure 5.1 and 5.2. This script should thus be used as an example to build your own scripts only. It is not a function or toolbox, not
% optimized for speed or functionality! Evaluate the code per Cell (using MATLAB's Cell Mode) and inspect the workspace to
% see what is going on.
%
% Posterior state probabilities - Figure 5.14

%% Data / Toolboxes / Scripts used
%
% * Fred's toolbox: https://github.com/FredHasselman/toolboxML
% * Data is available in the GithHub repository: <https://github.com/FredHasselman/BTB-Supplemental-Material-CHAPTER5>

%% Author / Version / License
%
% Repository: <https://github.com/FredHasselman/BTB-Supplemental-Material-CHAPTER5/ *BTB Chapter 5 on GitHub*>
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
load('PARAMSPACEsubj_2x21x21.mat')
load('paramspaceBDB_2x21x21.mat','BDB');

%% Get data

Nc = 10:1:30;
CS = 0:.05:1;
K = -1;
E= [-.5 0 .5];

s=tic;
for ss=1:length(LL)
 
 for c=1:6
  
  y=zeros(1,40);
  
   switch LL(ss).cond(c,1)
   case 1
    
    % Bak->DaK->BaK
    % Map 4th quadrant to resp=0 (BaK)
    % Map 2nd quadrant to resp=1 (DaK)
    y=2-LL(ss).resp(c,:);
    
   case 2
    
    % Dak->BaK->DaK
    % Map 4th quadrant to resp=1 (DaK)
    % Map 2nd quadrant to resp=0 (BaK)
    y=LL(ss).resp(c,:)+1;
    
   case 3
    
    % Random, sorted as Bak->DaK->BaK
    % Map 4th quadrant to resp=0 (BaK)
    % Map 2nd quadrant to resp=1 (DaK)
    y=2-LL(ss).resp(c,:);
    [O1,I1]=sort(LL(ss).stimnr(c,1:20),'ascend');
    [O2,I2]=sort(LL(ss).stimnr(c,21:40),'descend');
    y=y([I1 I2+20]);
    
  end    
  
     start=[0 1 0 0;0 0 0 1];
     nSTIMS=length(y);
     
     iE =find(LL(ss).llE(c)==E);
     iCS=find(LL(ss).llCS(c)==CS);
     iNc=find(LL(ss).llNc(c)==Nc);
     
     for n=2:nSTIMS
      
      EMISS=[0 0;BDB(iE,iCS,iNc).trial(n).EMISSmat+1e-5];
      TRANS=[0 start(y(n-1),:)+1e-5;zeros(4,1) BDB(iE,iCS,iNc).trial(n).TRANSmat+1e-5];
      [LL(ss).state(c).pp(:,n),~]= hmmdecode(y(n),TRANS,EMISS);
      
     end
     
     LL(ss).ppSUM(:,c)=sum(log(LL(ss).state(c).pp(2:5,:)+1e-05),2);
     LL(ss).ppMEAN(:,c)=mean(LL(ss).state(c).pp(2:5,:),2);
     LL(ss).ppMAX(:,c)=max(LL(ss).state(c).pp(2:5,:),[],2);
      
 end
 
end

en=toc(s);
disp(['total: ',num2str(en)]);

%% Get variables

Ac1_1=0;
Ac2_1=0;
Ac3_1=0;
Ac1_2=0;
Ac2_2=0;
Ac3_2=0;
Ac1_3=0;
Ac2_3=0;
Ac3_3=0;
AVE1bdb=1e-5;
AVE1dbd=1e-5;
AVE1rnd=1e-5;
AVE2bdb=1e-5;
AVE2dbd=1e-5;
AVE2rnd=1e-5;
AVE3bdb=1e-5;
AVE3dbd=1e-5;
AVE3rnd=1e-5;

Dc1_1=0;
Dc2_1=0;
Dc3_1=0;
Dc1_2=0;
Dc2_2=0;
Dc3_2=0;
Dc1_3=0;
Dc2_3=0;
Dc3_3=0;
DYS1bdb=1e-5;
DYS1dbd=1e-5;
DYS1rnd=1e-5;
DYS2bdb=1e-5;
DYS2dbd=1e-5;
DYS2rnd=1e-5;
DYS3bdb=1e-5;
DYS3dbd=1e-5;
DYS3rnd=1e-5;

A1=0;
A2=0;
A3=0;
D1=0;
D2=0;
D3=0;


for ss=1:length(LL)
 
 if LL(ss).age <= 125
  
  if LL(ss).dys==0
    A1=A1+1; 
    Acor1(A1,:) = [LL(ss).subjnr LL(ss).dys LL(ss).age LL(ss).klep LL(ss).dmt1 LL(ss).dmt2 LL(ss).dmt3 (LL(ss).ppSUM(1,:)+LL(ss).ppSUM(3,:)) LL(ss).llc ];

   for c=1:6
    
    switch LL(ss).llc(c)
     case 1
      AVE1bdb=AVE1bdb+log(LL(ss).state(c).pp(2:5,:)+1e-5);
      Ac1_1=Ac1_1+1;
     case 2
      tmp1=LL(ss).state(c).pp(3,:);
      tmp2=LL(ss).state(c).pp(5,:);
      tmp(1:4,:)=[LL(ss).state(c).pp(2,:);tmp2;LL(ss).state(c).pp(4,:);tmp1;];
      AVE1dbd=AVE1dbd+log(tmp+1e-5);
      Ac2_1=Ac2_1+1;
      clear tmp1 tmp2 tmp
     case 3
      AVE1rnd=AVE1rnd+log(LL(ss).state(c).pp(2:5,:)+1e-5);
      Ac3_1=Ac3_1+1;
    end
    
   end
   
  else
   
    D1=D1+1; 
    Dcor1(D1,:) = [LL(ss).subjnr LL(ss).dys LL(ss).age LL(ss).klep LL(ss).dmt1 LL(ss).dmt2 LL(ss).dmt3 (LL(ss).ppSUM(1,:)+LL(ss).ppSUM(3,:)) LL(ss).llc ];

   for c=1:6
    
    switch LL(ss).llc(c)
     case 1
      DYS1bdb=DYS1bdb+log(LL(ss).state(c).pp(2:5,:)+1e-5);
      Dc1_1=Dc1_1+1;
     case 2
      tmp1=LL(ss).state(c).pp(3,:);
      tmp2=LL(ss).state(c).pp(5,:);
      tmp(1:4,:)=[LL(ss).state(c).pp(2,:);tmp2;LL(ss).state(c).pp(4,:);tmp1;];
      DYS1dbd=DYS1dbd+log(tmp+1e-5)
      Dc2_1=Dc2_1+1;
      clear tmp1 tmp2 tmp
     case 3
      DYS1rnd=DYS1rnd+log(LL(ss).state(c).pp(2:5,:)+1e-5);
      Dc3_1=Dc3_1+1;
    end
    
   end
   
  end
 
 end
 
 
 if (LL(ss).age > 125) && (LL(ss).age <= 134)
  
  if LL(ss).dys==0
    A2=A2+1; 
    Acor2(A2,:) = [LL(ss).subjnr LL(ss).dys LL(ss).age LL(ss).klep LL(ss).dmt1 LL(ss).dmt2 LL(ss).dmt3 (LL(ss).ppSUM(1,:)+LL(ss).ppSUM(3,:)) LL(ss).llc];

   for c=1:6
    switch LL(ss).llc(c)
     case 1
     AVE2bdb=AVE2bdb+log(LL(ss).state(c).pp(2:5,:)+1e-5);
     Ac1_2=Ac1_2+1;
     case 2
      tmp1=LL(ss).state(c).pp(3,:);
      tmp2=LL(ss).state(c).pp(5,:);
      tmp(1:4,:)=[LL(ss).state(c).pp(2,:);tmp2;LL(ss).state(c).pp(4,:);tmp1;];
      AVE2dbd=AVE2dbd+log(tmp+1e-5)
      Ac2_2=Ac2_2+1;
      clear tmp1 tmp2 tmp
     case 3
      AVE2rnd=AVE2rnd+log(LL(ss).state(c).pp(2:5,:)+1e-5);
      Ac3_2=Ac3_2+1;
    end
    
   end
   
  else
    D2=D2+1; 
    Dcor2(D2,:) = [LL(ss).subjnr LL(ss).dys LL(ss).age LL(ss).klep LL(ss).dmt1 LL(ss).dmt2 LL(ss).dmt3 (LL(ss).ppSUM(1,:)+LL(ss).ppSUM(3,:)) LL(ss).llc];

   for c=1:6
    
    switch LL(ss).llc(c)
     case 1
      DYS2bdb=DYS2bdb+log(LL(ss).state(c).pp(2:5,:)+1e-5);
      Dc1_2=Dc1_2+1;
     case 2
      tmp1=LL(ss).state(c).pp(3,:);
      tmp2=LL(ss).state(c).pp(5,:);
      tmp(1:4,:)=[LL(ss).state(c).pp(2,:);tmp2;LL(ss).state(c).pp(4,:);tmp1;];
      DYS2dbd=DYS2dbd+log(tmp+1e-5);
      Dc2_2=Dc2_2+1;
      clear tmp1 tmp2 tmp
     case 3
      DYS2rnd=DYS2rnd+log(LL(ss).state(c).pp(2:5,:)+1e-5);
      Dc3_2=Dc3_2+1;
    end
    
   end
   
  end
  
 end
 
 if LL(ss).age > 134
  
  if LL(ss).dys==0
    A3=A3+1; 
    Acor3(A3,:) = [LL(ss).subjnr LL(ss).dys LL(ss).age LL(ss).klep LL(ss).dmt1 LL(ss).dmt2 LL(ss).dmt3 (LL(ss).ppSUM(1,:)+LL(ss).ppSUM(3,:)) LL(ss).llc];

   for c=1:6
 
    switch LL(ss).llc(c)
     case 1
      AVE3bdb=AVE3bdb+log(LL(ss).state(c).pp(2:5,:)+1e-5);
      Ac1_3=Ac1_3+1;
     case 2
      tmp1=LL(ss).state(c).pp(3,:);
      tmp2=LL(ss).state(c).pp(5,:);
      tmp(1:4,:)=[LL(ss).state(c).pp(2,:);tmp2;LL(ss).state(c).pp(4,:);tmp1;];
      AVE3dbd=AVE3dbd+log(tmp+1e-5);
      clear tmp1 tmp2 tmp
      Ac2_3=Ac2_3+1;
     case 3
      AVE3rnd=AVE3rnd+log(LL(ss).state(c).pp(2:5,:)+1e-5);
      Ac3_3=Ac3_3+1;
    end
    
   end
   
  else
    D3=D3+1; 
    Dcor3(D3,:) = [LL(ss).subjnr LL(ss).dys LL(ss).age LL(ss).klep LL(ss).dmt1 LL(ss).dmt2 LL(ss).dmt3 (LL(ss).ppSUM(1,:)+LL(ss).ppSUM(3,:)) LL(ss).llc];

   for c=1:6
    
    switch LL(ss).llc(c)
     case 1
      DYS3bdb=DYS3bdb+log(LL(ss).state(c).pp(2:5,:)+1e-5);
      Dc1_3=Dc1_3+1;
     case 2
      tmp1=LL(ss).state(c).pp(3,:);
      tmp2=LL(ss).state(c).pp(5,:);
      tmp(1:4,:)=[LL(ss).state(c).pp(2,:);tmp2;LL(ss).state(c).pp(4,:);tmp1;];
      DYS3dbd=DYS3dbd+log(tmp+1e-5)
      Dc2_3=Dc2_3+1;
      clear tmp1 tmp2 tmp
     case 3
      DYS3rnd=DYS3rnd+log(LL(ss).state(c).pp(2:5,:)+1e-5);
      Dc3_3=Dc3_3+1;
    end
    
   end
   
  end
  
 end
 
end


%% Figure 5.14
nr=8;
cmn=-1000;
cmx=0;

f0=figure;
maximize(f0)

C=gray(1200);
colormap(flipud(C(10:50:1150,:)))

stimnum={'','1','','','','','','','','','10','','','','','','','','','','20','','','','','','','','','','30','','','','','','','','','','40'};

subplot(nr,3,1)
imagesc(coll(AVE1bdb),[cmn cmx]) %+AVE2bdb+AVE3bdb))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('A1: Average readers - Sequential /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');



subplot(nr,3,2)
imagesc(coll(AVE1dbd),[cmn cmx]) %+(AVE2dbd)+(AVE3dbd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('A1: Average readers - Sequential /dAk/ >> /bAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

subplot(nr,3,3)
imagesc(coll(AVE1rnd),[cmn cmx]) %+(AVE2rnd)+(AVE3rnd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('A1: Average readers - Random sorted as /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

%---

subplot(nr,3,4)
imagesc(coll(DYS1bdb),[cmn cmx]) %+(DYS2bdb)+(DYS3bdb)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('D1: Dyslexic readers - Sequential /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

subplot(nr,3,5)
imagesc(coll(DYS1dbd),[cmn cmx]) %+(DYS2dbd)+(DYS3dbd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('D1: Dyslexic readers - Sequential /dAk/ >> /bAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

subplot(nr,3,6)
imagesc(coll(DYS1rnd),[cmn cmx]) %+(DYS2rnd)+(DYS3rnd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('D1: Dyslexic readers - Random sorted as /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

%-----

subplot(nr,3,10)
imagesc(coll(AVE2bdb),[cmn cmx]) %+AVE3bdb))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('A2: Average readers - Sequential /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');



subplot(nr,3,11)
imagesc(coll(AVE2dbd),[cmn cmx]) %+(AVE3dbd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('A2: Average readers - Sequential /dAk/ >> /bAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

subplot(nr,3,12)
imagesc(coll(AVE1rnd),[cmn cmx]) %+(AVE2rnd)+(AVE3rnd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('A2: Average readers - Random sorted as /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

%---

subplot(nr,3,13)
imagesc(coll(DYS2bdb),[cmn cmx]) %+(DYS3bdb)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('D2: Dyslexic readers - Sequential /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

subplot(nr,3,14)
imagesc(coll(DYS2dbd),[cmn cmx]) %+(DYS3dbd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('D2: Dyslexic readers - Sequential /dAk/ >> /bAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

subplot(nr,3,15)
imagesc(coll(DYS2rnd),[cmn cmx]) %+(DYS3rnd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('D2: Dyslexic readers - Random sorted as /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

%-----

subplot(nr,3,19)
imagesc(coll(AVE2bdb),[cmn cmx]) %+AVE3bdb))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('A3: Average readers - Sequential /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');



subplot(nr,3,20)
imagesc(coll(AVE2dbd),[cmn cmx]) %+(AVE3dbd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('A3: Average readers - Sequential /dAk/ >> /bAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

subplot(nr,3,21)
imagesc(coll(AVE1rnd),[cmn cmx]) %+(AVE2rnd)+(AVE3rnd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('A3: Average readers - Random sorted as /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

%---

subplot(nr,3,22)
imagesc(coll(DYS2bdb),[cmn cmx]) %+(DYS3bdb)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('D3: Dyslexic readers - Sequential /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

subplot(nr,3,23)
imagesc(coll(DYS2dbd),[cmn cmx]) %+(DYS3dbd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('D3: Dyslexic readers - Sequential /dAk/ >> /bAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');

subplot(nr,3,24)
imagesc(coll(DYS2rnd),[cmn cmx]) %+(DYS3rnd)))
shading faceted
axis equal
set(gca,'XTick',[0.5:1:40.5],'XTickLabel',stimnum,'YTick',[1:3],'YTickLabel',{'/dAk/ (Q2)','Spurious (Q1+Q3)','/bAk/ (Q4)'})
ylim([0.5 3.5])
xlim([1.5 40.5])
title('D3: Dyslexic readers - Random sorted as /bAk/ >> /dAk/')
xlabel('Stimulus number')
set(gca,'GridLineStyle','-','XGrid','on');



colorbar('Position',[0.05 .37 .015 .3]);
cbar_handle = findobj(f0,'tag','Colorbar');
set(cbar_handle,'YAxisLocation','left')
set(get(cbar_handle,'YLabel'),'String','Posterior -log Likelihood of States Q1-Q4')
set(gca,'GridLineStyle','-','XGrid','on');
% 
% 
% set(gcf, 'Color', 'w');
% export_fig 'Figure5_14' -eps -q101 -painters
