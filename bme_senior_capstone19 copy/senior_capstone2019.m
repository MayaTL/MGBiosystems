% File Name: senior_capstone2019.m
% Authors: Maya Laughton and Gonzalo Miranda
% Affiliation: Tufts University, Biomedical Engineering Dpt.
% Class: BME 7 & 8
% Date: Fall 2018, Spring 2019
% Project: Senior Capstone
% Study Title: Learning Mechanisms of the Pavlovian pathway in the
%              amygdala: a computational modeling and recording study
% Description: In Silico Fear Conditioning Model. Models the neural state 
%              of a virtual rodent undergoing a fear conditioning 
%              experimental procedure
% Associated files: 
%   network_math_v2.m - The calculations to determine neural state of brain
%                   under experimentation.
%   consolidate_v1.m - For memory consolidation. To test a post-conditioned
%                   rat. Increases LA-CEm connectivity. Set conductivity of
%                   cxt cells. 
%   subplot_code.m - Used for formating graphs in research report.

clear 

%%%%%%%%%%%%%%%%%%%%%%% Parameter Declarations %%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ExpData ActBL ActLA NLA NBL
global alphaLA betaLA etaLA zetaLA alphaBL betaBL etaBL zetaBL 
global a b c d
global psec ArMultLA ArMultBL E opThrsh pi maxsec
global Gi Gj Dp Dq HabitBL HabitLA  Habit
global Np Ncxi Ncxt Ncxp Nhc VpDist VqDist VpMax VqMax ViMax VjMax Gpq
global thrshRp thrshRi thrshRq thrshRj thrshXp endXp thrshXq endXq 
global thrshElig thrshXi gammaXi thrshXj gammaXj
global ilTau ilExp ilCmpndTau ilCmpndExp
global ilHxTau ilHxExp ilHxCmpndTau ilHxCmpndExp mx mxHx Gp_cemON
global Gpe Gie Gqe Gje Wp Ui Wq Uj Vq Vj Vp Bp Bq CaI Vi  
global is_run1 interrupt x Fsmrd Nsched nextstart Npoints interval Nrows
global recalcActs continueAt
global ilA ilB ilA1 ilB1 il
global ilHxA ilHxB ilHxA1 ilHxB1 ilHx
global conA conB conA1 conB1 con
global Fsmrd0 Aop AopPrv Gp_cem
global eligibleLA eligibleBL
global oldAq oldAu oldAx Ax oldRinput newRinput Au
global netErecipI
global  Hx Hs BLs naloxone Ishk FamAllCntxts
global is_consolidate
global CEMs LAs PFCs pcentHs

%Variables evaluated in network math
global inptCS1 inptA inptB  
global AcxpCS1 AcxpA AcxpB AcxpA1 AcxpB1
global rNcxpCS1 rNcxpA rNcxpB  rNcxpA1 rNcxpB1
global AcxiA AcxiB  
global rNcxiA rNcxiB  
global AcxtA AcxtB AcxtA1 AcxtB1
global rNcxtA rNcxtB rNcxtA1 rNcxtB1  
global AhcA AhcB AhcA1 AhcB1
global rNhcA rNhcB rNhcA1 rNcxt rNhcB1
global Ap Ai Ar Af Afhist
global rNcxp rNhc rNcxi 
global coefficientLA coefficientBL Act
global NT NP NH NI

%Variables used for graphs
global Wdata Eventdata Net_Actdata
global fig2on fig3on fig5on fig8on

%%%%%%%%%%%%%%%%%%%%%%% Parameter Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphaLA=.015; 
betaLA=.000003;
etaLA=.01;     
zetaLA=0;
alphaBL=.009;
betaBL=7e-7; 
etaBL=.01;  
zetaBL=.0000;

Gi=8;
Gj=8;
psec=.2;
ArMultLA=1;
ArMultBL=1;
Gpq=2;
Gp_cemON=1;
VpMax=80;  
VqMaxCa=80;
VqMax=66.67;
ViMax=30 ;
VjMax=30; 
opThrsh=0.2;
pi=1.25; 

NPcs=100;
NPcntxt=1; 
NPcntxt_cs=0;
NIcntxt=0;
NHcntxt=200;
NHcntxt_cs=250; 
NTcntxt=56;
NTcntxt_cs=250;
NPcntxt_0=0;     
NHcs=0;             
NIcs=0;
NTcntxt_0=0;
NIcntxt_cs=0;
NHcntxt_0=0;  
NIcntxt_0=0;  
Np=1;
                      
thrshRp=40;            
thrshRi=40;            
thrshRq=40;           
thrshRj=40;           
thrshXp=20; endXp=50; 
thrshXq=20; endXq=50;
thrshXi=0; gammaXi=13;
thrshXj=0; gammaXj=13;

a=.3;  
b=.66;
c=.9; %rise rate of freezing
d=.3; %fall rate of freezing (0 min; 1 max)
maxsec=0.4;
thrshElig=.05;
ilTau=.6; ilExp=2.5; ilCmpndTau=.1; ilCmpndExp=1;
ilHxTau=.2;
ilHxExp=2.5; ilHxCmpndTau=.1; ilHxCmpndExp=1;
mx=.9875; mxHx=.98;

E=100;
Ishk=100;

NP=[NPcntxt NPcntxt NPcs NPcntxt_cs NPcntxt_cs NPcntxt_0 NPcntxt_0 NPcntxt_0];
NH=[NHcntxt NHcntxt NHcs NHcntxt_cs NHcntxt_cs NHcntxt_0 NHcntxt_0 NHcntxt_0];
NI=[NIcntxt NIcntxt NIcs NIcntxt_cs NIcntxt_cs NIcntxt_0 NIcntxt_0 NIcntxt_0];
NT=[NTcntxt NTcntxt 0 NTcntxt_cs NTcntxt_cs NTcntxt_0 NTcntxt_0 NTcntxt_0];
NI=NI;
NH=NH;  
PFCs=0; LAs=0; BLs=0; CEMs=0;Hs=0;Hx=0;  
Dp=1; Dq=1;
Gp=0 ;Gq=0;
Kp=1/VpMax;
Kq1=1/VqMaxCa;
Kq2=1/VqMax;
Ki=1/ViMax;
Kj=1/VjMax;
ch=.2; %context representation habituation factor (<1)

%environmental contextual fear less effective not paired to CS
Habit=[ch ch 1 1 1 1 1 1];  
Habitcx=[Habit Habit Habit];                
HabitLA=Habitcx;          
HabitBL=[Habit 1];
Habit=[Habit Habit Habit Habit];  

Gp_cem=0;
naloxone=0;
pcentHs=100;
Nrows=10000;
Fsmrd0=zeros(Nrows,1);

FamAllCntxts=1;  
fig2on=0; fig3on=0; fig5on=0; fig8on=0;

%%%%%%%%%%%%%%%%%%%%%%% Read In Experiment Data %%%%%%%%%%%%%%%%%%%%%%%%%%%
ExpData = xlsread('test1.xlsx');

%%%%%%%%%%%%%%%%%%%%%%%% Initializations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_run1 = 1;
continueAt = 1;
interrupt=0;
Wdata=zeros(Nrows,6);
Eventdata=zeros(Nrows, 5);
Net_Actdata=zeros(Nrows,4);
Afhist=zeros(Nrows,1);
Fsmrd=zeros(Nrows, 1);
netErecipI = zeros(Nrows,4); 
Fsmrd0=zeros(Nrows,1);  
recalcActs=0;

if is_run1==1 %initialize all values
  Nsched = size(ExpData,1); 
  Npoints=0;
  is_consolidate=0;
  for j=1:Nsched   %Calculate Npoints for graphs.
    if ExpData(j,5)==1
      Npoints=Npoints+1 ; %Space in graph to place line indicating pause.
    else
      Npoints=Npoints+ExpData(j,6);  %Regular line.
    end  
  end
  oldAq=0; oldAu=0; oldAx=0; Ax=0; oldRinput=0; newRinput=0;Au=0; AopPrv=0;
  
  interval=0; nextstart=1;
  Act=zeros(1,8*4);  %Act is the full set of activations from which Vp and Vi are calculated.
  
  LAs=0;CEMs=0; BLs=0;Hx=0;Hs=0; PFCs=0;
  Gp_cem=0;
  ilA=0;ilB=0;ilA1=0;ilB1=0;
  ilHxA=0;ilHxB=0;ilHxA1=0;ilHxB1=0;
  conA=0;conB=0;conA1=0;conB1=0;
  con=[0 0 0 0 0 0 0 0];
  il= [0 0 0 0 0 0 0 0];  
  ilHx=il;  
  Aop=0; 
  Gpe=zeros(1,8*3); Gie=zeros(1,8*3);
  Gqe= zeros(1,8+1); Gje=zeros(1,8+1);
  inptCS1=0; inptA=0; inptB=0;  
  AcxpCS1=0; AcxpA=0; AcxpB=0; AcxpA1=0; AcxpB1=0; Scx=0; Lcx=1;
  rNcxpCS1=1; rNcxpA=1; rNcxpB=1; rNcxpA1=1; rNcxpB1=1; 
  AcxiA=0; AcxiB=0;  
  rNcxiA=1; rNcxiB=1;
  AcxtA=0; AcxtB=0; AcxttA1=0; AcxtB1=0; 
  rNcxtA=1; rNcxtB=1; rNcxtA1=1; rNcxtB1=1;
  AhcA=0; AhcB=0; AhcA1=0; AhcB1=0; Shc=0; Lhc=1;
  rNhcA=1; rNhcB=1; rNhcA1=1; rNhcB1=1;
  Ap=0; Ai=0;  Ar=0; rNr=1; rNus=1; Af=0; 
  eligibleLA=zeros(1,3*8);
  eligibleBL=zeros(1,8+1);
  Wp=zeros(1,3*8);
  Ui=zeros(1,3*8);
  Wq=zeros(1,8+1);
  Uj=zeros(1,8+1); 
  coefficientLA=zeros(1,3*8); coefficientBL=zeros(1,8+1);
  
  if FamAllCntxts==1
    if Hx==0 && Hs==0, ilA=1;il(1)=1;ilB=1;il(2)=1;il(3)=1;end
  end 
  
end

figure(7); openfigs(7)=1;
set(7,'Name','Progress')
set(7,'Units','normalized','Position',[.2 .2 .6 .1])
bars=[0 Npoints];
barh(bars)
if continueAt >1, close(4);close(2);close(3);close(8);close(5);close(9);end

%%%%%%%%%%%%%%%%%%%%%%%%%% Consolidation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%is_consolidate = 1;
if is_consolidate==1
    consolidate_v1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Network Calcuations %%%%%%%%%%%%%%%%%%%%%%%%%%%
network_math_v2

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if interval>=Npoints
  hh=findobj('Tag','run_continue');
  set(hh,'BackgroundColor','white');        
end    

BW=1; %~Black and White. Makes non-color graphs for publication.    
clear x
x=(1:interval);  
close(7)
MyColorOrder=...
  [1 0 0;... %red
  0 1 0;...  %green  
  0 0 1;...   %blue
  0 0 0;...   %black
  0.8  0  1;...   %purple
  0 1 1;...     %cyan
  0.8 0.8 0.8;];  %gray

figure(5) %Events
if BW==1,set(gcf,'DefaultAxesColorOrder',[0 0 0;.749 0 .749;1 0 0]),end
fig5on=1;
posit=get(5,'Position');
set(5,'Name','Events','Position',[posit(1) posit(2) posit(3) 75])
%To plot contexts as diff colors all on one line must displace, limit axes &
%to 0<1:
yA=Eventdata(1:interval,2)-1; 
yB=Eventdata(1:interval,2)-2;
y0=Eventdata(1:interval,2);
yCS1=Eventdata(1:interval,3)*.15+.2;
yUS=Eventdata(1:interval,4)*.15+.6;
plot(x,yA,'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5)
hold;
plot(x,yB,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);
plot(x,y0,'o',   'MarkerFaceColor','k','MarkerEdgeColor','k',   'MarkerSize',5);
plot(x,[yCS1 yUS]);
%Plot verticle lines at pauses (interrupts):
interruptsX=[];
interruptsY=[];
for p=1:interval
  if Eventdata(p,5)==1
    interruptsX=[interruptsX p p];
    interruptsY=[interruptsY 0 1];  
  end %if  
end %for  
plot(interruptsX,interruptsY)
ylim([0 .9]);
xlim([0 Npoints]);
savefig('Events.fig');

figure(2) %Vq and Vp
fig2on=1;
set(2,'Name','Vp (red) & Vq (blue)')
plot(x,Net_Actdata(1:interval,1));
plot(x,Net_Actdata(1:interval,1),'r');
hold
plot(x,Net_Actdata(1:interval,3),'b');
xlim([0 Npoints]);
savefig('Vp (red) & Vq (blue).fig');

figure(3)
fig3on=1;
set(3,'Name','Ap(red),Aq(blue),Af(green)')
plot(x,Net_Actdata(1:interval,2),'r');
hold
plot(x,Net_Actdata(1:interval,4),'b');
plot(x,Afhist(1:interval,1),'g');
xlim([0 Npoints]);
savefig('Ap(red),Aq(blue),Af(green).fig');

figure(8)
if BW==0
  USpoint='r';
  USpointEdge='r';
  Fcolor='b';
elseif BW==1 
  USpoint='w';
  USpointEdge='r';
  Fcolor='k';
end  
fig8on=1;
set(8,'Name','Graph of Fsmrd & Fcs')
clear nonusacts
nonusacts= Fsmrd(1:interval) - 2*(1-Eventdata(1:interval,4)); %don't know why orientation of Fsmrd is OK
hold;
plot(x,Fsmrd(1:interval),Fcolor) 
plot(x,nonusacts(1:interval),'o','MarkerFaceColor',USpoint,'MarkerEdgeColor',USpointEdge,...
  'MarkerSize',3)     
ylim([0 1]);
xlim([0 Npoints]);

savefig('Fsmrd.fig');

% Load saved figures
subplot1=hgload('Events.fig');
subplot2=hgload('Fsmrd.fig');
subplot3=hgload('Vp (red) & Vq (blue)');
subplot4=hgload('Ap(red),Aq(blue),Af(green)');

% Prepare subplots
figure

h(4)=subplot(4,1,1);
ylim([0 1]);
xlim([0 Npoints]);
ylabel('Firing Rate')

h(2)=subplot(4,1,2);
ylim([0 1]);
xlim([0 Npoints]);
ylabel('Freezing')

h(3)=subplot(4,1,3);
ylim([0 80]);
xlim([0 Npoints]);
ylabel('mV')

h(1)=subplot(4,1,4);
ylim([0 0.9]);
xlim([0 Npoints]);
ylabel('Events')
xlabel('Seconds')

% Paste figures on the subplots
copyobj(allchild(get(subplot1,'CurrentAxes')),h(1));

copyobj(allchild(get(subplot2,'CurrentAxes')),h(2));

copyobj(allchild(get(subplot3,'CurrentAxes')),h(3));

copyobj(allchild(get(subplot4,'CurrentAxes')),h(4));
