function consolidate_v1
% File Name: consolidate_v1
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
%   senior_capstone2019.m - Initializes parameter values. Calls
%                         network_math_v1, consolidate.m, subplot_code.m          
%   network_math_v2.m - The calculations to determine neural state of brain
%                       under experimentation.
%   subplot_code.m - Used for formating graphs in research report.

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

if PFCs==1, return,end
 
NLA=[NP NI NT];
NBL=[NH 1];

ilHx=il;
 
%Calculate Gpe's for cxt inputs of cntxts and conjunctions
phi=E/VpMax;
for j=1:8  %LA-BL not plastic
  if j==3, continue,end
  if il(j)>0 && BLs==0 && LAs==0
    nBL=NBL(j)*HabitBL(j)*il(j)^2;
    Aq=(E/VqMax)*nBL*Gqe(j)/(1+nBL*Gqe(j));
    nLA=NLA(8*2+j)*HabitLA(8*2+j)*ilHx(j)^2;
    Gpe(8*2+j)=(1/nLA)*Aq/(phi*(1+Gpq)-Aq*(1+Gpq*phi));
  end
end
 
%Calcuclate Gie's for cxt inputs cntxts and conjunctions
gam=E*Gi/ViMax;
for j=1:8  
  if j==3, continue,end
  if il(j)>0 && BLs==0 && LAs==0
    nBL=NBL(j)*HabitBL(j)*il(j)^2;
    nLA=NLA(8*2+j)*HabitLA(8*2+j)*ilHx(j)^2;
    Aj=(E/VjMax)*nBL*Gje(j)/(1+nBL*Gje(j)); 
    Aq=1/(1+Gj*Aj);
    Gie(8*2+j)=(1/nLA)*(1-Aq)*(1+Gpq)/(Aq*gam-(1-Aq)*(1+Gpq));
    Gje(j)=0;
  end
end
 
for j=1:8  
  if j==3, continue,end
  if il(j)>0 
    con(j)=1;
    il(j)=0;
  end  
end
end
