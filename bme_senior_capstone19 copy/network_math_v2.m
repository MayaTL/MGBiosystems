function network_math_v2
% File Name: network_math_v2
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
%   sneior_capstone2019.m - Initializes parameter values. Calls
%                         network_math_v1, consolidate.m, subplot_code.m          
%   consolidate_v1.m - For memory consolidation. To test a post-conditioned
%                   rat. Increases LA-CEm connectivity. Set conductivity of
%                   cxt cells. 
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

%Interprets lines of Experimental Schedule Array (ExpData)
for j=continueAt :Nsched 
  cntxt=ExpData(j,2); 
  cs1=ExpData(j,3); 
  us=ExpData(j,4); 
  interrupt=ExpData(j,5);
  dur=ExpData(j,6);
  if interrupt==1
    continueAt=j+1;
    interval=interval+1; 
    Eventdata(interval,1:5)=[interval cntxt cs1 us interrupt];
    Wdata(interval,1:5)=Wdata(interval-1,1:5); 
    break    %BREAK 
  end
  %Calculate coefficient. 
  %Transform ExpData info to inpt form for convenience in following calculations:
  if cs1==1, inptCS1=1; else inptCS1=0;end
  if cntxt==1, inptA=1; else inptA=0;end
  if cntxt==2, inptB=1; else inptB=0;end
  
  %Obtain consolidation status (0 or 1)and extent of implicit learning mediated by hippocampus, 
  %as well as implicit learning mediated by representation cortex (when hippocampus 
  %has been ablated)
  conA=con(1);conB=con(2);conA1=con(6);conB1=con(8);
  ilA=il(1);ilB=il(2);ilA1=il(6);ilB1=il(8);
  ilHxA=ilHx(1);ilHxB=ilHx(2);ilHxA1=ilHx(6);ilHxB1=ilHx(8);

  %The Calculation Cycle (W_INDEPENDENT EVALUATIONS - Steps 1 to 4)
        
  %Step 1: Calculate firing rate and ratio of neurons firing (to the
  %total number in the cortex) for principal cortical neurons
  if inptCS1==1, AcxpCS1=1; rNcxpCS1=1; else AcxpCS1=0;rNcxpCS1=0;end
  if inptA==1, AcxpA=1; rNcxpA=1; else AcxpA=0;rNcxpA=0;end
  if inptB==1, AcxpB=1; rNcxpB=1; else AcxpB=0;rNcxpB=0;end
  AcxpA1=AcxpA*AcxpCS1;rNcxpA1=AcxpA1; 
  AcxpB1=AcxpB*AcxpCS1;rNcxpB1=AcxpB1; 

  %Step 2: Calculate firing rate and ratio of neurons firing (to the
  %total number in the hippocampus) for hippocampal neurons
  AhcA=ilA*(1-conA)*(1-Hx)*(1-Hs*(pcentHs/100))*AcxpA*rNcxpA; 
  rNhcA=AhcA;
  AhcB=ilB*(1-conB)*(1-Hx)*(1-Hs*(pcentHs/100))*AcxpB*rNcxpB; 
  rNhcB=AhcB;
  AhcA1=ilA1*(1-conA1)*(1-Hx)*(1-Hs*(pcentHs/100))*sqrt(AcxpA*AcxpCS1)*sqrt(rNcxpA*rNcxpCS1); 
  rNhcA1=AhcA1;
  AhcB1=ilB1*(1-conB1)*(1-Hx)*(1-Hs*(pcentHs/100))*sqrt(AcxpB*AcxpCS1)*sqrt(rNcxpB*rNcxpCS1); 
  rNhcB1=AhcB1;

  %Step 3: Calculate firing rate and ratio of cortical inhibitory neurons (to
  %total number of cortical neurons)
  AcxiA=ilA*((1-conA)*AhcA + linsig(conA+ilHxA  ,0,1)*AcxpA)*((1-conA)*rNhcA + conA*rNcxpA); 
  rNcxiA=AcxiA;
  AcxiB=ilB*((1-conB)*AhcB + linsig(conB+ilHxB  ,0,1)*AcxpB)*((1-conB)*rNhcB + conB*rNcxpB); 
  rNcxiB=AcxiB;
  
  %Step 4: Caculate firing rate and ratio of neurons firing for transformed
  %cortical neurons
  AcxtA=linsig(conA+ilHxA,0,1)*AcxpA*rNcxpA; 
  rNcxtA=AcxtA;
  AcxtB=linsig(conB+ilHxB,0,1)*AcxpB*rNcxpB; 
  rNcxtB=AcxtB;
  AcxtA1=linsig(conA1+ilHxA1,0,1)*sqrt(AcxpA*AcxpCS1)*sqrt(rNcxpA*rNcxpCS1);
  rNcxtA1=AcxtA1;
  AcxtB1=linsig(conB1+ilHxB1,0,1)*sqrt(AcxpB*AcxpCS1)*sqrt(rNcxpB*rNcxpCS1);
  rNcxtB1=AcxtB1;

  if PFCs==1
    AcxtA=0;rNcxtA=0;
    AcxtB=0;rNcxtB=0;
    AcxtA1=0;rNcxtA1=0;
    AcxtB1=0;rNcxtB1=0;
  end

  rNcxp=[rNcxpA rNcxpB rNcxpCS1 rNcxpA1 rNcxpB1 0 0 0];
  Ncxp=rNcxp.*NP; 
  rNhc=[rNhcA rNhcB 0 rNhcA1 rNhcB1 0 0 0];
  Nhc=rNhc.*NH;
  rNcxi=[rNcxiA rNcxiB 0 0 0 0 0 0];
  Ncxi=rNcxi.*NI;
  rNcxt=[rNcxtA rNcxtB 0 rNcxtA1 rNcxtB1 0 0 0];
  Ncxt=rNcxt.*NT;

  NLAactive=[Ncxp Ncxi Ncxt]; 
  NBLactive=[Nhc Np];
   

  ActLA=(1-LAs)*[AcxpA AcxpB AcxpCS1 AcxpA1 AcxpB1 0 0 0 ...   
                 AcxiA AcxiB   0       0      0    0 0 0 ...
                 AcxtA AcxtB   0     AcxtA1 AcxtB1 0 0 0];   

  for k=1:3*8
    if NLAactive(k)==0,ActLA(k)=0;end
  end

  coefficientLA=ActLA.*NLAactive.*HabitLA;
  
  for i= 1 : dur  %This repeats the current line of ExpData dur times.
    interval=interval+1; %update calculation interval ("one-sec" time-step)
    %Write event data to data array:
    Eventdata(interval,1:5)=[interval cntxt cs1 us interrupt];

    %If implicit learning has occurred during the previous interval,we need
    %to do W_INDEPENDENT EVALUATIONS again before evaluating amygdala.
    if recalcActs==1 
        %Repeated - Obtain consolidation status (0 or 1)and extent of implicit learning mediated by hippocampus, 
        %as well as implicit learning mediated by representation cortex (when hippocampus 
        %has been ablated)
        conA=con(1);conB=con(2);conA1=con(4);conB1=con(5);
        ilA=il(1);ilB=il(2);ilA1=il(4);ilB1=il(5);
        ilHxA=ilHx(1);ilHxB=ilHx(2);ilHxA1=ilHx(4);ilHxB1=ilHx(5);

        %Repeated - W_INDEPENDENT EVALUATIONS - Steps 1 to 4
        
        %Step 1: Calculate firing rate and ratio of neurons firing (to the
        %total number in the cortex) for principal cortical neurons
        if inptCS1==1, AcxpCS1=1; rNcxpCS1=1; else AcxpCS1=0;rNcxpCS1=0;end
        if inptA==1, AcxpA=1; rNcxpA=1; else AcxpA=0;rNcxpA=0;end
        if inptB==1, AcxpB=1; rNcxpB=1; else AcxpB=0;rNcxpB=0;end
        AcxpA1=AcxpA*AcxpCS1;rNcxpA1=AcxpA1; 
        AcxpB1=AcxpB*AcxpCS1;rNcxpB1=AcxpB1; 

        %Step 2: Calculate firing rate and ratio of neurons firing (to the
        %total number in the hippocampus) for hippocampal neurons
        AhcA=ilA*(1-conA)*(1-Hx)*(1-Hs*(pcentHs/100))*AcxpA*rNcxpA; 
        rNhcA=AhcA;
        AhcB=ilB*(1-conB)*(1-Hx)*(1-Hs*(pcentHs/100))*AcxpB*rNcxpB; 
        rNhcB=AhcB;
  
        AhcA1=ilA1*(1-conA1)*(1-Hx)*(1-Hs*(pcentHs/100))*sqrt(AcxpA*AcxpCS1)*sqrt(rNcxpA*rNcxpCS1); 
        rNhcA1=AhcA1;
        AhcB1=ilB1*(1-conB1)*(1-Hx)*(1-Hs*(pcentHs/100))*sqrt(AcxpB*AcxpCS1)*sqrt(rNcxpB*rNcxpCS1); 
        rNhcB1=AhcB1;

        %Step 3: Calculate firing rate and ratio of cortical inhibitory neurons (to
        %total number of cortical neurons)
        AcxiA=ilA*((1-conA)*AhcA + linsig(conA+ilHxA,0,1)*AcxpA)*((1-conA)*rNhcA + conA*rNcxpA); 
        rNcxiA=AcxiA;
        AcxiB=ilB*((1-conB)*AhcB + linsig(conB+ilHxB,0,1)*AcxpB)*((1-conB)*rNhcB + conB*rNcxpB); 
        rNcxiB=AcxiB;
  
        %Step 4: Caculate firing rate and ratio of neurons firing for transformed
        %cortical neurons
        AcxtA=linsig(conA+ilHxA,0,1)*AcxpA*rNcxpA; 
        rNcxtA=AcxtA;
        AcxtB=linsig(conB+ilHxB,0,1)*AcxpB*rNcxpB; 
        rNcxtB=AcxtB;
        AcxtA1=linsig(conA1+ilHxA1,0,1)*sqrt(AcxpA*AcxpCS1)*sqrt(rNcxpA*rNcxpCS1);
        rNcxtA1=AcxtA1;
        AcxtB1=linsig(conB1+ilHxB1,0,1)*sqrt(AcxpB*AcxpCS1)*sqrt(rNcxpB*rNcxpCS1);
        rNcxtB1=AcxtB1;

        if PFCs==1
            AcxtA=0;rNcxtA=0;
            AcxtB=0;rNcxtB=0;
            AcxtA1=0;rNcxtA1=0;
            AcxtB1=0;rNcxtB1=0;
        end
        
        %Array of ratio of neurons active for different experiment
        %scenarios 
        rNcxp=[rNcxpA rNcxpB rNcxpCS1 rNcxpA1 rNcxpB1 0 0 0];
        Ncxp=rNcxp.*NP; 
        rNhc=[rNhcA rNhcB 0 rNhcA1 rNhcB1 0 0 0];
        Nhc=rNhc.*NH;
        rNcxi=[rNcxiA rNcxiB 0 0 0 0 0 0];
        Ncxi=rNcxi.*NI;
        rNcxt=[rNcxtA rNcxtB 0 rNcxtA1 rNcxtB1 0 0 0];
        Ncxt=rNcxt.*NT;
        
        NLAactive=[Ncxp Ncxi Ncxt]; %Number of neurons active in the LA is based on the number of neurons active in the cortex
        NBLactive=[Nhc Np]; % Number of neurons active in the BL is based on the number of neurons acrtive in the hippocampus and the number of neurons in the cortex

        ActLA=(1-LAs)*[AcxpA AcxpB AcxpCS1 AcxpA1 AcxpB1 0 0 0 ...   
                       AcxiA AcxiB   0       0      0    0 0 0 ...
                       AcxtA AcxtB   0     AcxtA1 AcxtB1 0 0 0];   
        %Loop through and check if the number of neurons active in the LA
        %in 
        for k=1:3*8
            if NLAactive(k)==0,ActLA(k)=0;end
        end

        coefficientLA=ActLA.*NLAactive.*HabitLA;
        
    end %end of if recalcActs == 1
    
    % Step 5: Evaluating the Amygdala
    
    %Step 5A: 
    % Lateral Amygdala Evaluations
    Wp=dot(coefficientLA,Gpe); %Gpe specifies the conductance produced in LA principal cells, numerator of Vpdist formula
    Ui=dot(coefficientLA,Gie); %Gie specifies the conductance produced in LA inhibitory cells
    Vi=Ui*E/(1+Ui); %Calculate depolarzation for inhibitory neurons
    Ai=linsig(Vi,0,ViMax);  % Calculate activity for inhbitory interneurons.
    VpDist=Wp*E/(1+Wp);   %Calculate V for distal compartment of principal cell.
    Vp = Dp*VpDist/(1+Ai*Gi); %Get V for proximal compartent
    Ap=linsig(Vp,0,VpMax); %Get Activations for proximal compartent
    
    if LAs==1, Vp=0; Vi=0;Ai=0;Ap=0; end  %Effect LA suppression
    
    %Basal Amygdala Evaluations
    ActBL=[AhcA  AhcB  0  AhcA1  AhcB1  0 0 0 Ap]; %Activation profile of Basal Amygdala
    %Note that Ap has been tacked on to the end of the above array and to BL cell
    %synaptic weight arrays because the LA principal cells provide one of the inputs to BL cells.
    %However, in this version of the program the LA-to-BL synapses are not plastic.
    
    %Set activation to zero if N=0, precaution for eligibility.
    for k=1:8  
      if NBLactive(k)==0,ActBL(k)=0;end
    end  
    coefficientBL=ActBL.*NBLactive.*HabitBL;
    
    Wq=dot(coefficientBL,Gqe); Uj=dot(coefficientBL,Gje);
    Vj=Uj*E/(1+Uj); Aj=linsig(Vj,0,VjMax);       
    VqDist=Wq*E/(1+Wq);                           
    Vq=Dq*VqDist/(1+Aj*Gj); Aq=linsig(Vq,0,VqMax); 
    if BLs==1, Vq=0;Vj=0;Aj=0;Aq=0; end
    
    %Step 5B: If the BL is supressed then CEm cells are driven by 
    % potentiated (post-learning/consolidation) LA principal neurons
    if BLs==0     
      Acem=Aq; %Medial central amygdala neurons are driven by principal neurons in basal amygdala
    else
      Acem=linsig(Ap*Gp_cem,0,1); 
    end 
    if CEMs==1, Acem=0; end
    
    %Step 6: Calcuclate Af and Freezing scores 
    %Note: Freezing occurs as a result of a conditioned stimulus and the
    %central medial amygdala activation profile (as CE is modeled as
    %"output" 
    if us==0, Af=linsig(Acem,a,b);else Af=0; end
    %Given an undconditioned stimulus (e.g. a shock not paired with a tone)
    %the mouse will physically react to the jolt by hopping around.
    
    Afhist(interval)=Af; %Array to keep track of freezing history
    
    %Fsmrd is a temporally smeared version of Af to mimic fading of memory
    %over time
    if interval==1
      Fsmrd0(interval)=0; %Initialize to zero
    else
      if Afhist(interval)> Fsmrd0(interval-1)
        Fsmrd0(interval)=Fsmrd0(interval-1)+c*(Afhist(interval)-Fsmrd0(interval-1)); %learning has occured
      else
        Fsmrd0(interval)=Fsmrd0(interval-1) - d*(Fsmrd0(interval-1) - Afhist(interval)); %forgetting has occured
      end  
    end
    %To calculate the smeared learning value, you add the previous smeared
    %learning value to the change in learning (some constant multiplied by
    %the difference between the freezing score on this interval and the
    %smeared freezing score from the last interval)
    
    %It is Fsmrd that we care about because that is what has been retained
    %long term
    %If unconditioned stimulus then freezing wouldn't happen
    if us==0, Fsmrd(interval)=Fsmrd0(interval); else Fsmrd(interval)=0;end
    
    %The activation profiles of the LA and BL to be used for figures
    %Activation calculated from depolarization values
    netErecipI(interval,1)=VpDist/E; %unitless method to quantify the excitatory input received by LA principal cells
    netErecipI(interval,2)= 1+Ai*Gi;
    netErecipI(interval,3)=VqDist/E;
    netErecipI(interval,4)=1+Aj*Gj;
    Net_Actdata(interval,1:4)=[Vp Ap Vq Aq];       
    
    %Step 7: Neuromodulation - Calculate the firing rate of R neurons and X neurons
    
    %Used from previous computation interval: 
    %oldRinput used to compute change (uses previous naloxone value)
    %AopPrv used to discount Au for prior CR
    Aop=linsig(Acem,opThrsh,1);
    As=linsig(Acem,0,maxsec);
    Au=us*(Ishk/100)*(1-((1-naloxone)*AopPrv)^pi); %naxolone = opiate receptor blocking agent
    newRinput=Au+psec*As;
    Ar=linsig(newRinput-oldRinput,0,1); %Ar will be zero if no change in CEm activity
    
    if (Ar>0)||(us==1),Ax=0;else Ax=(1-naloxone)*Aop;end    
    %If freezing occurred as a result of a conditioned stimulus, then X
    %neuron will not have an effect
    %If unconditioned stimulus then synaptic weights will decrease in BLA,
    %via X neuron effect
    %R neuron - reinforcement = neuromodulation increase synaptic weights
    % for excitatory input to principal and inihibitory neurons of BLA
    %X neuron = neuromodulation decreases synaptic weights for excitatory
    % input to principal and inhibitory neurons of BLA
    
    %Step 8: Calculate principal cell backpropagation activity....
    % and inhibitory interneuron Ca levels.
    
    %Backpropagation = graded spike-like signal, based on
    % (1) The depolarizationof proximal compartment of principal cells
    % (2) Direct (non-neuromodulatory) input from R 
    Bp=100*linsig(Vp+100*Ar*ArMultLA,0,100); %Backpropagation of principal cells in LA
    Bq=100*linsig(Vq+100*Ar*ArMultBL,0,100); %Backpropagation of principal cells in BL
    
    %Calcium levels effect on activation of inhibitory interneurons
    %For inhibitory interneurons, CA level determined by: 
    % (1) Recurrent input from principal cells
    % (2) Neuromodulator input from R neuron
    kappaLA=100;kappaBL=100;
    CaI=100*linsig(Ap*VpMax+ kappaLA*Ar, 0, 100); 
    CaJ=100*linsig(Aq*VqMax+ kappaBL*Ar, 0, 100);
    if LAs==1, Bp=0; CaI=0; end
    if BLs==1, Bq=0; CaJ=0; end
    
    %Step 9: Update amygdala synaptic weights...
    %Lateral Amygdala:
    % eligibleLA --> 0 or 1, neurons in absolute refractory period
    % shouldn't fire
    %To determine the change in synaptic strength of synapses on amygdala
    % principal cells, multiply activation of R neuron by the factor that is
    % a function of postsynaptic depolarization (ie backpropagation)
    %To determine the change in synaptic strength of synapses on amygdala
    %inhibitory neurons, r and x factors are a function of calcium levels 
    for k=1:8*3 
      if Ar>0 %if reinforcing connectivity, then synaptic weights increase
        deltaGpe= eligibleLA(k)*Ar*alphaLA*linsig(Bp,thrshRp,E);
        deltaGie=-eligibleLA(k)*Ar*etaLA*linsig(CaI,thrshRi,E);
      else
        deltaGie= eligibleLA(k)*Ax*betaLA*BumpCa(CaI,thrshXi,E,gammaXi);
        deltaGpe=-eligibleLA(k)*Ax*zetaLA*BumpB(Bp,thrshXp,endXp);        
      end  
      Gpe(k)=Gpe(k)+deltaGpe; 
      Gie(k)=Gie(k)+deltaGie;   
      if Gpe(k)<0, Gpe(k)=0;end
      if Gie(k)<0, Gie(k)=0;end
    end
    %Basal Amygdala
    for k=1:8
      if Ar>0 
        deltaGqe=eligibleBL(k)*Ar*alphaBL*linsig(Bq,thrshRq,E);
        deltaGje=-eligibleBL(k)*Ar*etaBL*linsig(CaJ,thrshRj,E);          
      else
        deltaGje= eligibleBL(k)*Ax*betaBL*BumpCa(CaJ,thrshXj,E,gammaXj);
        deltaGqe=-eligibleBL(k)*Ax*zetaBL*BumpB(Bq,thrshXq,endXq);
      end  
      Gqe(k)=Gqe(k)+deltaGqe; 
      Gje(k)=Gje(k)+deltaGje;       
      if Gqe(k)<0, Gqe(k)=0;end
      if Gje(k)<0, Gje(k)=0;end
    end  
    Gqe(9)=Gpq; Gje(9)=0; %synaptic connectivity between LA and BL not plastic
    
    %The activation profiles of the LA and BL to be used for figures
    netErecipI(interval,1)=VpDist/E;
    netErecipI(interval,2)= 1+Ai*Gi;
    netErecipI(interval,3)=VqDist/E;
    netErecipI(interval,4)=1+Aj*Gj;
    Net_Actdata(interval,1:4)=[Vp Ap Vq Aq]; 

    %Color Order: red grn blk purple cyan  
    %Variables to be graphed later
    Wdata(interval,1)=0;   %plotted red
    Wdata(interval,2)=0;   %plotted green
    Wdata(interval,3)=0;   %plotted black
    Wdata(interval,4)=0;   %plotted purple
    Wdata(interval,5)=0;   %plotted cyan 
    Wdata(interval,6)=0;   %plotted gray

    %Step 10: Re-evaluating the Amygdala
   
    % Lateral Amygdala Evaluations
    Wp=dot(coefficientLA,Gpe); %Gpe specifies the conductance produced in LA principal cells, numerator of Vpdist formula
    Ui=dot(coefficientLA,Gie); %Gie specifies the conductance produced in LA inhibitory cells
    Vi=Ui*E/(1+Ui); %Calculate depolarzation for inhibitory neurons
    Ai=linsig(Vi,0,ViMax);  % Calculate activity for inhbitory interneurons.
    VpDist=Wp*E/(1+Wp);   %Calculate V for distal compartment of principal cell.
    Vp = Dp*VpDist/(1+Ai*Gi); %Get V for proximal compartent
    Ap=linsig(Vp,0,VpMax); %Get Activations for proximal compartent
    
    if LAs==1, Vp=0; Vi=0;Ai=0;Ap=0; end  %Effect LA suppression
    
    %Basal Amygdala Evaluations
    %Activation profile of Basal Amygdala
    ActBL=(1-BLs)*[AhcA  AhcB 0 AhcA1 AhcB1 0 0 0 Ap];
    
    %ActBL=[AhcA  AhcB  0  AhcA1  AhcB1  0 0 0 Ap];
    %Set activation to zero if N=0, precaution for eligibility.
    for k=1:8  
      if NBLactive(k)==0,ActBL(k)=0;end
    end  
    coefficientBL=ActBL.*NBLactive.*HabitBL;
    Wq=dot(coefficientBL,Gqe); Uj=dot(coefficientBL,Gje);
    Vj=Uj*E/(1+Uj); Aj=linsig(Vj,0,VjMax);       
    VqDist=Wq*E/(1+Wq);                           
    Vq=Dq*VqDist/(1+Aj*Gj); Aq=linsig(Vq,0,VqMax); 
    if BLs==1, Vq=0;Vj=0;Aj=0;Aq=0; end
    
    %Step 11: Update connectivity between lateral amygdala and central
    %amygdala 
    %When the basal amygdala is supressed the connectivyt between LA and
    %CEm determines activity in CEm (this pathway strengthens with
    %strengthened LTP/mimics long term memory influence)
    if BLs==1 && Ap>0 && us==1 && CEMs==0, Gp_cem=Gp_cemON; end
    %IE. Conductance of central amygdala neurons get the values of that of
    %the conductance parameter specifying strength of fully potentiated
    %LAp-CEm connection
    
    %Step 12: Re-evaluate ACEm and  Rinpt with new synaptic weights
    Ar=0;
    if BLs==0
      Acem=Aq;
    else
      Acem=linsig(Ap*Gp_cem,0,1);
    end  
    if CEMs==1, Acem=0; end
    Aop=linsig(Acem,opThrsh,1);
    As=linsig(Acem,0,maxsec);        
    Au=us*(Ishk/100)*(1-((1-naloxone)*AopPrv)^pi);      
    oldRinput=Au+psec*As;
    AopPrv=Aop; 
    
    %Step 13: Update learning factors for hippocampus and cortex
    if Hx==0 && Hs==0      
      if AcxpA>.01 
        recalcActs=1; 
        ilo=ilA;
        if ilA<mx,  ilA=ilA+ ilTau*(1-ilA)^ilExp;  else ilA=1;end
        Gqe(1)=Gqe(1)*(ilo/ilA) ;  
        Gje(1)=Gje(1)*(ilo/ilA) ;
      end  
      if AcxpB>.01
        recalcActs=1;
        ilo=ilB;if ilB<mx,  ilB=ilB+ ilTau*(1-ilB)^ilExp;  else ilB=1;end
        Gqe(2)=Gqe(2)*(ilo/ilB) ;
        Gje(2)=Gje(2)*(ilo/ilB) ;
      end  
      if AcxpA>.01 && AcxpCS1>.01 
        recalcActs=1;
        ilo=ilA1;if ilA1<mx,ilA1=ilA*linsig(ilA1+ilCmpndTau*(1-ilA1)^ilCmpndExp,0,ilA);else ilA1=1;end 
        Gqe(4)=Gqe(4)*(ilo/ilA1) ;
        Gje(4)=Gje(4)*(ilo/ilA1) ;
      end     
      if AcxpB>.01 && AcxpCS1>.01
        recalcActs=1;
        ilo=ilB1;if ilB1<mx, ilB1=ilB*linsig(ilB1+ilCmpndTau*(1-ilB1)^ilCmpndExp,0,ilB);else ilB1=1;end  %Try this
        Gqe(5)=Gqe(5)*(ilo/ilB1) ;
        Gje(5)=Gje(5)*(ilo/ilB1) ;
      end 
      il=[ilA ilB 0 ilA1 ilB1 0 0 0];
      
    elseif Hx==1 && PFCs<1
      if AcxpA>.01  
        recalcActs=1;
        ilHxo=ilHxA;if ilHxA<mxHx,ilHxA=ilHxA+ilHxTau*(1-ilHxA)^ilHxExp;else ilHxA=1;end
        Gpe(17)=Gpe(17)*(ilHxo/ilHxA) ;
        Gie(17)=Gie(17)*(ilHxo/ilHxA) ;
      end  
      if AcxpB>.01  
        recalcActs=1;
        ilHxo=ilHxB;if ilHxB<mxHx,ilHxB=ilHxB+ilHxTau*(1-ilHxB)^ilHxExp;else ilHxB=1;end
        Gpe(18)=Gpe(18)*(ilHxo/ilHxB) ;
        Gie(18)=Gie(18)*(ilHxo/ilHxB) ;
      end
      
      if AcxpA>.01 && AcxpCS1>.01 
        recalcActs=1;
        ilHxo=ilHxA1;if ilHxA1<mxHx, ilHxA1=ilHxA*linsig(ilHxA1+ ilHxCmpndTau*(1-ilHxA1)^ilHxCmpndExp,0,ilHxA);else ilHxA1=1;end  %Try this
        Gpe(20)=Gpe(20)*(ilHxo/ilHxA1) ;
        Gie(20)=Gie(20)*(ilHxo/ilHxA1) ;
      end     
      if AcxpB>.01 && AcxpCS1>.01 
        recalcActs=1;
        ilHxo=ilHxB1;if ilHxB1<mxHx, ilHxB1=ilHxB*linsig(ilHxB1+ilHxCmpndTau*(1-ilHxB1)^ilHxCmpndExp,0,ilHxB);else ilHxB1=1;end  %Try this
        Gpe(21)=Gpe(21)*(ilHxo/ilHxB1) ;
        Gie(21)=Gie(21)*(ilHxo/ilHxB1) ;
      end        
      ilHx=[ilHxA ilHxB 0 ilHxA1 ilHxB1 0 0 0];
    end  % End of "if Hx==0 && Hs==0".      
    
    %Step 14: Update neurons to reflect refractory period (eligibility)
    for k=1:8*3
      eligibleLA(k)=linsig(ActLA(k),thrshElig,1); 
    end  
    for k=1:8;
      eligibleBL(k)=linsig(ActBL(k),thrshElig,1);  
    end  
    if coefficientBL(9)>0, eligibleBL(9)=1;else eligibleBL(9)=0;end;
    
    %Progress graph code continued
    if mod(interval,50)==0  %For Progress graph during computations
      pause(.1) %If don't put a pause in here, graph doesn't get made until end.
      bars=[interval Npoints];
      barh(bars);
    end    
  end  %End of "for i= 1 : dur" (within-line for loop)
end   %End of  "for j=continueAt :Nsched" (lines of ExpData table)

if interval>=Npoints
  hh=findobj('Tag','run_continue');
  set(hh,'BackgroundColor','white');        
end    

%Graphing Code ...??


%Consolidation: ...??

    
%%%%%%%%%%%%%%%%%%%%%%%%%% Functions Called %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=BumpCa(x,ExponentThrsh,ExponentMxAt,phi)
th=ExponentThrsh;
ma=ExponentMxAt;
ph=phi;
y=1-exp(-1*ph*linsig(x,th,ma));
end

function y=BumpB(x,Thrsh, EndAt)
th=Thrsh;
ea=EndAt;
y=0;
if x>th && x<ea, y=1;end
end

function y=linsig(x, start, stop)
    if x<start 
        y=0;
    elseif x>=stop  
        y=1;
    else
        y=(x-start)/(stop-start);
    end 
end
end