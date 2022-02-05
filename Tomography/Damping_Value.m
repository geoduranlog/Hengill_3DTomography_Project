%=================Study and Choose Best Damping Value after 3D Inversion Tomo ===================
%This requires to do the L-plot:  Var(data) vs Variace[model]
% I do it using weigths, based on Tobias paper "High-resolution 3-D P-wave model of the Alpine crust"  https://core.ac.uk/download/pdf/85209659.pdf
clear all 
close all

%--Parameters--
nEQ=130;  % Nb events 


%%  MODEL VARIANCE  - Computed by myself. But Better consider values from SIMULPS output file!! (see bellow)

%--- SELECT INPUT FILES---
%Select label of tests you want to study
%----------------------------------
name='Coseismiq_3x3x1p5_d100';  %  %Coseismiq_2x2x1_d100
folder=['/Users/alejandro/ICELAND_Project/Tomography/SIMULPS14/TEST_simulps/',name,'/'];

%{
%--Ouput Vel Model after inversion (SimulPS)---
FILE=[folder,'velomod.out']; 


% Load 
fileid1 = fopen(FILE,'rt');
%  Nb columns is the nb nodes in x-direction (LON points)
%C=textscan(fileid1,' %f %f %f %f %f    %f %f %f %f %f    %f %f %f %f %f   %f %f ','headerLines',5); %MOD2x2x1
C=textscan(fileid1,' %f %f %f %f %f    %f %f %f %f %f    %f %f %f %f %f  ','headerLines',5); %MOD3x3x1p5
fclose(fileid1);


%---Total Output Vmodel values - not in order--
vel=cell2mat(C);   %From cell to matrix 
vel = vel'; vel = vel(:)'; % Tranform all into array
nnodes=length(vel);  %Nb of nodes
mean_v=mean(vel);

%----Variance of Model (km/s)^2----
var_mod=var(vel)
%----------------------

%a=(vel-mean(vel)).^2;
%var=sum(a)/(length(a)-1)



%%  DATA VARIANCE
%Output file from SIMULPS is "residuals",  But in order to read it easily
%in MATLAB I create a new file "residuals_label2matlab" in which I put a
%label "#" to comment the lines that I don't need

%--Ouput Residuals (and weigth classes) after inversion (SimulPS). Modified file labeled with # ---
FILE=[folder,'residuals_label2matlab']; 

% Load 
fileid1 = fopen(FILE,'rt');
R=textscan(fileid1,' %s %f %f %f %f    %f %f %f %f %f    %f %f ','headerLines',2,'TreatAsEmpty',{'NA','na'},'CommentStyle','#');
fclose(fileid1);


%Time Residuals per each observation (s)
tres=R{3}; % => Nb of phases=length(tres)

%Weigths (Picking Classes)
for i=1:length(R{1})
w(i)=str2double( R{1}{i}(end) );

end

%Replace with real weigths  Classes=[0 1 2 3 4] -> Weigths=1/(2^classes), i.e. [1 0.5 0.25 0.125 0.625] 
%Better start from the lower weigth (replace 0 by 1 only at the end!!)
w(w==4)=0.0625; w(w==3)=0.125; w(w==2)=0.25;
w(w==1)=0.5; w(w==0)=1;
 
%---Weighted number of degrees of freedom  DOBLE CHECK THIS...
nnodes_inv=529;  %Vnodes that were actually inverted (Check in file "output"  last it where it says "nodes inverted for=") 
wndf=sum(w)-4*nEQ-nnodes_inv;



%---Variance of data (s^2) -----
temp=w'.*(tres.^2);
var_data=sum(temp)./wndf
%----------------------


%----- Weithed RMS (s) -------
rms_w=sqrt(  sum(temp)/sum(w) )
%}

%===============================================================
%%   VALUES OBTAINED DIRECTLY FROM SIMULPS  "output" FILE
%Better use these ones
%===============================================================
%Avg Picking Error P-phases  ~ ± 0.0215s  (Phases used in Min1D Vmodel)
%Avg Picking Error S-phases  ~ ± 0.0622s
avg_pickerr=0.0215;

%-----------MOD 2x2x1-------------
%Damping values d=1, 10, 100, 1000, 10000   
d=[1 10 100]; %mil y 10mil got errors, I need to check
nit=7; %Nb iterations tested

var_m=zeros(length(d),nit); %Model Variance  [dampingParam1iteration1 d1it2 ...d1it7; d2it1 d2it2...d2it7; ....  ;d4it1 d4it2...d4it7]
var_d=var_m; %Data Variance

var_m=[ 0.00558	0.01737	0.02399	0.02932	0.03333	0.03682	0.04021; 0.0023	0.00684	0.00913	0.01121	0.01293	0.01442	0.01569; 0.00014	0.0005	0.00077	0.00108	0.00144	0.00181	0.0022];
var_d=[0.003317	0.00141	0.000808	0.000656	0.000603	0.00058	0.000564; 0.003763	0.001736	0.001249	0.001036	0.000924	0.000858	0.000806; 0.005937	0.003315	0.002975	0.002701	0.002483	0.002288	0.002132];

%Weighted RMS
w_rms=[0.06865	0.04718	0.02821	0.02134	0.01928	0.0185	0.01814; 0.06865	0.05034	0.03137	0.02659	0.02422	0.02287	0.02204; 0.06865	0.0635	0.04394	0.04163	0.03971	0.03808	0.03658 ]

%{
figure(1)
for i=1:length(d)
 plot(var_m(i,:) , var_d(i,:),'-o','LineWidth',3)  
     hold on
end
     title('MOD 2x2x1')
    ylabel('Var[data] (s^2)','fontsize',16)
    xlabel('Var[model] (km/s)^2','fontsize',16)
    legend('d=1','d=10','d=100') 
   set(gca,'fontsize',25)
   grid on
      % xlim([0 16]) 
  
      figure(2)
      for i=1:length(d)
 plot(w_rms(i,:),'-o','LineWidth',3)  
     hold on
      end
      plot(avg_pickerr*ones(1,nit),'k--','LineWidth',1)
title('MOD 2x2x1')
    ylabel('Weighted RMS (s)','fontsize',16)
    xlabel('iteration','fontsize',16)
    legend('d=1','d=10','d=100') 
   set(gca,'fontsize',25)
   grid on
    xlim([1 nit]) 
    set(gca,'XTick',[0:1:nit]) 
  %}  
    
    
% %-----------MOD 3x3x1.5-------------
% %Damping values d=1, 10, 100, 1000, 10000   
% d=[1 10 100 1000]; %mil y 10mil got errors, I need to check
% nit=7; %Nb iterations tested
% 
% var_m=zeros(length(d),nit); %Model Variance  [dampingParam1iteration1 d1it2 ...d1it7; d2it1 d2it2...d2it7; ....  ;d4it1 d4it2...d4it7]
% var_d=var_m; %Data Variance
% 
% var_m=[0.00481	0.01621	0.02398	0.03214	0.03785	0.04229	0.04672;
% 0.0028	0.00848	0.01163	0.01438	0.01654	0.01819	0.01953;
% 0.00038	0.00122	0.00176	0.00238	0.00302	0.00368	0.00432;
% 0.00001	0.00002	0.00004	0.00005	0.00008	0.00008	0.00008];
% 
% 
% 
% var_d=[0.003475	0.001504	0.000879	0.000721	0.00068	0.000661	0.000647;
% 0.003326	0.001507	0.001064	0.000907	0.000836	0.000799	0.000775;
% 0.004919	0.002555	0.002226	0.001984	0.001793	0.001647	0.00153;
% 0.006112	0.003603	0.003527	0.003441	0.003378	0.003384	0.003382];
% 
% 
% %Weighted RMS
% w_rms=[0.07207	0.05335	0.03277	0.02506	0.02268	0.02205	0.02173;
% 0.07207	0.05225	0.03283	0.02758	0.02546	0.02445	0.02389;
% 0.07207	0.0636	0.04295	0.04009	0.03784	0.03596	0.03446;
% 0.07207	0.07095	0.05107	0.05054	0.04993	0.04947	0.04951];
% 
% figure(3)
% for i=1:length(d)
%  plot(var_m(i,:) , var_d(i,:),'-o','LineWidth',3)  
%      hold on
% end
%      title('MOD 3x3x1.5')
%     ylabel('Var[data] (s^2)','fontsize',16)
%     xlabel('Var[model] (km/s)^2','fontsize',16)
%     legend('d=1','d=10','d=100','d=1000') 
%    set(gca,'fontsize',25)
%    grid on
%     %set(gca,'YScale','log')
%       % xlim([0 16]) 
%   
%       figure(4)
%       for i=1:length(d)
%  plot(w_rms(i,:),'-o','LineWidth',3)  
%      hold on
%       end
%       plot(avg_pickerr*ones(1,nit),'k--','LineWidth',1)
% title('MOD 3x3x1.5')
%     ylabel('Weighted RMS (s)','fontsize',16)
%     xlabel('iteration','fontsize',16)
%     legend('d=1','d=10','d=100','d=1000') 
%    set(gca,'fontsize',25)
%    grid on
%     xlim([1 nit]) 
%     set(gca,'XTick',[0:1:nit]) 
%       % set(gca,'YScale','log')
  
      
 % %-----------MOD 3x3x2-------------
%Damping values d=1, 10, 100, 1000, 10000   
%d=[1 8 10 100 1000]; %mil y 10mil got errors, I need to check
d=[1 8 100 1000]; %mil y 10mil got errors, I need to check
nit=8; %Nb iterations tested

var_m=zeros(length(d),nit); %Model Variance  [dampingParam1iteration1 d1it2 ...d1it7; d2it1 d2it2...d2it7; ....  ;d4it1 d4it2...d4it7]
var_d=var_m; %Data Variance

var_m=[0.00  0.00480  0.01620   0.02203   0.02765  0.03208 0.03589 0.03589 ;  %0.006966 it=0 varw.(ssqrw/wndof) = 0.005856
       0.00     0.00279    0.00875  0.01195  0.01454   0.01631  0.01762  0.01873 ;  %d=8
       %0.00  0.00258   0.00799  0.01089    0.01326   0.01498   0.01629  0.01737     ; %d=10
       0.00  0.00039   0.00124  0.00177    0.00235   0.00295   0.00353   0.00410      ;
       0.00  0.00001   0.00003  0.00004    0.00006   0.00006   0.00006   NaN      ];



var_d=[0.004961    0.003194  0.001386   0.000826   0.000692  0.000659    0.000646   0.000649   ;
       0.004961     0.003062 0.001390  0.000927     0.000803 0.000756    0.000739   0.000730  ; %d=8
     %  0.004961        0.003083  0.001419    0.000957   0.000832    0.000781   0.000753  0.000736     ; 
       0.004961     0.004493     0.002232     0.001932  0.001707    0.001556   0.001439   0.001351          ;
       0.004961       0.0036   0.003184     0.003091   0.003033    0.003031   0.003034    NaN        ];

%Weighted RMS
w_rms=[0.07043    0.05186   0.03196   0.02472     0.02262   0.02207     0.02188    0.02191;
       0.07043   0.05080    0.03205   0.02622      0.02439   0.02368   0.02341     0.02326;    %d=8 
     %  0.07043    0.05097    0.03238    0.02664     0.02484   0.02407   0.02364   0.02338     ; 
       0.07043    0.06163    0.04078    0.03792   0.03566    0.03404   0.03272   0.03170      ;
       0.07043    0.07043    0.06928    0.04878     0.04808    0.04762   0.04762    0.04764  ];
     

   
   
figure(3)
for i=1:length(d)
 plot(var_m(i,:) , var_d(i,:),'-o','LineWidth',3)  
     hold on
end
     title('MOD 3x3x2')
    ylabel('Var[data] (s^2)','fontsize',16)
    xlabel('Var[model] (km/s)^2','fontsize',16)
    legend('d=1','d=8','d=100','d=1000') 
   set(gca,'fontsize',25)
   grid on
    %set(gca,'YScale','log')
      % xlim([0 16]) 
  
      figure(6)
  for i=1:length(d)
 %plot(w_rms(i,:),'-o','LineWidth',3)  
 plot([0:nit-1],w_rms(i,:),'-o','LineWidth',3)  
     hold on
      end
      plot([0:nit-1],avg_pickerr*ones(1,nit),'k--','LineWidth',1)
title('MOD 3x3x2')
    ylabel('Weighted RMS (s)','fontsize',16)
    xlabel('iteration','fontsize',16)
    legend('d=1','d=8','d=100','d=1000') 
   set(gca,'fontsize',25)
   grid on
     xlim([0 nit-1]) 
    set(gca,'XTick',[0:1:nit-1])
    ylim([0.015 0.08]) 
      % set(gca,'YScale','log')



      
      %-----------MOD 5x5x3-------------
%Damping values d=1, 10, 100, 1000, 10000   
d=[1 10 100 1000]; %mil y 10mil got errors, I need to check
nit=8; %Nb iterations tested starting from it=0

var_m=zeros(length(d),nit); %Model Variance  [dampingParam1iteration1 d1it2 ...d1it7; d2it1 d2it2...d2it7; ....  ;d4it1 d4it2...d4it7]
var_d=var_m; %Data Variance

var_m=[ 0 0.00608   0.01946  0.03173 0.04657 0.05966    0.07180 0.08255  ; %varw.(ssqrw/wndof) = 0.005476 =0.003405
        0 0.00351   0.01133  0.01695  0.02264  0.02679  0.02996 0.03258     ; 
        0 0.00104   0.00312  0.00437  0.00568  0.00692   0.00810  0.00918    ;
        0 0.00003   0.00013  0.00019  0.00027  0.00037   0.00048  0.00060   ];



var_d=[0.005012  0.003405 0.001833 0.001283 0.001063 0.000963  0.000915  0.000882  ;
       0.005012  0.003190 0.001734  0.001267   0.001114  0.001064  0.001042  0.001026     ;
       0.005012  0.003767  0.002152 0.001858   0.001691  0.001579  0.001500  0.001444   ;
       0.005012  0.003190  0.003033  0.002928   0.002837  0.002755  0.002670  0.002608    ]; %0.005126 


%Weighted RMS
w_rms=[0.07079   0.05579  0.03853 0.03228 0.02937  0.02796  0.02725  0.02675  ;
       0.07079  0.05398   0.03752  0.03207   0.03006 0.02938  0.02908 0.02885;
       0.07079  0.05870   0.04182   0.03885  0.03706   0.03581  0.03490  0.03424 ;
       0.07079  0.06848   0.04968   0.04881  0.04805   0.04735   0.04661  0.04607       ];
     
figure(7)
for i=1:length(d)
 plot(var_m(i,:) , var_d(i,:),'-o','LineWidth',3)  
     hold on
end
     title('MOD  5x5x3')
    ylabel('Var[data] (s^2)','fontsize',16)
    xlabel('Var[model] (km/s)^2','fontsize',16)
    legend('d=1','d=10','d=100','d=1000') 
   set(gca,'fontsize',25)
   grid on
    %set(gca,'YScale','log')
      % xlim([0 16]) 
  
      figure(8)
      for i=1:length(d)
 %plot(w_rms(i,:),'-o','LineWidth',3)  
 plot([0:nit-1],w_rms(i,:),'-o','LineWidth',3)  
     hold on
      end
      plot([0:nit-1],avg_pickerr*ones(1,nit),'k--','LineWidth',1)
title('MOD  5x5x3')
    ylabel('Weighted RMS (s)','fontsize',16)
    xlabel('iteration','fontsize',16)
    legend('d=1','d=10','d=100','d=1000') 
   set(gca,'fontsize',25)
   grid on
    xlim([0 nit-1]) 
    set(gca,'XTick',[0:1:nit-1])
    ylim([0.015 0.08]) 
      % set(gca,'YScale','log')
  