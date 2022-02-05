%=================Study and Choose Best Damping Value after 3D Inversion S-wave Tomo ===================
%This requires to do the L-plot:  Var(data) vs Variace[model]
% I do it using weigths, based on Tobias paper "High-resolution 3-D P-wave model of the Alpine crust"  https://core.ac.uk/download/pdf/85209659.pdf
clear all 
close all

%--Parameters--
nEQ=99;  % Nb events 


%===============================================================
%%   VALUES OBTAINED DIRECTLY FROM SIMULPS  "output" FILE
%Better use these ones
%===============================================================
%Avg Picking Error P-phases  ~ ± 0.0215s  (Phases used in Min1D Vmodel)
%Avg Picking Error S-phases  ~ ± 0.0622s
avg_pickerr=0.062;

%Damping values d=1, 10, 100, 1000, 10000   
d=[1 10 100 1000]; %10mil got errors, I need to check
nit=7; %Nb iterations tested


%-----------MOD 2x2x1-------------
var_m=zeros(length(d),nit); %Model Variance  [dampingParam1iteration1 d1it2 ...d1it7; d2it1 d2it2...d2it7; ....  ;d4it1 d4it2...d4it7]
var_d=var_m; %Data Variance

var_m=[0.00441 0.01201  0.01518  0.01812  0.02096  0.02096  0.02096  ;  %d1 stop in 5it
       0.00216  0.00545  0.00675  0.00779   0.00856  0.00920  0.00920 ; %d2 stop in 6it 
       0.00033   0.00098  0.00129  0.00164  0.00201  0.00237  0.00271;  %d3 normal behaviour for all it
       0.00001   0.00003   0.00004  0.00004  0.00004   nan      nan ]; %d2 stop in 3it 

var_d=[ 0.021625  0.009663  0.006329  0.005445  0.005200   0.005197 0.005264;  % d1 start having backups in var_d after 5it
        0.022916 0.011221  0.007702  0.006562  0.006210  0.005884  0.005965;   % d2 start having backups in var_d after 6it
        0.039207  0.021592  0.018290  0.015941  0.014307  0.012874  0.011829 ;  
        0.052026  0.031296  0.030588  0.030322   0.030380  nan      nan ];      % d4 No more data after 3it?
    
%Weighted RMS
w_rms=[ 0.17545  0.10761  0.06240  0.05039  0.04701  0.04572  0.04568 ;  
        0.17545  0.11113  0.06745  0.05635  0.05239  0.05094  0.04946 ;  
        0.17545  0.14645  0.09753  0.08994  0.08392  0.07928   0.07516;
        0.17545  0.17087  0.12112  0.11975  0.11928   nan        nan]; %After 5it I see backups

figure(1)
for i=1:length(d)
 plot(var_m(i,:) , var_d(i,:),'-o','LineWidth',3)  
     hold on
end
     title('MOD 2x2x1')
    ylabel('Var[data] (s^2)','fontsize',16)
    xlabel('Var[model] (km/s)^2','fontsize',16)
    legend('d=1','d=10','d=100','d=1000') 
   set(gca,'fontsize',25)
   grid on
      % xlim([0 16]) 
  
      figure(2)
      for i=1:length(d)
 plot(w_rms(i,:),'-o','LineWidth',3)  
     hold on
      end
    %  plot(avg_pickerr*ones(1,nit),'k--','LineWidth',1)
title('MOD 2x2x1')
    ylabel('Weighted RMS (s)','fontsize',16)
    xlabel('iteration','fontsize',16)
    legend('d=1','d=10','d=100','d=1000') 
   set(gca,'fontsize',25)
   grid on
    xlim([1 nit]) 
    set(gca,'XTick',[0:1:nit]) 
    
    
    
%-----------MOD 3x3x1.5-------------
%Damping values d=1, 10, 100, 1000, 10000   

var_m=zeros(length(d),nit); %Model Variance
var_d=var_m; %Data Variance

var_m=[0.00423  0.01272  0.01665  0.02049 0.02423  0.02794   0.02794  ;
        0.00244  0.00652  0.00788 0.00879  0.00947  0.01001  0.01001 ;
        0.00058  0.00167  0.00217  0.00271  0.00322  0.00369  0.00412 ;
        0.00002  0.00007  0.00009  0.00013  0.00016  0.00021  0.00021 ];



var_d=[ 0.018821  0.008297  0.005713  0.005188  0.004960  0.005003  0.004867;
        0.017317  0.008345  0.006144  0.005764   0.005473 0.005369  0.005334 ;
        0.026455  0.013322  0.010863   0.009224  0.008252 0.007576  0.007213 ;
        0.039480  0.021818  0.020866  0.020111   0.019291  0.018657  0.018735 ];


%Weighted RMS
w_rms=[ 0.17848 0.11835  0.07136  0.05916 0.05643  0.05518    nan; %backups after 6it
        0.17848  0.11336 0.07144  0.06151  0.05960  0.05808  0.05756   ;
        0.17848  0.14040  0.09127  0.08253  0.07615  0.07193  0.06885   ;
        0.17848  0.17210   0.11859  0.11586  0.11374   0.11129  0.10941 ];

figure(3)
for i=1:length(d)
 plot(var_m(i,:) , var_d(i,:),'-o','LineWidth',3)  
     hold on
end
     title('MOD 3x3x1.5')
    ylabel('Var[data] (s^2)','fontsize',16)
    xlabel('Var[model] (km/s)^2','fontsize',16)
    legend('d=1','d=10','d=100','d=1000') 
   set(gca,'fontsize',25)
   grid on
    %set(gca,'YScale','log')
      % xlim([0 16]) 
  
      figure(4)
      for i=1:length(d)
 plot(w_rms(i,:),'-o','LineWidth',3)  
     hold on
      end
  %    plot(avg_pickerr*ones(1,nit),'k--','LineWidth',1)
title('MOD 3x3x1.5')
    ylabel('Weighted RMS (s)','fontsize',16)
    xlabel('iteration','fontsize',16)
    legend('d=1','d=10','d=100','d=1000') 
   set(gca,'fontsize',25)
   grid on
    xlim([1 nit]) 
    set(gca,'XTick',[0:1:nit]) 
      % set(gca,'YScale','log')
      
      %-----------MOD 5x5x3-------------
%Damping values d=1, 10, 100, 1000, 10000   
var_m=zeros(length(d),nit); %Model Variance
var_d=var_m; %Data Variance

var_m=[0.00559 0.01545  0.02284  0.03174  0.03974  0.03974  0.03974   ;
       0.00314  0.00919  0.01202  0.01400   0.01530  0.01530   nan      ;
        0.00131 0.00363  0.00464   0.00566  0.00652   0.00726 0.00791  ;
        0.00011  0.00035  0.00046  0.00060  0.00074  0.00089  0.00105];



var_d=[0.019305  0.008699   0.005916  0.005274  0.005172  0.005150  0.005102  ;
       0.018157  0.008229   0.006126  0.005775  0.005662  0.005590  nan    ;
        0.021011 0.010892   0.008740  0.007804  0.007208  0.006936  0.006725        ;
        0.033471  0.017896  0.016809  0.015924  0.015192  0.014357  0.013825 ];


%Weighted RMS
w_rms=[ 0.18768  0.13061  0.08104   0.06677  0.06299  0.06240  0.06223   ;
        0.18768   0.12667  0.07888   0.06804  0.06804   0.06538 0.06498  ;
        0.18768   0.13632   0.09083   0.08134   0.07687  0.07383  0.07242  ;
        0.17211   0.11669    0.11306   0.11001   0.10743  0.10443  0.10248   ];

figure(5)
for i=1:length(d)
 plot(var_m(i,:) , var_d(i,:),'-o','LineWidth',3)  
     hold on
end
     title('MOD 5x5x3')
    ylabel('Var[data] (s^2)','fontsize',16)
    xlabel('Var[model] (km/s)^2','fontsize',16)
    legend('d=1','d=10','d=100','d=1000') 
   set(gca,'fontsize',25)
   grid on
    %set(gca,'YScale','log')
      % xlim([0 16]) 
  
      figure(6)
      for i=1:length(d)
 plot(w_rms(i,:),'-o','LineWidth',3)  
     hold on
      end
   %   plot(avg_pickerr*ones(1,nit),'k--','LineWidth',1)
title('MOD 5x5x3')
    ylabel('Weighted RMS (s)','fontsize',16)
    xlabel('iteration','fontsize',16)
    legend('d=1','d=10','d=100','d=1000') 
   set(gca,'fontsize',25)
   grid on
    xlim([1 nit]) 
    set(gca,'XTick',[0:1:nit]) 
      % set(gca,'YScale','log')
  