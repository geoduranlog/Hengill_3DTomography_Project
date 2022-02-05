%==============Travel time Residual for EQ relocation using various Min 1D ===================
%Tobias used different Vel Model to relocate EQ and test which Vmodel gave the best (lowest) RMS residual

clear all
close all
clc
   
%---------------------------------------    
% File 1 - Relocation with HU-P2 ISOR model
  file='Relocation_HUP2_EQ20200617_02.55.25.csv';
 [Ph1 NET1 Res1 Dist1] = csvimport( file, 'columns', {'PHASE', 'NET', 'TimeRes(s)', 'Dist(km)'});
 Res1=str2double(Res1); %Convert string to double
 
%--Phase selection ---
idx_P1=find(contains(Ph1,'Pg'));   %Test itPh(idx_P)
idx_S1=find(contains(Ph1,'Sg')); 


%---------------------------------------    
% File 2 - Relocation with P and S are the individual derived models using EDT
  file='Relocation_PandSindividualModels_EQ20200617_02.55.25.csv';
 [Ph2 NET2 Res2 Dist2] = csvimport( file, 'columns', {'PHASE', 'NET', 'TimeRes(s)', 'Dist(km)'});
 %Res2=str2double(Res2); %Convert string to double
 
%--Phase selection ---
idx_P2=find(contains(Ph2,'Pg'));   %Test itPh(idx_P)
idx_S2=find(contains(Ph2,'Sg')); 


%---------------------------------------    
% File 3 - Relocation with P+S models (simultaneous inversion) using EDT
  file='Relocation_P+SModel_EQ20200617_02.55.25.csv';
 [Ph3 NET3 Res3 Dist3] = csvimport( file, 'columns', {'PHASE', 'NET', 'TimeRes(s)', 'Dist(km)'});
 % Res3=str2double(Res3); %Convert string to double
 
%--Phase selection ---
idx_P3=find(contains(Ph3,'Pg'));   %Test itPh(idx_P)
idx_S3=find(contains(Ph3,'Sg')); 



%==Itersection of conditions== 
%idx=intersect( idx_region,idx_mag, idx_z);
%idx=intersect(intersect(idx_region,idx_mag,'stable'),idx_z,'stable');  %3 cases
 
 sz=10; %Markersize
 lw=2; % %Linewidht marker
 xlt=[0,max(Dist1)+2];
 ylt=[-0.6 0.6]; %Yaxis
 ytk=[-0.4 0 .4];%[-.8:.4:max(ylim) ];

 
 figure(1)
subplot(3,1,1)
 plot(Dist1(idx_P1),Res1(idx_P1),'o','Markersize',sz,'LineWidth',lw)
 hold on
  plot(Dist1(idx_S1),Res1(idx_S1),'sq','Markersize',sz,'LineWidth',lw)
  plot(0*ones(1,length(Res1)),'k')
 set(gca,'Fontsize',25)
legend('P-phases','S-phases','fontsize',20, 'Location','NorthEast');
%xlabel('Distance (km)')
ylabel('Residual (s)')
title('HUP2')
ylim([ylt])
xlim([xlt])
grid on
set(gca,'YTick',ytk)

subplot(3,1,2)
 plot(Dist2(idx_P2),Res2(idx_P2),'o','Markersize',sz,'LineWidth',lw)
 hold on
  plot(Dist2(idx_S2),Res2(idx_S2),'sq','Markersize',sz,'LineWidth',lw)
  plot(0*ones(1,length(Res2)),'k')
 set(gca,'Fontsize',25)
legend('P-phases','S-phases','fontsize',20, 'Location','NorthEast');
%xlabel('Distance (km)')
ylabel('Residual (s)')
title('P & S')
ylim([ylt])
xlim([xlt])
grid on
set(gca,'YTick',ytk)

subplot(3,1,3)
 plot(Dist3(idx_P3),Res3(idx_P3),'o','Markersize',sz,'LineWidth',lw)
 hold on
  plot(Dist3(idx_S3),Res3(idx_S3),'sq','Markersize',sz,'LineWidth',lw)
  plot(0*ones(1,length(Res3)),'k')
 
 set(gca,'Fontsize',25)
legend('P-phases','S-phases','fontsize',20, 'Location','NorthEast');
xlabel('Distance (km)')
ylabel('Residual (s)')
title('P + S')
ylim([ylt])
xlim([xlt])
grid on
set(gca,'YTick',ytk)

 fig.PaperUnits = 'inches';
 fig.PaperPosition =[0 5 14 17]  %[0 0 x y];
% print(gcf,['Hist_Depth_MLabove1_H1_2'],'-dpng','-r500');  %dpi set at 800 to increase image resolution 
 



 %%  PLot Min 1D (OUTPUT VELEST inversion) when using each of these Vmodels for relocation
 
 %1) HU-P2 ISOR model as used in SeiscomP3 for routinery relocation
 %   /home/sysop/nll/data/NLL.iceland_HU-P2_1d.conf

 vp1=[ 3.22 4.13 4.96 5.77 6.32 6.41 6.5 6.53 6.6 6.66];
 z1=[-.2 .8 1.8 2.8 3.6 3.9 5 5.8 7.8 9.8];
 
 vs1=round(vp1/1.78,2);   %Simply they used vp/vs=1.78 cte
 
 
 %2) P and S model from individual VELEST inversions
 %CS_VMGralTrendMiddle_invertratio1_P_CNVFeb2020_2'; %PREFERED VP Vmodel
%CS_VMGralTrendMiddle_invertratio1_S_ratio_rm'; %PREFERED VS model

 vp2=[2.7200    3.2300    3.7800    4.3000    4.8100    5.6600    6.3000    6.6100    6.8300 ... 
    6.8300    6.8300    7.0200    7.0200    7.0200    7.0200    7.0200    7.0800    7.1700 7.2600];

vs2=[1.6 1.67 1.91 2.16 2.71 3.25 3.59 3.82 3.88 3.88 3.91 3.95 3.95 3.95 3.95 3.95 3.97 4.02 4.07]; 
 
z2=  [-1.0000         0    0.5500    1.1000    1.6000    2.2000    2.9000    4.2000    5.3300...
6.4700    7.6000    9.4700   11.3300   13.2000   15.0600   17.0000   19.0000   21.0000 25.0000];
 
%3) Model from simultaneous P+S  VELEST inversion
 %/P+S_Inversion/CS_invertratio1_PS

 vp3=[2.69 3.17 3.72 4.26 4.85 5.77 6.39 6.79 7.0 7.0 7.01 7.40 7.40 7.40 7.40 7.40 7.40 7.40 7.40 ];
 
 vs3=[1.63 1.69 1.89 2.15 2.77 3.35 3.61 3.8 3.9 3.98 4.06 4.06 4.06 4.07  4.07  4.07  4.07 4.07  4.07 ];
 
 z3=z2;
 
 %--Vp/Vs Ratio--
 VRatio=vp3./vs3;
 
 
  L_width=2;
% %==== Plot all Models Together Vp and Vs===
% figure(2)
%  ax4=subplot(2,2,[1,2,3]);
% %--Vp--
% Try=stairs(vp1,z1,'k--','LineWidth',L_width); %HUP1  (V,Z)
% hold on
% stairs(vp2,z2,'k:','LineWidth',L_width);   %Min1D Vp individual
% stairs(vp3,z3,'k','LineWidth',L_width); %Min1D Vp from P+S inversion
% 
% %--Vs--
% stairs(vs1,z1,'k--','LineWidth',L_width); %HUP1  (V,Z)
% PS=stairs(vs2,z2,'k:','LineWidth',L_width); %Min1D Vs individual
% PplusS=stairs(vs3,z3,'k','LineWidth',L_width);  %Min1D Vs from P+S inversion
% 
% set(ax4,'yaxislocation','left')
% 
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.3, 0.9]);% [0, 0.04, 1, 0.96]);
% set(gca,'Fontsize',25)
% legend([Try(1),PS(1),PplusS(1)],'HUP2','P&S','P+S','fontsize',20,'Location','Southwest');
% 
% view([0 -90 ])
% xlabel('Velocities (km/s)')
% ylabel('Depth (km)')
% ylim([0 15])
% xlim([1.5 8])
% grid on
% set(gca,'YTick',[0:2:max(ylim) ])
% set(gca,'XTick',[floor(min(xlim)):1:max(xlim) ])
% box on
% 
% % 

%%  Final Models (from P+S) compared with Literature 
%--Max Depth for Plotting
mxdepth=10; %max depth to plot (km)

%===Wagner et al. 2019===
vp_wg=[3.0603
3.2766
3.3904
3.9937
4.4148
5.2686
5.6442
6.1337
6.1934
6.2845
6.3129
6.4837
6.5463
6.6316
6.6829
6.6883
6.7736
6.8192
6.8761];

% z_wg=[-1 
% -0.5322
% -0.0179
% 0.4842
% 0.974
% 1.4884
% 1.9782
% 2.4925
% 2.9915
% 3.4936
% 3.9895
% 4.4915
% 4.9875
% 5.4956
% 5.9916
% 6.4867
% 6.9887
% 7.4969
% 8.0051]; %8.0051 last depth  I skipped
%  

z_wg=[-1:0.5:8];



%Vs 
d=load('/Users/alejandro/ICELAND_Project/Min1D_Model/Pilar_files/alex/Vs_Wagner.mat')  %Extracted from the (red line) Fig8 in  Wagner et al 2019

D=d.Vs_Wagner(1:2:end,:);  %Actual points for the VELEST Vmodel  .mod file

vs_wg=D(:,1); %km/s 

%===Tryggvason et al 2002=====
vp_try=[3.40 3.60 4.70 5.60 6.10 6.40 6.50 6.60 6.70 6.80 6.90 7.00 7.10 7.20 7.40 7.40];

z_try=[-2 0 1 2 3 4 5 6 9 14 15 16 23 24 25 32  ];

vs_try=[1.95 2.05 2.65 3.15 3.45 3.60 3.65 3.75 3.80 3.85 3.90 3.95 4 4.05 4.20 4.20]; %km/s -> From Table2 in Trygvasson et al 2002


  
%==== Plot all Models Together Vp and Vs and vp/vs===
figure(3)
 ax4=subplot(2,2,[1,2,3]);
%--Vp--
Try=stairs(vp_try,z_try,'k--','LineWidth',L_width); %Tryggvason  (V,Z)
hold on
stairs(vp_wg,z_wg,'k:','LineWidth',L_width);   %Wagner
stairs(vp3,z3,'k','LineWidth',L_width); %Min1D Vp from P+S inversion

%--Vs--
stairs(vs_try,z_try,'k--','LineWidth',L_width); 
Wg=stairs(vs_wg,z_wg(1:end-1),'k:','LineWidth',L_width); 
PplusS=stairs(vs3,z3,'k','LineWidth',L_width);  %Min1D Vs from P+S inversion


set(ax4,'yaxislocation','left')

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.3, 0.9]);% [0, 0.04, 1, 0.96]);
set(gca,'Fontsize',25)
%legend([Try(1),Wg(1),PplusS(1)],'Tryggvason et al. 2002','Wagner et al. 2019','Final Model P+S','fontsize',20,'Location','Southwest');
legend([Try(1),Wg(1),PplusS(1)],'Tryggvason et al. 2002','Wagner et al. 2019','Final Model P+S','fontsize',20,'Location','Southwest');

view([0 -90 ])
xlabel('Velocities (km/s)')
ylabel('Depth (km)')
ylim([0 mxdepth]) %15
xlim([1.5 8])
grid on
set(gca,'YTick',[0:2:max(ylim) ])
set(gca,'XTick',[floor(min(xlim)):1:max(xlim) ])
box on



 figure(4)
%--Vp/Vs Ratio--
rat=stairs(VRatio,z3,'k','LineWidth',L_width+1);  

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.12, 0.9]);% [0, 0, 0.3, 0.9]); %[0, 0.04, 1, 0.96]);
set(gca,'Fontsize',25)
legend([rat(1)],'V_{P}/V_{S}','fontsize',20, 'Location','Southwest');
view([0 -90 ])
xlabel('Velocity Ratio')
ylabel('Depth (km)')
ylim([0 mxdepth])
xlim([1.6 2.1])
grid on
set(gca,'YTick',[0:2:max(ylim) ])
%set(gca,'XTick',[floor(min(xlim)):0.1:max(xlim) ])
%set(gca,'XTick',[0:2:max(ylim) ])
%     axes('YAxisLocation','right');

box on
 

%%   Plot P+S Models with Focal Depth Histogram in Tobias Style

%-- Relocated Manual events (Vp and Vs relocation are the same because it is a P+S inversion)--
EQFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/P+S_Inversion/CS_invertratio1_PS/velest2mat_files/velest2mat.latlondepmag';

fid=fopen(EQFile);
EQ=textscan(fid,'%f %f %f %f'); % "%*f" --> skip a value
fclose(fid);
DEP=EQ{3}; %Vp
DEP_s=DEP; %


%-- Plot parameters--
Fsz=30%20; %Font Size
histBW=0.2; %Focal depth histogram bin width in %km



%-------
%close all
figure (5)
%--Vs--
ax1=subplot(1,3,1);
Try=stairs(vs_try,z_try,'k--','LineWidth',L_width); 
hold on
Wg=stairs(vs_wg,z_wg(1:end-1),'k:','LineWidth',L_width); 
PplusS=stairs(vs3,z3,'k','LineWidth',L_width);  %Min1D Vs from P+S inversion

set(ax1,'xaxislocation','top') %left
%set(ax1,'yaxislocation','left')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.9]);% [0, 0.04, 1, 0.96]);
set(gca,'Fontsize',Fsz)
legend([Try(1),Wg(1),PplusS(1)],'Tryggvason et al. 2002','Wagner et al. 2019','Final Model P+S','fontsize',20,'Location','Southwest');
view([0 -90 ])
xlabel('Vs (km/s)')
ylabel('Depth (km)')
ylim([0 mxdepth])
xlim([1 5])
grid on
set(gca,'YTick',[0:2:max(ylim) ])
set(gca,'XTick',[floor(min(xlim)):2:max(xlim) ])
box on


%--Vp--
ax1=subplot(1,3,2);
Try=stairs(vp_try,z_try,'k--','LineWidth',L_width); %Tryggvason  (V,Z)
hold on
Wg=stairs(vp_wg,z_wg,'k:','LineWidth',L_width);   %Wagner
PplusS=stairs(vp3,z3,'k','LineWidth',L_width); %Min1D Vp from P+S inversion

set(ax1,'xaxislocation','top') %left
%set(ax1,'yaxislocation','left')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.9]);% [0, 0.04, 1, 0.96]);
set(gca,'Fontsize',Fsz)
%legend([Try(1),Wg(1),PplusS(1)],'Tryggvason et al. 2002','Wagner et al. 2019','Final Model P+S','fontsize',20,'Location','Southwest');
view([0 -90 ])
xlabel('Vp (km/s)')
%ylabel('Depth (km)')
ylim([0 mxdepth])
xlim([3 7.5])
grid on
set(gca,'YTick',[0:2:max(ylim) ])
set(gca,'XTick',[floor(min(xlim)):2:max(xlim) ])
box on



%--Focal Depth Vp---
ax2=subplot(1,3,3);
%hp=histogram(DEP)
hp=histogram(DEP,'FaceColor',rgb('Silver')) %Using rgb Function -Type "rgb chart" to see all colors 
hp.BinWidth= histBW;  
view([90 90 ])
ylabel('No. Events')
xlim([0 mxdepth])
ylim([0 13])
%set(gca,'xtick',[])


%ax = gca;
set(ax2,'yaxislocation','right')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.9]);% [0, 0.04, 1, 0.96]);
set(gca,'Fontsize',Fsz)



%==== Plot all Models Together Vp and Vs===
figure(5)
 ax4=subplot(2,2,[1,2,3]);
%--Vp--
Try=stairs(vp1,z1,'k--','LineWidth',L_width); %HUP1  (V,Z)
hold on
stairs(vp2,z2,'k:','LineWidth',L_width);   %Min1D Vp individual
Wg=stairs(vp_wg,z_wg,'k-.','LineWidth',L_width);   %Wagner
stairs(vp3,z3,'k','LineWidth',4); %Min1D Vp from P+S inversion


%--Vs--
stairs(vs1,z1,'k--','LineWidth',L_width); %HUP1  (V,Z)
Wg=stairs(vs_wg,z_wg(1:end-1),'k-.','LineWidth',L_width); 
PS=stairs(vs2,z2,'k:','LineWidth',L_width); %Min1D Vs individual
PplusS=stairs(vs3,z3,'k','LineWidth',4);  %Min1D Vs from P+S inversion


set(ax4,'yaxislocation','left')

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.3, 0.9]);% [0, 0.04, 1, 0.96]);
set(gca,'Fontsize',30)
%legend([Try(1),Wg(1),PS(1),PplusS(1)],'Tryggvason','Wagner','P&S','P+S','fontsize',20,'Location','Southwest');
%legend([Try(1),Wg(1),PS(1),PplusS(1)],'Tryggvason et al. 2002','Wagner et al. 2019','P&S','Final Model P+S','fontsize',20,'Location','Southwest'); %P_{ind},S_{ind},=P&S
lg=legend([Try(1),Wg(1),PS(1),PplusS(1)],'Tryggvason et al. 2002','Wagner et al. 2019','P_{ind},S_{ind}','P+S','fontsize',20,'Location','Southwest'); %P_{ind},S_{ind},=P&S
lg.FontSize=30

view([0 -90 ])
xlabel('Velocities (km/s)')
ylabel('Depth (km)')
ylim([0 mxdepth]) %15
xlim([1.5 8])
grid on
set(gca,'YTick',[0:2:max(ylim) ])
set(gca,'XTick',[floor(min(xlim)):1:max(xlim) ])
box on

% HUP2

