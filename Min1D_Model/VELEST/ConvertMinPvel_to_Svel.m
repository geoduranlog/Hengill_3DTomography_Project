%================= Convert output of Min 1D P Vmodel to S ===================
clear all 
close all


%--- SELECT INPUT FILES---
%Select label of Vmodels (tests) you want to compare
%----------------------------------
name='CS_VMGralTrendMiddle_invertratio1_P_CNVFeb2020_2';  %This was your best Solution for P-phases
folder=['/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/',name,'/']

%--Station Residual FILES  - Output From VELEST---
FILE=[folder,name,'.mod']; 


% Load 
fileid1 = fopen(FILE,'rt');
C=textscan(fileid1,' %f %f  %f ','headerLines',2); %  |vel|z(km)|damping|
fclose(fileid1);

%--Vp/Vs Ratio--
v_ratio=1.78; %  Wagner 1.77, Tryggvason (2002) 1.78 -> Using joint Wadati-Inversion

%Convert vp to vs
C{1}=C{1}/v_ratio;

%{
%---Print S Vmodel---
fileID = fopen('S_Vmodel.mod','w');
fprintf(fileID,'S-Vmodel built from Final output vel. P-model using vp/vs cte Ratio., |Vel.| depth(top)|damp.|\n');
fprintf(fileID,'%2.0f \n', length(C{1}) ); %Nb Layers
%fprintf(fileID,'%2.4f %2.4f %2.2f\n',[C{1}  C{2} C{3}]);
fprintf(fileID,'%6.3f   %6.2f %6.2f\n',[C{1}';  C{2}'; C{3}']); %Keep the spaces
fclose(fileID);
%}


%% High Low Test  S-waves.



%--- SELECT INPUT FILES---
%Select label of Vmodels (tests) you want to compare
%----------------------------------
name1='CS_VMGralTrendMiddle_invertratio1_S_ratio_rm'; %Prefix of the VELEST output files
name2='CS_VMGralTrendMiddle_invertratio1_S_ratio_lowTest';   %Low test
name3='CS_VMGralTrendMiddle_invertratio1_S_ratio_highTest';   %High test

folder1=['/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/S_Waves/',name1,'/']
folder2=['/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/S_Waves/',name2,'/']
folder3=['/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/S_Waves/',name3,'/']



%-- Select Models models
MOD_IN1=[folder1,'S_Vmodel.mod'];
MOD_OUT1=[folder1,name1,'.mod'];

MOD_IN2=[folder2,'S_LowTest.mod'];
MOD_OUT2=[folder2,name2,'.mod'];

MOD_IN3=[folder3,'S_HighTest.mod'];
MOD_OUT3=[folder3,name3,'.mod'];

%----- Load Vmodels -- INPUTS
fileid1 = fopen(MOD_IN1,'rt');
ModIn1=textscan(fileid1,'%f %f %*f','headerLines',2); %
fclose(fileid1);

fileid2 = fopen(MOD_IN2);
ModIn2=textscan(fileid2,'%f %f %*f','headerLines',2); %
fclose(fileid2);

fileid3 = fopen(MOD_IN3);
ModIn3=textscan(fileid3,'%f %f %*f','headerLines',2); %
fclose(fileid3);


%----- Load Vmodels -- OUTPUTS
fileid1 = fopen(MOD_OUT1);
ModOut1=textscan(fileid1,'%f %f %*f','headerLines',2); %
fclose(fileid2)

fileid2 = fopen(MOD_OUT2);
ModOut2=textscan(fileid2,'%f %f %*f','headerLines',2); %
fclose(fileid2)

fileid3 = fopen(MOD_OUT3);
ModOut3=textscan(fileid3,'%f %f %*f','headerLines',2); %
fclose(fileid3);




%% Create Low and High Vmodels
%v_low = [0.5618    0.8427    1.2921    1.7416    2.1348    2.5281    2.9213    3.1180    3.2584    3.3708 3.5112    3.7079    3.8390     ModIn1{1}(14:end)'];
v_low = [0.8    1    1.4    1.7416    2.1348    2.5281    2.9213    3.34    3.52    3.59    3.65    3.72    3.8390     ModIn1{1}(14:end)'];

v_high =  [  2.5356    2.8864    3.0687    3.2509    3.4894    3.5843    3.75    3.8    3.86    3.88    3.9         3.944       3.944   ModIn1{1}(14:end)'];


%-- Plot Input S-Vmodel and its low and high Vmodels
figure(100)
 ax4=subplot(2,2,[1,2,3]);
  stairs(ModIn1{1}, ModIn1{2},'-o')
 hold on
  stairs(v_low, ModIn1{2},'-o')
 hold on
  stairs(v_high, ModIn1{2},'-o')
set(ax4,'yaxislocation','left')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.3, 0.9]); %[0, 0.04, 1, 0.96]);
set(gca,'Fontsize',25)
legend({'In. S Model' ,'Model low','Model high'},'Location','southWest','FontSize',16)
view([0 -90 ])
xlabel('Velocities (km/s)')
ylabel('Depth (km)')
ylim([0 25]) %25
xlim([0 4.5]) %25
grid on


%-- Plot Results after inversion
L_width=2;
figure(9)
ax4=subplot(2,2,[1,2,3]);
%Low
stairs(ModIn2{1},ModIn2{2},'r--','LineWidth',L_width);
hold on
stairs(ModOut2{1},ModOut2{2},'r-','LineWidth',L_width);
hold on
%High
 stairs(ModIn3{1},ModIn3{2},'b--','LineWidth',L_width);
 hold on
stairs(ModOut3{1},ModOut3{2},'b-','LineWidth',L_width);
 hold on
 %Reference
stairs(ModOut1{1},ModOut1{2},'g-','LineWidth',L_width); %New Model Out

set(ax4,'yaxislocation','left')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.3, 0.9]); %[0, 0.04, 1, 0.96]);
set(gca,'Fontsize',25)
legend({'In. Low','Out. Low','In. High','Out. High','Out. S Min1D'},'Location','southWest','FontSize',16)
view([0 -90 ])
xlabel('Velocities (km/s)')
ylabel('Depth (km)')
ylim([0 15])%25
set(gca,'YTick',[0:2:15 ])
grid on