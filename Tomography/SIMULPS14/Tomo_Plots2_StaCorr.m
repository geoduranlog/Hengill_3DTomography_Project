%==== Plot TOMOGRAPHY good style + StaCorrections====
%I used SIMULPS for the inversion. Then, Tobias codes to get the Hoz slides
%with all the info.  THey are saved per each test in the folder called VP 
clear all
close all

%--P phases---
name='2ndRound/Coseismiq_3x3x2_d8_it3_pmtsJun12'; % Name of inversion Test  Coseismiq_3x3x1p5_d100_it4  2ndRound/Coseismiq_5x5x3_d10_it2
Folder=['/Users/alejandro/ICELAND_Project/Tomography/SIMULPS14/TEST_simulps/',name,'/VP/'];

% %--S phases---
%name='CoseismiqSw_3x3x2_d8_it3'; % Name of inversion Test
%Folder=['/Users/alejandro/ICELAND_Project/Tomography/SIMULPS14/TEST_simulps/2ndRound/S_Waves/',name,'/VP/'];

%3x3x1.5
% VFile1=[Folder,'plane0000.0.xyz'];  %Velocity Profile to plot
% VFile2=[Folder,'plane0001.5.xyz'];  %Velocity Profile to plot
% VFile3=[Folder,'plane0003.0.xyz'];  %Velocity Profile to plot
% VFile4=[Folder,'plane0004.5.xyz'];  %Velocity Profile to plot

%3x3x2
 VFile1=[Folder,'plane0000.0.xyz'];  %Velocity Profile to plot
 VFile2=[Folder,'plane0002.0.xyz'];  %Velocity Profile to plot
 VFile3=[Folder,'plane0004.0.xyz'];  %Velocity Profile to plot
 VFile4=[Folder,'plane0006.0.xyz'];  %Velocity Profile to plot

%5x5x3
%VFile1=[Folder,'plane0000.0.xyz'];  %Velocity Profile to plot
%VFile2=[Folder,'plane0003.0.xyz'];  %Velocity Profile to plot
%VFile3=[Folder,'plane0006.0.xyz'];  %Velocity Profile to plot
%VFile4=[Folder,'plane0009.0.xyz'];  %Velocity Profile to plot


%VFile=[Folder,planes(1)];  %Velocity Profile to plot


%--Synthetic test----
% name='modplot_ADNtemp_coseismiq'; % Name of inversion Test
% Folder=['/Users/alejandro/ICELAND_Project/Tomography/SIMULPS14/Tobi_Codes/SIMULPS14/Testcase_CH_SyntheticTest/',name,'/VP/'];
% 
% VFile1=[Folder,'plane_in00-1.0.xyz'];  %Velocity Profile to plot
% VFile2=[Folder,'plane_in0000.0.xyz'];  %Velocity Profile to plot
% VFile3=[Folder,'plane_in0001.5.xyz'];  %Velocity Profile to plot
% VFile4=[Folder,'plane_in0003.0.xyz'];  %Velocity Profile to plot


fid=fopen(VFile1);  % [x,y,lon,lat,z,vp,pvp(%),dvp(absol.),phit,pkhit,pdws,presol.,spread_p,vpvs,pvpvs(%),skhit,sdws,sresol,spread_s ]
Vel1=textscan(fid,'%*f %*f   %f %f %f %f %f   %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f','headerLines',2); % %"%*f" --> skip a value
fclose(fid);
LON1=Vel1{1};    
LAT1=Vel1{2}; 
DEP1=Vel1{3};
v_ref1=Vel1{4}(1);
dv1=Vel1{5};

fid=fopen(VFile2);  % [x,y,lon,lat,z,vp,pvp(%),dvp(absol.),phit,pkhit,pdws,presol.,spread_p,vpvs,pvpvs(%),skhit,sdws,sresol,spread_s ]
Vel2=textscan(fid,'%*f %*f   %f %f %f %f %f   %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f','headerLines',2); % %"%*f" --> skip a value
fclose(fid);
LON2=Vel2{1};    
LAT2=Vel2{2}; 
DEP2=Vel2{3};
v_ref2=Vel2{4}(1);
dv2=Vel2{5};

fid=fopen(VFile3);  % [x,y,lon,lat,z,vp,pvp(%),dvp(absol.),phit,pkhit,pdws,presol.,spread_p,vpvs,pvpvs(%),skhit,sdws,sresol,spread_s ]
Vel3=textscan(fid,'%*f %*f   %f %f %f %f %f   %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f','headerLines',2); % %"%*f" --> skip a value
fclose(fid);
LON3=Vel3{1};    
LAT3=Vel3{2}; 
DEP3=Vel3{3};
v_ref3=Vel3{4}(1);
dv3=Vel3{5};

fid=fopen(VFile4);  % [x,y,lon,lat,z,vp,pvp(%),dvp(absol.),phit,pkhit,pdws,presol.,spread_p,vpvs,pvpvs(%),skhit,sdws,sresol,spread_s ]
Vel4=textscan(fid,'%*f %*f   %f %f %f %f %f   %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f','headerLines',2); % %"%*f" --> skip a value
fclose(fid);
LON4=Vel4{1};    
LAT4=Vel4{2}; 
DEP4=Vel4{3};
v_ref4=Vel4{4}(1);
dv4=Vel4{5};

%--Load Stations--
rcv_file=[Folder,'sta.xyz'];  
fid=fopen(rcv_file); 
rcv=textscan(fid,'%f %f  %*f %*f %*f   %*s'); % %"%*f" --> skip a value
fclose(fid);
LON_rcv=rcv{1};    
LAT_rcv=rcv{2}; 

%}

%% %--  Shape Files from ISOR ---
 dir0='/Users/alejandro/ICELAND_Project/Tomography/Pilar_files_tomo_CheckHDISK_MyPassport/ShapeFiles_Iceland';
  dir_shape=[dir0,'/shapefiles'];

%lat=[63.7 64.25]; lon=[-22 -20.7];
% lat=[63.98 64.14]; lon=[-21.48 -21.17]; % Zoom central area


roads_q='yes'; faults_q='yes'; faults_prob_q='yes';
faults_uncer_q='no'; eruptive_q='no';

 
%-- Load shape Files---
if strcmp(roads_q,'yes')
    disp('Adding roads...');
    roads=shaperead([dir_shape,'/roads/roads_latlon/Hengill_vegir.shp'],'UseGeoCoords',true);
    mapshow(roads,'Color','b')
end

if strcmp(faults_q,'yes')
    disp('Adding faults...');
    faults=shaperead([dir_shape,'/faults/faults_latlon/hengill_brot.shp'],'UseGeoCoords',true);
    mapshow(faults,'Color','k')
end

if strcmp(faults_prob_q,'yes')
    disp('Adding probable faults...');
    faults_prob=shaperead([dir_shape,'/faults_probable/faults_probable_latlon/LiI-kleg-brot.shp'],'UseGeoCoords',true);
    mapshow(faults_prob,'Color','k')
end

if strcmp(faults_uncer_q,'yes')
    disp('Adding uncertain faults...');
    faults_uncer=shaperead([dir_shape,'/faults_uncertain/faults_uncertain_latlon/oviss_brot.shp'],'UseGeoCoords',true);
    mapshow(faults_uncer,'Color','k')
end


%--Geothermal power plants---
%Hellisheiði
%lat_Helli=64.0373; lon_Helli=-21.4010;  %From google - but this is just an exhibition center
lat_Helli=64.0372; lon_Helli=-21.4008;  
%Nesjavellir
%lat_Nesj=64.1082; lon_Nesj=-21.2571;  %google
lat_Nesj=64.1082; lon_Nesj=-21.2560;    %Pili
lat_geoplant=[lat_Helli lat_Nesj]; lon_geoplant=[lon_Helli lon_Nesj];

% Zoom area
latz=[63.96 64.16]; lonz=[-21.55 -21.1]; 

%latz=[63.9 64.15]; lonz=[-21.5 -21.1])  %other, bigger map

%--Volcanos--
%Coordinates from http://icelandicvolcanos.is/#
%Hengill Volcano 
lat_HeVolc=64.083;
lon_HeVolc=-21.330; %-21.416;

%Hrómundartindur Volcano
lat_HrVolc=64.0667; %64.0769;
lon_HrVolc=-21.202; %-21.2167; %-21.1997;

%Grensdalur Volcano
lat_GrVolc=64.0464; %64.0458;
lon_GrVolc=-21.1905; %-21.1918;

 

%%  DEFINE COLOR MAP
%From red (low) to blue (high)

%close all
%L = 10;             %number of datapoints
%data = 3.6*rand(L); % create example data set with values ranging from 0 to 3.6

% L=50;
% a = min(dv); b = max(dv);
% data = (b-a).*rand(L,L) + a;  %r_range = [min(data) max(data)]


L=length(dv1);
data=dv1;

indexValue = 0;     % value for which to set a particular color
topColor = [0 0 1];         % color for maximum data value (blue = [0 0 1])
indexColor = rgb('Gainsboro')  % color for indexed data value (white = [1 1 1])
bottomcolor = [1 0 0];      % color for minimum data value (red = [1 0 0])

% Calculate where proportionally indexValue lies between minimum and
% maximum abs values        max( [max(dv1),max(dv2),max(dv3),max(dv4)] )
                      %     min( [min(dv1),min(dv2),min(dv3),min(dv4)] )
abs_max=max( [max(abs(dv1)),max(abs(dv2)),max(abs(dv3)),max(abs(dv4))] ); %ABS MAX
abs_max=15;

largest =abs_max; %abs( min(min(data)) );    %max(max(data)); 
smallest =-abs_max; % min(min(data));
%largest =max(max(data)); 
%smallest = min(min(data));
index = L*abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
num=100;
customCMap1 = [linspace(bottomcolor(1),indexColor(1),num*index)',...
            linspace(bottomcolor(2),indexColor(2),num*index)',...
            linspace(bottomcolor(3),indexColor(3),num*index)'];
% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topColor(1),num*(L-index))',...
            linspace(indexColor(2),topColor(2),num*(L-index))',...
            linspace(indexColor(3),topColor(3),num*(L-index))'];
customCMap = [customCMap1;customCMap2];  % Combine colormaps
%customCMap = [customCMap2;customCMap2];
%colormap(customCMap)
%colormap(flipud(customCMap));
%psudo = pcolor(data);
%colorbar


%[val, idx] = max(dv);  LON(idx)  LAT(idx)


%% STATION CORRECTIONS
% RCV WITH CORRECTIONS
A=200;
%P=waves
rcvCorrFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/CS_VMGralTrendMiddle_invertratio1_P_CNVFeb2020_2/velest2mat_files/velest2mat.statcorr';

%S=waves
%rcvCorrFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/S_Waves/CS_VMGralTrendMiddle_invertratio1_S_ratio/velest2mat_files/velest2mat.statcorr';

%P+S inversion
%rcvCorrFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/P+S_Inversion/CS_invertratio1_PS_iusestacorr1/velest2mat_files/velest2mat.statcorr';


fid1=fopen(rcvCorrFile);
STAT=textscan(fid1,'%s %f %f %*d %*d %d %f %f %d %d'); % "%*f" --> %* skip a value
fclose(fid1);

NAME=STAT{1};
LAT_rcv=STAT{2}; LON_rcv=STAT{3};
RCV_IDNb=STAT{4};
PTCOR=STAT{5}; %P-wave time station correction
STCOR=STAT{6}; %S-wave time station correction.  unused for now
POBS=STAT{7};  %Nb observations. Works for both P or S inversion

%---P+S inversion--
%POBS=STAT{8}; % ONLY if you used P+S inversion and you want to consider Nb observations in S
%PTCOR=STCOR;  % ONLY if you used P+S inversion and you want to plot S station corr

 %Reference Station index  (labeled as 999)
ref_rcvidx=find(STAT{4}==max(STAT{4})); 

idx_maj=find(PTCOR>0);
idx_min=find(PTCOR<0);
idx_zero=find(PTCOR==0.00);

% Condition Plot only rcv with POBS>=idx_obs   %ADN 2020
%idx_obs=find(POBS>=10);
idx_obs=10; 
idx_obs=find(POBS>=idx_obs); 
idx_maj=intersect(idx_maj,idx_obs);
idx_min=intersect(idx_min,idx_obs);
idx_zero=intersect(idx_zero,idx_obs);

%--Size --
size_symb=2e5; %4.  %Size of the symbols used to plot rcv correction.  If 'geom' mode  %ADN


%%  FIGURES
%--Symbol sizes within map--
A=100; %Size symbol stations
A_gplant=320;%Size symbol GethermalPlants

scalnb=1e6; %Scale weigths of colormap (Densityplot) so it matches with dv values
posi=[0.5703    0.1100    0.2844    0.3412]; 

r_smooth=1500; % Original=1.2549e+03 %Radius of interpolation densityplot

%======FIGURE TOMOGRAPHY=======
close all
figure(3)
%----Plane 1----
%subplot(2,2,1)
colormap(customCMap);
%colormap(flipud(redblue)); %Flip the order of the colormap
%dp=geodensityplot(LAT1,LON1,dv1*scalnb,'FaceColor','interp')  %*1e6    %Radius of influence on density calculation 1.2549e+03m
dp=geodensityplot(LAT1,LON1,dv1*1e6,'FaceColor','interp','Radius', r_smooth)  %*1e6    %Radius of influence on density calculation 1.2549e+03m

hold on
%-----Faults----
for i=1:584
geoplot(faults(i).Lat,faults(i).Lon,'k')
end
%-----Roads----
for i=1:567
geoplot(roads(i).Lat,roads(i).Lon,'color',rgb('Gray'),'Linewidth',2.5)
end
%--- Stations ---
s=geoscatter(LAT_rcv,LON_rcv,A,'^','filled')
s.MarkerFaceColor = rgb('Orange');
s.MarkerEdgeColor = rgb('DarkRed');  %Using RGB Function
%--- Geothermal Plants----
s=geoscatter(lat_geoplant,lon_geoplant,A_gplant,'h','filled')
s.MarkerFaceColor = rgb('Lime');
s.MarkerEdgeColor = rgb('Black');

% %---Volcanos----
% Hengill
s=geoscatter(lat_HeVolc,lon_HeVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black')
%Hrómundartindur
s=geoscatter(lat_HrVolc,lon_HrVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black')
%Grensdalur Volcano
s=geoscatter(lat_GrVolc,lon_GrVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black');

%--Station Corrections ---  
%Negative Correction
             hs_min=geoscatter(LAT_rcv(idx_min),LON_rcv(idx_min), ...
                abs(PTCOR(idx_min))*size_symb,'ob','Linewidth',1);
%Positive Corretion
  hs_max=geoscatter(LAT_rcv(idx_maj),LON_rcv(idx_maj), ...
                 PTCOR(idx_maj)*size_symb,'+r','Linewidth',1);    
             
                 %Plot Ref Station
           hs_rcv_ref=geoplot(LAT_rcv(ref_rcvidx),LON_rcv(ref_rcvidx),'v','MarkerSize',18)
           hs_rcv_ref.MarkerFaceColor = rgb('Green'); 
           hs_rcv_ref.MarkerEdgeColor = rgb('Black'); 

colb=colorbar;
set(gca,'fontsize',25)
%set(gx,'TickLabelForma','dd') % Use decimal degrees in axis
set(dp.Parent,'TickLabelForma','dd') % Use decimal degrees in axis
geobasemap none  %Erase Basemap
%colormapeditor  
%gx.Scalebar.Visible = 'on';
%geolimits([63.9 64.15],[-21.5 -21.1]);
geolimits(latz,lonz);
caxis([-abs_max  abs_max])
title('z=0.0km')
%dp.Parent.LongitudeAxis.Visible='off'
grid off
set(get(colb,'label'),'string','V_{P} (%)','Rotation',90);

%{
%----Plane 2----
  subplot(2,2,2)
colormap(customCMap);
dp=geodensityplot(LAT2,LON2,dv2*1e6,'FaceColor','interp','Radius', r_smooth)  %*1e6    %Radius of influence on density calculation 1.2549e+03m
hold on
%-----Faults----
for i=1:584
geoplot(faults(i).Lat,faults(i).Lon,'k')
end
%-----Roads----
for i=1:567
geoplot(roads(i).Lat,roads(i).Lon,'color',rgb('Gray'),'Linewidth',2.5)
end
%--- Stations ---
s=geoscatter(LAT_rcv,LON_rcv,A,'^','filled')
s.MarkerFaceColor = rgb('Orange');
s.MarkerEdgeColor = rgb('DarkRed');  %Using RGB Function
%--- Geothermal Plants----
s=geoscatter(lat_geoplant,lon_geoplant,A_gplant,'h','filled')
s.MarkerFaceColor = rgb('Lime');
s.MarkerEdgeColor = rgb('Black');
% %---Volcanos----
% Hengill
s=geoscatter(lat_HeVolc,lon_HeVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black')
%Hrómundartindur
s=geoscatter(lat_HrVolc,lon_HrVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black')
%Grensdalur Volcano
s=geoscatter(lat_GrVolc,lon_GrVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black');

colb=colorbar;
set(gca,'fontsize',25)
set(dp.Parent,'TickLabelForma','dd') % Use decimal degrees in axis
geobasemap none  %Erase Basemap
geolimits(latz,lonz);
caxis([-abs_max  abs_max])
title('z=2.0km')
dp.Parent.LatitudeAxis.Visible='off'
dp.Parent.LongitudeAxis.Visible='off'
ylabel(colb,'dv_{P} (%)')
set(get(colb,'label'),'string','V_{S} (%)','Rotation',90);
grid off

%----Plane 3----
 subplot(2,2,3)
colb=colormap(customCMap);
%dp=geodensityplot(LAT3,LON3,dv3*scalnb,'FaceColor','interp')  %*1e6    %Radius of influence on density calculation 1.2549e+03m
dp=geodensityplot(LAT3,LON3,dv3*1e6,'FaceColor','interp','Radius', r_smooth)  %*1e6    %Radius of influence on density calculation 1.2549e+03m
hold on
%-----Faults----
for i=1:584
geoplot(faults(i).Lat,faults(i).Lon,'k')
end
%-----Roads----
for i=1:567
geoplot(roads(i).Lat,roads(i).Lon,'color',rgb('Gray'),'Linewidth',2.5)
end
%--- Stations ---
s=geoscatter(LAT_rcv,LON_rcv,A,'^','filled')
s.MarkerFaceColor = rgb('Orange');
s.MarkerEdgeColor = rgb('DarkRed');  %Using RGB Function
%--- Geothermal Plants----
s=geoscatter(lat_geoplant,lon_geoplant,A_gplant,'h','filled')
s.MarkerFaceColor = rgb('Lime');
s.MarkerEdgeColor = rgb('Black');
% %---Volcanos----
% Hengill
s=geoscatter(lat_HeVolc,lon_HeVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black')
%Hrómundartindur
s=geoscatter(lat_HrVolc,lon_HrVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black')
%Grensdalur Volcano
s=geoscatter(lat_GrVolc,lon_GrVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black');

colorbar;
set(gca,'fontsize',25)
set(dp.Parent,'TickLabelForma','dd') % Use decimal degrees in axis
geobasemap none  %Erase Basemap
%geolimits([63.9 64.15],[-21.5 -21.1]);
geolimits(latz,lonz);
caxis([-abs_max  abs_max])
title('z=4.0km')
grid off

%----Plane 4----
 subplot(2,2,4)
colb=colormap(customCMap);
%dp=geodensityplot(LAT3,LON3,dv3*scalnb,'FaceColor','interp')  %*1e6    %Radius of influence on density calculation 1.2549e+03m
dp=geodensityplot(LAT4,LON4,dv4*1e6,'FaceColor','interp','Radius', r_smooth)  %*1e6    %Radius of influence on density calculation 1.2549e+03m

hold on
%-----Faults----
for i=1:584
geoplot(faults(i).Lat,faults(i).Lon,'k')
end
%-----Roads----
for i=1:567
geoplot(roads(i).Lat,roads(i).Lon,'color',rgb('Gray'),'Linewidth',2.5)
end
%--- Stations ---
s=geoscatter(LAT_rcv,LON_rcv,A,'^','filled')
s.MarkerFaceColor = rgb('Orange');
s.MarkerEdgeColor = rgb('DarkRed');  %Using RGB Function
%--- Geothermal Plants----
s=geoscatter(lat_geoplant,lon_geoplant,A_gplant,'h','filled')
s.MarkerFaceColor = rgb('Lime');
s.MarkerEdgeColor = rgb('Black');
% %---Volcanos----
% Hengill
s=geoscatter(lat_HeVolc,lon_HeVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black')
%Hrómundartindur
s=geoscatter(lat_HrVolc,lon_HrVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black')
%Grensdalur Volcano
s=geoscatter(lat_GrVolc,lon_GrVolc,A_gplant,'s','filled')
s.MarkerFaceColor = rgb('Fuchsia');
s.MarkerEdgeColor = rgb('Black');

colb=colorbar;
set(gca,'fontsize',25)
set(dp.Parent,'TickLabelForma','dd') % Use decimal degrees in axis
geobasemap none  %Erase Basemap
%geolimits([63.9 64.15],[-21.5 -21.1]);
geolimits(latz,lonz);
caxis([-abs_max  abs_max])
title('z=6.0km')
dp.Parent.LatitudeAxis.Visible='off'
set(get(colb,'label'),'string','V_{S} (%)','Rotation',90);
grid off
%}

%Size of the Window to show and print later
    % set(gcf, 'Position',  [100, 100, 1416.04, 1081.92])
    set(gcf, 'Position',  [100, 100, 1416.04/2, 1081.92/2])

     %--Check actual limits plotted by Matlab--
%[latlim, lonlim] = geolimits

%--Save---
        %- To be able to save as real pdf
        fig1=figure(3);
fig1.Renderer='Painters';
save_folder='/Users/alejandro/ICELAND_Project/Tomography/Figures/Tomography_Vp/';
 %print(gcf,[save_folder,'Tomo_HengilVs_prmtJun12'],'-dpdf','-r800');  %dpi set at 800 to increase image resolution 
