%-----MAP OF EVENTS USED FOR CONSTRUCTION OF MIN 1D VMODEL --------
% USE .CNV Output file from VELEST, in velest2mat


%%
clear all
close all
clc
   
 %P-waves  
%EQFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/CS_VMGralTrendMiddle_invertratio1_P_CNVFeb2020_2/velest2mat_files/velest2mat.latlondepmag';

%S-waves
EQFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/S_waves/CS_VMGralTrendMiddle_invertratio1_S_ratio_rm/velest2mat_files/velest2mat.latlondepmag';


fid=fopen(EQFile);
EQ=textscan(fid,'%f %f %f %f'); % "%*f" --> skip a value
fclose(fid);
LAT=EQ{1};   
LON=EQ{2}; LON=-LON; % Longuitude 'West' will be read as Negative in Matlab =>  Make it a negative number
DEP=EQ{3};
MAG=EQ{4};
    

% Convert into a Table to use geobubble function
EQ=table(LAT,LON,DEP,MAG);


Severity = discretize(MAG, [ min(MAG) 1.5  2  max(MAG)],...
'categorical', {'[0, 1.4]','[1.5, 1.9]', '[2.0, 3.0]'});
EQ.ML=Severity;

%---Select BaseMap
url = 'a.tile.opentopomap.org';
name = 'opentopomap';
copyright = char(uint8(169));
%attribution = ["OpenTopoMap"]
attribution = [ ...
      copyright + "OpenTopoMap"];
     %"map style: " + copyright + "OpenTopoMap (CC-BY-SA)"];
  displayName = 'Open Topo Map';
  addCustomBasemap(name,url,'Attribution',attribution,'DisplayName',displayName)

 
 
%% FIGURES
%{
    figure(1)
 gb = geobubble(EQ,'LAT','LON', ...    %ISORtest            %ColorVariable      'colorterrain' 'bluegreen'
        'SizeVariable','MAG', ...
        'Basemap','opentopomap','ColorVariable','ML') %'color_scale')       %'colordata'  from 1 to 7. ->you can create a vector which is normalized based on ML values
gb.BubbleWidthRange =[5 60] %Range of sizes of the circles
title('Events used Min 1D Vmodel'); %,'fontsize',14);
set(gca,'fontsize',30) %23
%geolimits([63.9 64.2],[-21.6 -21.1])
 
 %}
 

% %--Test Stations Location
RCVFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/CS_VMGralTrendMiddle_invertratio1_P_CNVFeb2020_2/velest2mat_files/velest2mat.statcorr';

fid=fopen(RCVFile);
rcv=textscan(fid,'%s %f %f %f %f %*f %*f  %*f  %*f  %*f'); % "%*f" --> skip a value
fclose(fid);
NAME=rcv{1}; 
LAT_rcv=rcv{2};   
LON_rcv=rcv{3}; 


%Separete RCV and EQ TEST 
LAT1=LAT(1:50); LAT2=LAT(51:end);
LON1=LON(1:50); LON2=LON(51:end);

LAT_rcv1=LAT_rcv(1:40); LAT_rcv2=LAT_rcv(41:end);
LON_rcv1=LON_rcv(1:40); LON_rcv2=LON_rcv(41:end);



% figure(2)
% geobasemap('opentopomap')
% %lat=rcv.LAT; lon=rcv.LON; %-> Most be changed to the REAL locations
% geoscatter(LAT,LON,'^')



 % ===BASEMAPS===
 %{
 %--usgstopo--
url = "https://basemap.nationalmap.gov/ArcGIS/rest/services";
fullurl = url + "/USGSTopo/MapServer/tile/${z}/${y}/${x}";
nm = 'usgstopo';
att = 'Credit: US Geological Survey';


%--- World_Topo_Map from Esri---
url = "http://services.arcgisonline.com/ArcGIS/rest/services";
fullurl = url + "/World_Topo_Map/MapServer/tile/${z}/${y}/${x}";
nm = 'World_Topo_Map';
att = 'Credit: Esri';

addCustomBasemap(nm,fullurl,'Attribution',att)
  %}


 
 %--  Shape Files from ISOR ---
 dir0='/Users/alejandro/ICELAND_Project/Tomography/Pilar_files_tomo/ShapeFiles_Iceland';
  dir_shape=[dir0,'/shapefiles'];

lat=[63.7 64.25]; lon=[-22 -20.7];
% lat=[63.98 64.14]; lon=[-21.48 -21.17]; % Zoom central area


roads_q='yes'; faults_q='yes'; faults_prob_q='yes';
faults_uncer_q='no'; eruptive_q='yes';%'yes';


% rgb=imread(['hengill_large_co.png']);
% [lx, ly, lz]=size(rgb);
% inix=-22; finix=-20-42/60; dxa=(abs(inix-finix))/lx;
% iniy=64+16/60; finiy=63+40/60; dya=(abs(iniy-finiy))/ly;
% lons=[inix:dxa:finix-dxa]; lats=[iniy:-dya:finiy-dya];
% image(lons,lats,rgb); view([0 -90]); hold on
% xlim(lon); ylim(lat);
% xlabel('Longitude (\circE)')
% ylabel('Latitude (\circN)')
% set(gca,'FontSize',14)

 
 

 
  

A=320; %Size symbol stations
Mag_sz=200*(MAG+1);   %Size EQ symbols

%close all
figure(2)
geobasemap('World_Topo_Map')%('usgstopo') %('opentopomap')
hold on
%--- Stations ---
geoscatter(LAT_rcv1,LON_rcv1,A,'c^','filled') 
geoscatter(LAT_rcv2,LON_rcv2,A,'k^','filled') 
%--- EQ ---
geoscatter(LAT1,LON1,Mag_sz(1:50),'bo') 
geoscatter(LAT2,LON2,Mag_sz(51:end),'ro')
%geoplot(faults.Lat,faults.Lon,'k')
geolimits([63.9 64.15],[-21.4 -21.1])
%Read from the map the dms coordinates u like and convert it to decimal decimal_degrees=deg + mins/60 + secs/3600;
%lon_dm=[-(21 +35/60)  -(21 +05/60)  ] ;   lat_dm=[(63 +56/60)  (64 +08/60)  ] ;     
%geolimits([lat_dm], [lon_dm])
legend('a','b','c','d')
set(gca,'fontsize',30)


if strcmp(faults_q,'yes')
    disp('Adding faults...');
    faults=shaperead([dir_shape,'/faults/faults_latlon/hengill_brot.shp'],'UseGeoCoords',true);
    %mapshow(faults,'Color','k')
    %geoshow(faults,'Color','k')
end





%--Print--  HighRes, takes some time
%{
 fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 13 10]; %[0 0 16 13];  %[0 0 15 12];
print(gcf,'EQ_Min1DModel2','-dpdf','-r800');  %dpi set at 800 to increase image resolution 
 %}



%Find Largest Seismic Events (ML)
%max(Data_C.Mag)
[fil]=find(MAG>2.7);
EQ(fil,:);  %22560     29724     29725


