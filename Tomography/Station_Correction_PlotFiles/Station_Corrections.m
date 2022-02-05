% Satellital map with faults and roads
close all
clear all
clc
%% ------------------ HEADER ------------------%
dir0='/Users/alejandro/ICELAND_Project/Tomography/Pilar_files_tomo/ShapeFiles_Iceland';
dir_shape=[dir0,'/shapefiles'];

sav='no';

roads_q='no'; faults_q='no'; faults_prob_q='no';
faults_uncer_q='no'; eruptive_q='no';%'yes';

% ------------------ END HEADER ------------------%

%% ------ STATION INFO -------
%  RCVFile='/Users/alejandro/ICELAND_Project/All_Events_Sept2018_Feb2020/Manual_events_Sept18_Feb20/data_20182019/station.csv';
%  [LAT_rcv LON_rcv z net rcv_name] = csvimport( RCVFile, 'columns', {'latitude', 'longitude', 'elevation','networkCode','stationCode'});
% % 
% % A=200;


% RCV WITH CORRECTIONS
A=200;
%---Hengill Regional Model
%P=waves
%rcvCorrFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/CS_VMGralTrendMiddle_invertratio1_P_CNVFeb2020_2/velest2mat_files/velest2mat.statcorr';

%S=waves
%rcvCorrFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/S_Waves/CS_VMGralTrendMiddle_invertratio1_S_ratio/velest2mat_files/velest2mat.statcorr';

%P+S inversion
%rcvCorrFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/P+S_Inversion/CS_invertratio1_PS_iusestacorr1/velest2mat_files/velest2mat.statcorr';


%---Nesjavellir Local Model
%P-waves
%rcvCorrFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/Data_Summer2020/TESTS/Nesjavellir_Iceland/CS_NES_VM_PS_invertratio1_P_Origin_NESRCV_Nearby/velest2mat_files/velest2mat.statcorr';

%S-waves
rcvCorrFile='/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/Data_Summer2020/TESTS/Nesjavellir_Iceland/S-waves/CS_NES_VM_PS_invertratio1_S_Origin_NESRCV_Nearby/velest2mat_files/velest2mat.statcorr';


fid1=fopen(rcvCorrFile);
STAT=textscan(fid1,'%s %f %f %*d %*d %d %f %f %d %d'); % "%*f" --> %* skip a value
fclose(fid1);

NAME=STAT{1};
LAT_rcv=STAT{2}; LON_rcv=STAT{3};
RCV_IDNb=STAT{4};
PTCOR=STAT{5}; %P-wave time station correction
STCOR=STAT{6}; %S-wave time station correction.  unused for now
POBS=STAT{7};  %Nb observations. Works for both P or S inversion

%---P+S inversion--S Info only
%POBS=STAT{8}; % ONLY if you used P+S inversion and you want to consider Nb observations in S
%PTCOR=STCOR;  % ONLY if you used P+S inversion and you want to plot S station corr

 %Reference Station index  (labeled as 999)
ref_rcvidx=find(STAT{4}==max(STAT{4})); 

idx_maj=find(PTCOR>0);
idx_min=find(PTCOR<0);
idx_zero=find(PTCOR==0.00);

% Condition Plot only rcv with POBS>=idx_obs   %ADN 2020
idx_obs=10; 
idx_obs=find(POBS>=idx_obs); 
idx_maj=intersect(idx_maj,idx_obs);
idx_min=intersect(idx_min,idx_obs);
idx_zero=intersect(idx_zero,idx_obs);

%--Size --
size_symb=2e4; %4.5e4; %2.4^14;  %2.1^10  %Size of the symbols used to plot rcv correction.  If 'geom' mode  %ADN


%% Gelogical image -snapshot
%--For small map--
%geoimage=imread(['GeoMap_Hengil.png']);
geoimage=imread(['GeoMap_Hengil_opaco.png']);
[lx, ly, lz]=size(geoimage);
inix=-21.77; finix=-20.60; dxa=(abs(inix-finix))/lx;
iniy=64.18; finiy=63.85; dya=(abs(iniy-finiy))/ly;    %y starts at top left! 


%--For larger map-- 
% geoimage=imread(['GeoMap_Hengil_larger.png']);
% [lx, ly, lz]=size(geoimage);
% inix=-22.31; finix=-20.80; dxa=(abs(inix-finix))/lx;
% iniy=64.16; finiy=63.84; dya=(abs(iniy-finiy))/ly;


%--For BIG map--
%geoimage=imread(['GeoMap_Hengil_Big.png']);
%[geoimage, map, alpha] = imread('GeoMap_Hengil_Big.png');

% f=imshow(geoimage);
% set(f, 'AlphaData', 0.6) %Save as png

%geoimage=imread(['GeoMap_Hengil_Big.png']);
%geoimage=imread(['GeoMap_Hengil_Big_opaco2.png']);
% [lx, ly, lz]=size(geoimage);
% inix=-22.54; finix=-20.03; dxa=(abs(inix-finix))/lx;
% iniy=64.35; finiy=63.77; dya=(abs(iniy-finiy))/ly;



%Get LON LAT from Image
 lons=[inix:dxa:finix-dxa]; lats=[iniy:-dya:finiy-dya];


%Area to plot
%lat_plot=[63.77 64.3]; lon_plot=[-22.2 -20.75];
%lat_plot=[63.88 64.14]; lon_plot=[-21.68 -20.9];

lat_plot=[63.9055   64.1390];  lon_plot=[-21.64  -21]; %Actual limits used in Map_ISOR_SEDManual_faults  


%Legend Size
%Find corrections close to the <P_picking error>=0.04 ->This will be your min
%error to see
min_sybl=min(find(PTCOR==-0.04)); % I just take the 1st case. Only one is needed as a reference
max_sybl=min(find(PTCOR==0.04)); % PTCOR(min_sybl) PTCOR(max_sybl)

%--Just to visualize correction size that I want
figure(1)
ax = gca();
image(ax,lons,lats,geoimage);  view([0 -90]);
hold on
hs_min=scatter(ax,LON_rcv(min_sybl),LAT_rcv(min_sybl), ...
                abs(PTCOR(min_sybl))*size_symb,'ob','Linewidth',1);
          
            hs_max=scatter(ax,LON_rcv(max_sybl),LAT_rcv(max_sybl), ...
                 PTCOR(max_sybl)*size_symb,'+r','Linewidth',1);
             [h,icons,plots,legend_text] =legend([hs_min(1),hs_max(1),],'-0.04','+0.04','fontsize',20);
for k =3:4 %Increase Symbol size of legend
icons(k).Children.MarkerSize = 30;
end
xlim(lon_plot); ylim(lat_plot);
xlabel('Longitude (\circW)')
ylabel('Latitude (\circN)')
set(gca,'FontSize',22)
set(gca,'position',[  0.1124    0.1100    0.6698    0.8150]);



%close all
figure(2)
ax = gca();
image(ax,lons,lats,geoimage);  view([0 -90]);
hold on
%--- Stations ---
% s=scatter(ax, LON_rcv,LAT_rcv,100,'^','filled');
%  s.MarkerFaceColor = rgb('Black'); %DimGray;
%  s.MarkerEdgeColor = rgb('Snow');  %Using RGB Function

%Negative Correction
    hs_min=scatter(ax,LON_rcv(idx_min),LAT_rcv(idx_min), ...
                abs(PTCOR(idx_min))*size_symb,'ob','Linewidth',1);
                %abs(PTCOR(idx_min))*(2^10),'vb'); % abs for minor
                               
                %Plor Rcv Symbol
            hs_min_rcvsym=plot(LON_rcv(idx_min),LAT_rcv(idx_min),'v','MarkerSize',8, ...
                'MarkerFaceColor',[0.9 0.9 0.9],'MarkerEdgeColor','k'); 
            
 %Positive Corretion
  hs_max=scatter(ax,LON_rcv(idx_maj),LAT_rcv(idx_maj), ...
                 PTCOR(idx_maj)*size_symb,'+r','Linewidth',1);
                 %PTCOR(idx_maj)*(2^10),'^r');
                 
           %Plor Rcv Symbol
            hs_max_rcvsym=plot(ax,LON_rcv(idx_maj),LAT_rcv(idx_maj),'v','MarkerSize',8, ...
           'MarkerFaceColor',[0.9 0.9 0.9],'MarkerEdgeColor','k'); %Plot rcv symbol
       
       %Plot Ref Station
           hs_rcv_ref=plot(ax,LON_rcv(ref_rcvidx),LAT_rcv(ref_rcvidx),'v','MarkerSize',18)
           hs_rcv_ref.MarkerFaceColor = rgb('Green'); 
           hs_rcv_ref.MarkerEdgeColor = rgb('Black'); 
 %hold(ax, 'on');
% hold(ax, 'off')
[h,icons,plots,legend_text] =legend([hs_min(1),hs_max(1),hs_rcv_ref(1)],'-0.04s','+0.04s','Ref. Station','fontsize',20,'Location','Northwest');
for k =4:5 %Increase Symbol size of legend
icons(k).Children.MarkerSize = 25;  % This represent the size for +-0.04s
end
%xlim(lon_plot); ylim(lat_plot);
set(gca,'YTick',[63.95:0.05:64.1 ])
xlabel('Longitude (\circW)')
ylabel('Latitude (\circN)')
set(gca,'FontSize',22)
set(gca,'position',[  0.1124    0.1100    0.6698    0.8150]);

get(gca,'position');




%{
% Shape files
if strcmp(eruptive_q,'yes')
    disp('Adding roads...');
    eruptive=shaperead([dir_shape,'/eruptive_info/eruptive_info_latlon/eldv_polygon.shp'],'UseGeoCoords',true);
    mapshow(eruptive,'FaceColor','r')
end

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



% Saving
if strcmp(sav,'yes')
    disp('Saving...');
    set(gcf,'PaperPositionMode','auto')
    %%print(gcf,[dir0,'/maps/satel_shapes.png'],'-dpng','-r500')
   % print(gcf,[dir0,'/satel_shapes.png'],'-dpng','-r500')
end
 %}


mean(PTCOR(idx_min)) %Avg of negative correction 
mean(PTCOR(idx_maj)) %Avg of positive correction 
