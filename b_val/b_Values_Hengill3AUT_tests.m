%====== Compute b-value Hengill   ======
%== from Francesco's reprocessed data (using EDT norm) ====
%01.October.2020
clear all
close all
clc

% ------ Tobias Automatic Data (Reliable(HQ) or Full Data) Period 01.12.2018-20.08.2020  -------
name='/Users/alejandro/ICELAND_Project/DATA_FILES/QC_AUT_Data/Data_LucaandTobias'

%HQ Selection
DataFile=[name,'/tobias-data/reliableCat/SC3DB_ExtractEvents_AutoCat_PB_EDT.2C.Events.out.Reliable_S1.0.csv']; 

%Full catalogue
%DataFile=[name,'/tobias-data/fullCat/SC3DB_ExtractEvents_AutoCat_PB_EDT.2C.Events.csv']; 

[lon lat z mag yr mth day hr min sec gap EQscore Lmeth ] = csvimport( DataFile, 'columns', {'Lon', 'Lat', 'Depth','Mag','YYYY','MM','DD','HH','MI','SS','GAP','EVScore','LMethod',});

%----Manual Data SED---

%Convert to dateTime format
time=[yr , mth, day, hr, min, sec];
time=datetime(time);


%-- Select Region ---
%Hengill Region
lat_Hengill=[63.9 64.15];  lon_Hengill=[-21.6 -21.1 ]; %According to Tobias selected Region

%Index of the region 
  idx_regAUT=find(lat>=lat_Hengill(1) & lat<=lat_Hengill(2) & lon>=lon_Hengill(1) & lon<=lon_Hengill(2) ); 

%Take values only within selected region 
time=time(idx_regAUT); 
lat=lat(idx_regAUT);  lon=lon(idx_regAUT); 
z=z(idx_regAUT); 
mag=mag(idx_regAUT); 




figure(1)
hist2d(z,mag,50) %50
 xlabel('Depth (km)')
 ylabel('Magnitude (ML)')
 zlabel('Counts')
  set(gca,'fontsize',25)
  grid on
  xlim([0 10])
  ylim([-1 4.5])

 
%-Hist Nb EQ Coseismiq-
figure (2)
histogram(time)
ylabel('No. Events')
xlabel('Date')
box on
set(gca,'fontsize',25)
  
%%  B-Value per Clusters

%Define CLusters

%Hengill Region
lat_Hengill=[63.9 64.15];  lon_Hengill=[-21.6 -21.1 ]; %According to Tobias selected Region

lat_h1=[64.03  64.07 ]; lon_h1=[-21.42 -21.35]; %Husmuli
lat_h2=[64.04  64.07 ]; lon_h2=[-21.29 -21.15]; %Skardsnyrarfjall
lat_h3=[63.98  64.02 ]; lon_h3=[-21.39 -21.31]; %Hverahlid
lat_h4=[63.99 64.03 ];  lon_h4=[-21.45 -21.39 ]; %Grauhnjukar - I used a larger area compared to the one in the siesmicity map
lat_h5=[64.09 64.13 ];  lon_h5=[-21.3  -21.22 ]; %Nesjavellir
lat_h6=[64.09 64.13 ];  lon_h6=[-21.4  -21.32 ]; %Left from Nesj.
lat_h7=[63.93 63.98 ];  lon_h7=[-21.46 -21.27 ]; %CentralSouth, contains 4.5ML EQ
lat_h8=[63.94 64 ];  lon_h8=[-21.21 -21.14 ]; %EastSouth, contains also a large, ML~3 EQ


%Index per clusters
 idx_h1=find(lat>=lat_h1(1) & lat<=lat_h1(2) & lon>=lon_h1(1) & lon<=lon_h1(2) ); 
 idx_h2=find(lat>=lat_h2(1) & lat<=lat_h2(2) & lon>=lon_h2(1) & lon<=lon_h2(2) ); 
 idx_h3=find(lat>=lat_h3(1) & lat<=lat_h3(2) & lon>=lon_h3(1) & lon<=lon_h3(2) ); 
 idx_h4=find(lat>=lat_h4(1) & lat<=lat_h4(2) & lon>=lon_h4(1) & lon<=lon_h4(2) ); 
 idx_h5=find(lat>=lat_h5(1) & lat<=lat_h5(2) & lon>=lon_h5(1) & lon<=lon_h5(2) ); 
 idx_h6=find(lat>=lat_h6(1) & lat<=lat_h6(2) & lon>=lon_h6(1) & lon<=lon_h6(2) );
 
 idx_h7=find(lat>=lat_h7(1) & lat<=lat_h7(2) & lon>=lon_h7(1) & lon<=lon_h7(2) );
 idx_h8=find(lat>=lat_h8(1) & lat<=lat_h8(2) & lon>=lon_h8(1) & lon<=lon_h8(2) );

        %Entire Hengill
 idx_h9=find(lat>=lat_Hengill(1) & lat<=lat_Hengill(2) & lon>=lon_Hengill(1) & lon<=lon_Hengill(2) ); 

 
 %All names in one string array
 idx_h=  ["idx_h1" , "idx_h2", "idx_h3", "idx_h4", "idx_h5", "idx_h6", "idx_h7", "idx_h8" ];
 cluster_names=["Húsmúli" , "Skarsmýrarfjall", "Hveragerði", "Grauhnjukar", "Nesjavellir", "H6", "H7", "H8" ];
 

%%  B-value using Max Likelihood Estimation (Tobias' code)

cd 'b_value_TobiasCode'

%--Save ML values of each cluster in an ASCII file
file_names=  ["H1_Mag" , "H2_Mag", "H3_Mag", "H4_Mag", "H5_Mag", "H6_Mag", "H7_Mag", "H8_Mag"];
%mag( eval(idx_h(i)) );
for i=1: length(file_names)
writematrix(mag( eval(idx_h(i)) ),file_names(i)) 
end

%writematrix(mag( idx_h9 ),"H9_Mag") ; %Hengill


%--Parameters
mc=1.9;  %All clusters had Mc~0, I got it from various tests (see bellow)
mbin=0.01
catalogID='H cluster';
magtype='ML';

%--- B-value Fitting and FIGURE
%OJO check fBinning value on calc_bmemagMag Code!
close all
[fBValue, fStdDev, fAValue] = bvalue_maxlike("H7_Mag.txt",mc,mbin,catalogID,magtype);
xlim([-1  4.5])
%ylim([10^-0.5 10^4])
set(gca,'fontsize',25)

fBValue
fStdDev

%{
%--Get Mc (Magnitude of Completeness by using various cutt off values and
%study fluctuations on b--- Do it for each Cluster
clear Mco fBValue b_McTest std_McTest min max
Mco=[min(mag(idx_h9)):0.1:max(mag(idx_h9))]  ;% cut-off to tests for Mc on a given cluster
for i=1:length(Mco)
    
    [fBValue, fStdDev, fAValue] = bvalue_maxlike("H9_Mag.txt",Mco(i),mbin,catalogID,magtype);

    b_McTest(i)=fBValue;
    std_McTest(i)=fStdDev;
end


close all
b_valH9=  0.9361;  %for deltaM=0.01
figure(5)
errorbar(Mco(1:length(b_McTest)),b_McTest,std_McTest,'ko','Markersize',8);
hold on
plot(Mco(1:length(b_McTest)),ones(1,length(b_McTest)).*b_valH9,'--k')
xlabel('Mco')
ylabel('b-value' )
box on
legend('b-value estimates','correct b-value','Location','Northwest')
set(gca,'fontsize',25)
xlim([-1.5 4])  %2.7
ylim([0 3]) %6.2
%ylim([0 1.5])



%}

cd ../