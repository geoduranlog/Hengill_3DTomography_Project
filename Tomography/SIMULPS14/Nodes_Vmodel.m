%% ===========Model Parameterization for SimulPS =================
% Set a node distance of the 3D initial Vmodel

clear all 
close all


%% Case 1x1x0.5   %dist always in km
%Here I assume that node region could be < Large Coseismiq region as defined in
%Simulps.  I have to test if it runs

x=0:1:7;
z=0:0.5:5.5;   


nz_pls=length(z)-1;  %Hyp, Should be <19 (19 was Nb layers in Min1DVmodel)
nx_pls=length(x)-1;

x=[x 14 ]
z=[-0.5 z 7 15] %7 10 15]

% Nb nodes
nx=2*length(x)-1
nz=length(z)   

%Total Nb Nodes + Hpar < NbPhases
Hpar=130*4; %Hypocenter parameters ->NbEQ*4
nx*nx*nz+Hpar  %Should be <4900  .  Nb P phases is ~5100

%% Case 1x1x0.5  
%Here I assume that node region could be < Large Coseismiq region as defined in
%Simulps.  I have to test if it runs

x=0:1:6;
z=0:0.5:5.5;   


nz_pls=length(z)-1;  %Hyp, Should be <19 (19 was Nb layers in Min1DVmodel)
nx_pls=length(x)-1;

x=[x 14 50 ]
z=[-1 z 7 15] %7 10 15]

% Nb nodes
nx=2*length(x)-1
nz=length(z)   

%Total Nb Nodes + Hpar < NbPhases
Hpar=130*4; %Hypocenter parameters ->NbEQ*4
nx*nx*nz+Hpar  %Should be <4900  .  Nb P phases is ~5100
%% Case 2x2x1   %%%%%%%%NEW
%Large COSEISMIQ Region as defined in Simulps LAT [60 , 70] LON [-25 , -15]
%To cover my entire area with the nodes, starting form the center, I’d need to reach
%xmax~277.5km and ymax~555.5km. However, to incluse all stations I only
%need xmax~40km and ymax~25km. 

%Limits of my area of interest (km)
x=0:2:12;
y=0:2:14;
z=0:1:6 ;

nz_pls=length(z)-1  %Hyp, Should be <19 (19 was Nb layers in Min1DVmodel)
nx_pls=length(x)-1

%Extend (Last point is boundary, Penultimo point is to include all stations)
x=[x 40 200]
y=[y 25 120 ]
z=[-10 -1 z 8 10 50]   %Inlude -10 as your Boundary

% Nb nodes
nx=2*length(x)-1
ny=2*length(y)-1
nz=length(z)

%Total Nb Nodes + Hpar < NbPhases
Hpar=130*4; %Hypocenter parameters ->NbEQ*4
nx*ny*nz+Hpar  %Should be <4900  .  Nb P phases is ~5100. But only 3773P were used in Min1DVmodel


%% Case 3x3x1.5
%Large COSEISMIQ Region as defined in Simulps LAT [60 , 70] LON [-25 , -15]
%To cover my entire area with the nodes, starting form the center, I’d need to reach
%xmax~277.5km and ymax~555.5km.  However, to incluse all stations I only
%need xmax~40km and ymax~25km. 


%Limits of my area of interest (km)
x=0:3:12;
y=0:3:15;
z=0:1.5:6 ;

nz_pls=length(z)-1  %Hyp, Should be <19 (19 was Nb layers in Min1DVmodel)
nx_pls=length(x)-1

%Extend (Last point is boundary, Penultimo point is to include all stations)
x=[x 15 40 200]
y=[y 25 120]
z=[-10 -1 z 8 10 15 50]

% Nb nodes
nx=2*length(x)-1
ny=2*length(y)-1
nz=length(z)

%Total Nb Nodes + Hpar < NbPhases
Hpar=130*4; %Hypocenter parameters ->NbEQ*4
nx*ny*nz+Hpar  %Should be <4900  .  Nb P phases is ~5100


%% Case 5x5x3
%Large COSEISMIQ Region as defined in Simulps LAT [60 , 70] LON [-25 , -15]
%To cover my entire area with the nodes, starting form the center, I’d need to reach
%xmax~277.5km and ymax~555.5km.  However, to incluse all stations I only
%need xmax~40km and ymax~25km. 


%Limits of my area of interest (km)
x=0:5:15;
y=0:5:15;
z=0:3:6 ;

nz_pls=length(z)-1  %Hyp, Should be <19 (19 was Nb layers in Min1DVmodel)
nx_pls=length(x)-1

%Extend (Last point is boundary, Penultimo point is to include all stations)
x=[x 40 200]
y=[y 25 120]
z=[-10 -1 z 9 12 50]

% Nb nodes
nx=2*length(x)-1
ny=2*length(y)-1
nz=length(z)

%Total Nb Nodes + Hpar < NbPhases
Hpar=130*4; %Hypocenter parameters ->NbEQ*4
nx*ny*nz+Hpar  %Should be <4900  .  Nb P phases is ~5100   Nb S phases is ...





%% Input to convert Min1D into a 3DVmodel
% To Copy on make3DsimulMod.cmn  file

x=[-x x(2:end)];
x=sort(x,'ascend')


y=[-y y(2:end)];
y=sort(y,'ascend')

format bank  %Show2 digits after the decimal point.
z


format short %Go back to standard