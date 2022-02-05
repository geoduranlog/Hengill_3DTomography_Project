tomo%==== Max Resolution for Tomography====
clear all
close all


%% MOD 3x3x2
%Diagonal of such cube (l)
x=3; y=3; z=2; %km
l=sqrt( x.^2+y.^2 + z.^2 );

%Reduce the number - Most of rays won't have the diagonal lenght (which is the max ray length within the cube)
l=4;  


% Parameters
vo=6.4; %km/s  -> Background velocity (from Min1D) it varies with depth, select a meaningful value
vo=vo/1.78  %For S-waves Hengill
ter=0.0215; %<Picking Error P> (s)
ter_s=0.0622; %<Picking Error S> (s) for Phases used in the tomo inversion (S-min1D data)

%--You can resolve dv (or pvp) only when |dt|>|ter| --

%dt=ter;  %
%dt=0.02

%Percent-velocity-perturbation  pvp=100*dv/vo  (%)
%pvp=-100*vo*dt/(l+vo*dt)  
pvp=6

%Travel time residual between background model and perturbed model
dt=(l/vo)*(-1/(100/pvp+1)) 


%% Vp/Vs Tomo Resolution
vp=6.4; %km/s
ko=1.78;  %Vp/Vs Ratio Hengill

ter_ps=ter+ter_s


%Percent-velocity-perturbation of vp/vs  pkp  (%)
pkp=8;

dt_ps=(l*ko*pkp)/(vp*100)
