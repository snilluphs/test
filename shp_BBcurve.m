function [lambda, specrad] = shp_BBcurve(Tc, lowerWL, upperWL)
% INPUTS
% Tc = temperature (C)
% lowerWL = starting wavelength (microns)
% upperWL = ending wavelength (microns)
% OUTPUTS
% lambda  = array of wavelengths (microns)
% specrad = array of blackbody spectral rad. (W m-2 sr-1 um-1)
% NOTES
% copyright 2021, Steve Pullins

%% Temperature
Tsk = Tc + 273.15; %temperature in Kelvin
%% Wavelength
wlstep_um = 0.01; %10nm
dlambda = upperWL-lowerWL;
lambda = (linspace(lowerWL, upperWL, (dlambda/wlstep_um)))'; %wavelength array (microns)
%% Constants
h = 6.62607015e-34;  % J*s (Planck's constant)
csp = 299792458;     %m/s (speed of light)
kB = 1.38064853e-23; %J/K (Boltzmann's constant)
C1 = 2*csp*csp*h; %2*csp*csp*h; %W*m2
C2 = h*csp/kB;       %m*K
%% Blackbody  
lambdam = lambda*1e-6; %convert wavelength from microns to meter
specrad  = (C1./(lambdam.^5.0)).*(1./(exp(C2./(lambdam*Tsk)) -1.0)); %W sr-1 m-3
specrad = specrad*(1e-6); %convert FROM W sr-1 m-3 TO: W sr-1 m-2 um-1
