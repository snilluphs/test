function [lambda, specrad, bbcm, bbph] = shp_BB_Integrated_Irradiance(Tc, lowerWL, upperWL, omega, to, ta)
% INPUTS
% Tc = black body temperature (deg. C)
% lowerWL = lower wavelength (microns)
% upperWL = upper wavelength (microns)
% omega = solid angle subtended (steradians)
% to = optical system transmission spectrum (from lowerWL to upperWL)
% ta = atmospheric transmission spectrum (from lowerWL to upperWL)
% OUTPUTS
% bbcm = integrated irradiance in units of W/cm^2
% bbph = integrated irradiance in units of photons/cm^2/sec
% NOTES
% requires function shp_BBcurve.m
% uses trapezoidal rule for integration
% copyright 2021, Steve Pullins


%% required constants
h = 6.62607015e-34;  % J*s (Planck's constant)
csp = 299792458;     %m/s (speed of light)
%% black body spectral radiance calc
[lambda, specrad] = shp_BBcurve(Tc, lowerWL, upperWL);
%% Integrate the BB curve
% Integrated irradiance in W/cm2
bbcm = 1e-4*omega*trapz(lambda, to.*ta.*specrad); % integrated blackbody (W/cm2)
% Integrated irradiance in photons/cm2/sec
bbph = 1e-4*omega*trapz(lambda, (to.*ta.*specrad.*1e-6.*lambda/(h*csp)));
