%% Set parameters of the system
function [Parameters]=setParameters()
WaveLengthMin=200e-9; %lower bound on wavelength
WaveLengthMax=900e-9; %upper bound on wavelength
NumberOfWaveLength=100; %number of wavelenght
WaveLengthToPlot=500e-9; %give wavelength for which the intensity distribution will be plotted
d=600e-9; %period of the grating
h=400e-9; %depth of the grating
Polarization={'TM'}; %set polarization
%Asymetrical profile
a_fun=@(x) (0.1*h)*cos((2*pi/d).*x) + (0.02*h)*cos((4*pi/d).*x - 5*pi/9);
a_diff_fun=@(x) -(0.1*h*2*pi/d).*sin((2*pi/d).*x) - (0.02*h*4*pi/d).*(sin((4*pi/d).*x - 5*pi/9));
%a_fun=@(x) 0.*x;
%a_diff_fun=@(x) 0.*x;
%Symmetrical profile
%a_fun=@(x) (h/2)*cos((2*pi/d).*x);
%a_diff_fun=@(x) -((h/2)*2*pi/d).*sin((2*pi/d).*x);
Profile={a_fun, a_diff_fun};
ThetaInc=15; %set incident angle
n1=1; %refraction index in incident media
%n2=1.5; %refraction index in transmission media
n2=1+5i;
numTr=20; %Truncation number
tol=10e-10; %error tolerance
ParameterField=[WaveLengthMin WaveLengthMax NumberOfWaveLength d...
    ThetaInc n1 n2 numTr tol h WaveLengthToPlot];
field1='Param';
field2='Polar';
field3='Prof';

Parameters=struct(field1,ParameterField,field2,Polarization,field3,Profile);
