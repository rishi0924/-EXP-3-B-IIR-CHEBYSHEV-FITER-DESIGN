# EXP 3 : IIR-CHEBYSHEV-FITER-DESIGN

## AIM: 

 To design an IIR Chebyshev filter  using SCILAB. 

## APPARATUS REQUIRED: 
PC installed with SCILAB. 

## PROGRAM (LPF): 
```python
clc ;
close ;
wp=input('Enter the pass band frequency (Radians )= ' );
ws=input('Enter the stop band frequency (Radians )= ' );
alphap=input( ' Enter the pass band attenuation (dB)=' );
alphas=input( ' Enter the stop band attenuation(dB)=' );
T=input('Enter the Value of sampling Time=');
//Pre warping- Bilinear Transformation
omegap=(2/T)*tan(wp/2);
disp(omegap,'omegap=');
omegas=(2/T)*tan(ws/2);
disp(omegas,'omegas=');
//Order of the filter
N=acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1)))/(acosh(omegas/omegap));
disp(N,'N=');
N=ceil(N);
disp(N,'Round off value of N=');
//Cut off frequency
omegac=omegap/(((10^(0.1*alphap)) -1)^(1/(2* N)));
disp(omegac,'omegac=');
Epsilon = sqrt ((10^(0.1*alphap))-1);
disp(Epsilon,'Epsilon=');
[pols ,gn] = zpch1(N, Epsilon,omegap );
disp(gn,'Gain');
disp(pols,'Poles');
hs=poly(gn,'s','coeff')/real(poly(pols,'s'));
disp(hs,'Analog Low pass Chebyshev Filter Transfer function');
z=poly(0,'z');//Defining variable z
Hz=horner(hs,(2/ T)*((z -1)/(z+1)))// Bilinear Transformation
disp(Hz,'Digital LPF Transfer function H(Z)=');
HW=frmag(Hz,512); // Frequency response
w=0:%pi/511:%pi ;
plot(w/%pi,abs(HW));
xlabel(' Normalized Digital Frequency w');
ylabel('Magnitude ');
title(' Frequency Response of Chebyshev IIR LPF');
```
## PROGRAM (HPF): 
```python
clc;
close;

// ============================
// User Inputs
// ============================
wp = input('Enter the pass band frequency (Radians )= ');  // Passband edge > Stopband edge
ws = input('Enter the stop band frequency (Radians )= ');
alphap = input('Enter the pass band attenuation (dB)= ');
alphas = input('Enter the stop band attenuation (dB)= ');
T = input('Enter the Value of sampling Time= ');

// ============================
// Pre-warping (Bilinear Transformation)
// ============================
omegap = (2/T)*tan(wp/2);   // Passband edge (analog)
disp(omegap,'omegap=');
omegas = (2/T)*tan(ws/2);   // Stopband edge (analog)
disp(omegas,'omegas=');

// ============================
// Order of the HPF
// ============================
// For HPF: use acosh(omegap/omegas)
N = acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1))) / acosh(omegap/omegas);
disp(N,'N=');
N = ceil(N);
disp(N,'Round off value of N=');

// ============================
// Cutoff frequency
// ============================
omegac = omegap / (((10^(0.1*alphap))-1)^(1/(2*N)));
disp(omegac,'omegac=');

// ============================
// Chebyshev Prototype (Normalized LPF)
// ============================
Epsilon = sqrt((10^(0.1*alphap))-1);
disp(Epsilon,'Epsilon=');

[pols, gn] = zpch1(N, Epsilon, 1);   // normalized LPF prototype at 1 rad/s
disp(gn,'Gain');
disp(pols,'Poles');

s = poly(0,'s');   // Laplace variable
hs = poly(gn,'s','coeff') / real(poly(pols,'s'));
disp(hs,'Analog Normalized Chebyshev LPF');

// ============================
// LPF → HPF Transformation (s -> omegac/s)
// ============================
sh = horner(hs, omegac/s);
disp(sh,'Analog Chebyshev High-Pass Filter');

// ============================
// Bilinear Transformation: s -> (2/T)*((z-1)/(z+1))
// ============================
z = poly(0,'z');   // z-domain variable
Hz = horner(sh, (2/T) * ((z-1)/(z+1)));
disp(Hz,'Digital HPF Transfer function H(Z)=');

// ============================
// Frequency Response
// ============================
HW = frmag(Hz, 512);
w = 0:%pi/511:%pi;

plot(w/%pi, abs(HW));
xlabel('Normalized Digital Frequency ×π rad/sample');
ylabel('Magnitude');
title('Frequency Response of Chebyshev IIR High-Pass Filter');
```
## OUTPUT (LPF) : 
<img width="1919" height="1079" alt="image" src="https://github.com/user-attachments/assets/533de03d-dff7-4876-9de3-88d00fd84ecb" />

<img width="1882" height="1068" alt="image" src="https://github.com/user-attachments/assets/7c6c949b-d150-4616-aaa4-f4ce57e5d924" />

## OUTPUT (HPF) : 
<img width="1919" height="1077" alt="image" src="https://github.com/user-attachments/assets/d2413a83-14ca-4610-a271-4afe055799c4" />

<img width="1919" height="1079" alt="image" src="https://github.com/user-attachments/assets/b1f27989-c16e-429d-b175-e1f3e2021a71" />

## RESULT: 
IIR-CHEBYSHEV-FITER-DESIGN OF LPF AND HPF IS SUCCESSFULLY EXECUTED.
