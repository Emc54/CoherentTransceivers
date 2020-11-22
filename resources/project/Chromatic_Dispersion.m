%%Chromatic Dispersion modelling



% For a linear optical fiber we can model the effect of chromatic 
% dispersion on the envelope A(z,t) of a pulse by the following 
% partial differential equation:
%
% par_d(A(z,t))/par_d(z) =
% j*(D*lambda^2/(4*pi*c))*[par2_d(A(z,t))/par2_d(t)]
%
% where par_d = partial derivative, par2_d = second partial derivative,
% c=speed of light, lambda = wavelength, z=length, t=time, D = Dispersion coeff. 
%
% Taking Fourier Transform and solving the above equation:

% G(z,omega) = exp(-[j*z*omega^2*(D*lambda^2/[4*pi*c])])
%
% and setting beta_2 = D*lambda^2/(2*pi*c) we get a
% Chromatic Dispersion Transfer Function in Frequency Domain:
%
% G(z,omega) = exp(-(j*omega^2*beta_2*z)/2)
%
% where: z = length, beta_2 = group delay dispersion parameter
%
%FROM:
%Digital Signal Processing for Coherent Transceivers Employing Multilevel
%Formats, Md. Saifuddin Faruk, Member, OSA and Seb J. Savory, Fellow, IEEE,
%Fellow, OSA

%% Constants %%

no_of_symbols = 2048;
symbol_rate = 100e9; %Hz
sampling_rate =  2*symbol_rate; %Hz
time_step = 1/sampling_rate; %s
z = 50e3; %m
D = 20*10^-6; %s/m/m % Fiber dispersion in ps/nm/km (For non-dispersion-shifted fiber near 1550 nm this is typically 17.)
lambda = 1550*10^-9; %m
c = 299792458; %m/s

%% Time and Freq Vectors %%

time = linspace(0,no_of_symbols/symbol_rate,no_of_symbols);
w = linspace(-symbol_rate/2,symbol_rate/2,no_of_symbols);

symbols = randi([0 3],1,no_of_symbols);

CD_symbols = chrom_disp(z,w,D,lambda,c);

CD_symbols_2 = chrom(symbols,z,D,lambda);

function Chrom_Dispersed_Array = chrom_disp(distance,omega,dispersion_coeff,lambda,c)

    Chrom_Dispersed_Array = exp(-(1j*distance*omega.^2*(dispersion_coeff*lambda^2/(4*pi*c))));

end

function Chromatically_Dispersed_Signal = chrom(signal,distance,dispersion_coeff,lambda)

        c = 299792458;

        spectrum = fftshift(fft(signal));
        
        chrom_dispersed_spectrum = exp(-(1j*distance*spectrum.^2*(dispersion_coeff*lambda^2/(4*pi*c))));
        
        Chromatically_Dispersed_Signal = ifftshift(ifft(chrom_dispersed_spectrum));
end