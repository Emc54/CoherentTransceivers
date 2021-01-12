function chromatically_dispersed_signal = Chromatic_Dispersion(signal,sample_rate,D,z,l)
%CHROMATIC_DISPERSION A function that models chromatic dispersion and
%applies it to a provided input signal

%   Hard-coded parameters instead of inputs for easier use,
%   lambda=1550e-9, z=5000e3, D=34e-6

no_of_samples = length(signal);
T = 1/sample_rate; % s Sampling Period
%z = 5000e3; % m
%D = 16*10^-6; % s/m/m % Fiber dispersion in ps/nm/km (For non-dispersion-shifted fiber near 1550 nm this is typically 17.)
lambda = l; %1550*10^-9; % m
c = 299792458; % m/s
K = (D*lambda^2*z)/(4*pi*c*T^2); 
% constant defined in #Optimal Least-Squares FIR Digital Filters
% for Compensation of Chromatic Dispersion
% in Digital Coherent Optical Receivers#

%FFT frequency array should start at 0 go slightly below N/2 ((N-1)/2) then wrap to -N/2 and go to
%zero again
f_up = (1:(no_of_samples)/2);
f_down = (-(no_of_samples)/2+1:0);
f = [f_down f_up]*sample_rate/no_of_samples;
w = 2*pi*f;

chrom_dispersion_model = exp(-1j*K.*(w.^2)*T^2);
%chrom_dispersion_model = ifftshift(chrom_dispersion_model); %Center around origin
                       
chromatically_dispersed_signal = ifft((chrom_dispersion_model.*fftshift(fft(signal)))); 
end

