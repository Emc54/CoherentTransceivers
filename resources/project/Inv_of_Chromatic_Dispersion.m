function chromatically_dispersed_signal = Inv_of_Chromatic_Dispersion(signal,sample_rate,D,z,l)
%CHROMATIC_DISPERSION A function that models chromatic dispersion and
%applies it to a provided input signal

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

%%%Revert to -fs/2:fs/2 approach -- this gives 238 GB of an array, too large
%%%to implement 
%%% Recheck implementation of FD_Filter
%%% Apply Reverse CD Model in case filter doesn't work as a proof of
%%% concept!!
% 
f_up = (0:(no_of_samples)/2-1);
f_down = (-(no_of_samples)/2:-1);
fs = [f_up f_down]*sample_rate/no_of_samples;
w = 2*pi*fs;

% df = sample_rate/no_of_samples;
% fs = (-sample_rate/2:df:sample_rate/2);
% w = 2*pi*fs(1:end-1);

chrom_dispersion_model = exp(1j*K.*(w.^2)*T^2);
%chrom_dispersion_model = ifftshift(chrom_dispersion_model); %Center around origin
                       
chromatically_dispersed_signal = ifft(chrom_dispersion_model.*fft((signal))); 
end

