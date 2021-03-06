function CD_signal = Chromatic_Dispersion(signal,num_symbols,symbol_rate,D,z,l)
%CHROMATIC_DISPERSION A function that models chromatic dispersion and
%applies it to a provided input signal

no_of_samples = length(signal);
%T = 1/sample_rate; % s Sampling Period
lambda = l; %1550*10^-9; % m
c = 3e8; % m/s
%z = 5000e3; % m
%D = 16*10^-6; % s/m/m % Fiber dispersion in ps/nm/km (For non-dispersion-shifted fiber near 1550 nm this is typically 17.)

%FFT frequency array should start at 0 go slightly below N/2 ((N-1)/2) then wrap to -N/2 and go to
%zero again


f_up = (0:(no_of_samples)/2-1);
f_down = (-(no_of_samples)/2:-1);
fs = [f_up f_down]*symbol_rate/num_symbols;
w = 2*pi*fs;


 
%         %%%OLD APPROACH%%%
% % constant defined in #Optimal Least-Squares FIR Digital Filters
% % for Compensation of Chromatic Dispersion
% % in Digital Coherent Optical Receivers#
% K = (D*lambda^2*z)/(4*pi*c*T^2); 
% chrom_dispersion_model = exp(-1j*K.*(w.^2)*T^2);
% CD_signal = ifft(chrom_dispersion_model.*fft(signal)); 
%        %%%END OLD APPROACH%%%

beta_2 = -lambda*lambda*D/(2*pi*c);

CD_signal = ifft(fft(signal).*exp(-1j*beta_2/2*w.^2*z));

end

