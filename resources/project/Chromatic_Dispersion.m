%%/ Chromatic Dispersion modeling /%%
%
%
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

%%/ Constants /%%

no_of_symbols = 2048;
no_of_samples = 2 * no_of_symbols;
symbol_rate = 100e9; % Baud rate
sampling_rate =  2*symbol_rate; % Hz
time_step = 1/sampling_rate; % s
z = 5e3; % m
D = 17*10^-6; % s/m/m % Fiber dispersion in ps/nm/km (For non-dispersion-shifted fiber near 1550 nm this is typically 17.)
lambda = 1550*10^-9; % m
c = 299792458; % m/s

%%/ Time and Freq Vectors /%%

symbols = pskmod(randi([0 3],1,no_of_symbols),4,pi/4,'gray');
samples = kron(symbols,ones(1,no_of_samples/no_of_symbols));

time = linspace(0, no_of_samples/sampling_rate, no_of_samples);
w = linspace(-sampling_rate/2, sampling_rate/2, no_of_samples);

cd_samples = chrom(samples,z,D,lambda,sampling_rate);

[cd_filt,trunc_cd_filt] = chrom_filt(lambda,z,D,time_step);

filtered_signal = conv2(cd_samples,cd_filt,'same');

figure
plot(time*10^9,samples);
hold on
plot(time*10^9,CD_samples);
hold off
figure
scatter(real(CD_samples),imag(CD_samples));
figure
plot(time*10^9,filtered_signal);
xlabel('time (ns)');
ylabel('signal amplitude');
figure
scatter(real(filtered_signal),imag(filtered_signal));
%legend('samples','cd samples','filtered samples');

function chromatically_dispersed_signal = chrom(analog_signal,z,D,lambda,sampling_rate)

        c = 299792458;
        
        w = linspace(-sampling_rate/2, sampling_rate/2, length(analog_signal))*2*pi;
        
        chrom_dispersion = exp(-1j*D*lambda^2*z.*(w.^2)/(4*pi*c));
        
        chrom_dispersion = ifftshift(chrom_dispersion); %Center around origin
        
        AS = fft(analog_signal); %Freq domain
       
        Chrom_Dispersed_Spectrum = chrom_dispersion.*AS; %Freq domain
       
        chromatically_dispersed_signal = ifft(Chrom_Dispersed_Spectrum); %Time domain
end


function [chromatic_dispersion_filter,trunc_filter] = chrom_filt(lambda,z,D,T)
    
    c = 299792458;

    filter_taps = 2*floor((abs(D)*lambda^2*z/(2*c*T^2))) + 1; %upper_bound
    
    %%% Use this for optimisation later if needed
    truncated_filter_taps = floor(0.6*filter_taps); %lower_bound
    
    if (mod(truncated_filter_taps,2) == 0)
        truncated_filter_taps = truncated_filter_taps + 1;
    end
    %%%
    
    k = linspace(-floor(filter_taps/2),floor(filter_taps/2),filter_taps); %arranged array from -L to L where L is half the taps
    
    k_trunc = linspace(-floor(truncated_filter_taps/2),floor(truncated_filter_taps/2),truncated_filter_taps);
    
    %taps weight calculation
    a_k_amplitude = sqrt((1j*c*T^2)/(D*lambda^2*z));  
    a_k_exp = exp((-1j*pi*c*T^2).*(k.^2)/(D*lambda^2*z));
    
    a_trunc_exp = exp((-1j*pi*c*T^2).*(k_trunc.^2)/(D*lambda^2*z));
    
    [chromatic_dispersion_filter,trunc_filter] = deal(a_k_amplitude*a_k_exp, a_k_amplitude*a_trunc_exp);
    
end
    
    
    
