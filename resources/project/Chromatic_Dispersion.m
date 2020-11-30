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

%% Constants %%

no_of_symbols = 2048;
no_of_samples = 8 * no_of_symbols;
symbol_rate = 10e9; % Baud rate
fs = 8*symbol_rate; % Hz %Sampling Rate / Sampling Frequency
T = 1/fs; % s Sampling Period
z = 5000e3; % m
D = 17*10^-6; % s/m/m % Fiber dispersion in ps/nm/km (For non-dispersion-shifted fiber near 1550 nm this is typically 17.)
lambda = 1550*10^-9; % m
c = 299792458; % m/s
K = (D*lambda^2*z)/(4*pi*c*T^2); 
% constant defined in #Optimal Least-Squares FIR Digital Filters
% for Compensation of Chromatic Dispersion
% in Digital Coherent Optical Receivers#
N_c = 32;

%% Signal Vectors %%

symbols = pskmod(randi([0 3],1,no_of_symbols),4,pi/4,'gray');

%samples = kron(symbols,ones(1,no_of_samples/no_of_symbols));

time = linspace(0, no_of_samples/fs, no_of_samples);
 
%Testing for pulse shaping
tx_ps_filter = comm.RaisedCosineTransmitFilter("FilterSpanInSymbols",10,"OutputSamplesPerSymbol",8,"RolloffFactor",0.5,"Shape","Square root","Gain",1);
rx_ps_filter = comm.RaisedCosineReceiveFilter("FilterSpanInSymbols",10,"InputSamplesPerSymbol",8,"RolloffFactor",0.5,"Shape","Square root","Gain",1);

rc_samples =  tx_ps_filter(symbols);
rc_samples = reshape(rc_samples,1,[]);

%FFT frequency array should start at 0 go slightly below N/2 ((N-1)/2) then wrap to -N/2 and go to
%zero again
f_up = (1:(no_of_samples)/2);
f_down = (-(no_of_samples-1)/2:0);
f = [f_up f_down]*fs/no_of_samples;
w = 2*pi*f;
       
%% Modelling Chromatic Dispersion %%

chrom_dispersion_model = exp(-1j*K.*(w.^2));
        
chrom_dispersion_model = ifftshift(chrom_dispersion_model); %Center around origin
        
analog_signal_spectrum = fft(rc_samples); %Freq domain
        
chrom_dispersed_Spectrum = chrom_dispersion_model.*analog_signal_spectrum; %Freq domain
       
chromatically_dispersed_signal = ifft(chrom_dispersed_Spectrum); %Time domain


%% Modelling the Simple Chromatic Dispersion Filter %%

N = 2*floor((abs(D)*lambda^2*z/(2*c*T^2))) + 1; % Number of Filter Taps (upper bound) 
     
k = linspace(-floor(N/2),floor(N/2),N); %arranged array from -L to L where L is N/2
    
%taps weight calculation
a_k_amplitude = sqrt((1j*c*T^2)/(D*lambda^2*z));  
a_k_exp = exp((-1j*pi*c*T^2).*(k.^2)/(D*lambda^2*z));
   
chromatic_dispersion_filter = a_k_amplitude*a_k_exp;   %% This is h(n)

%fvtool(chromatic_dispersion_filter);

%% Modelling the Least Squares Filter %%

% Assumes full bandwidth from -pi to pi, can be reduced

n = linspace(-floor(N/2),floor(N/2),N);
         
D = exp(-1j*((n.^2)/(4*K)+3*pi/4))/(4*sqrt(pi*K)).*(...
       erfi(exp(3j*pi/4)*(2*K*pi-n/(2*sqrt(K))))...
       +erfi(exp(3j*pi/4)*(2*K*pi+n/(2*sqrt(K))))...
       );
        
%% Modelling the Adaptive NLMS Filter %%

lms = dsp.LMSFilter('Length',N_c,'StepSize',0.08,'Method','Normalized LMS');


%% Filtering the Signals %%

cd_filtered_signal = fftfilt(chromatic_dispersion_filter,chromatically_dispersed_signal);

ls_filtered_signal = fftfilt(D,chromatically_dispersed_signal);

[lms_filtered_signal,lms_error,lms_weights] = lms(chromatically_dispersed_signal',rc_samples');

%% Figures %%
figure
scatter(real(symbols),imag(symbols),'x');
hold on
scatter(real(rc_samples),imag(rc_samples),'+');
grid on
title('Signal');
xlabel('In-phase amplitude');
ylabel('Quadrature amplitude');

figure
scatter(real(rc_samples),imag(rc_samples),'+');
grid on
title('Signal');
xlabel('In-phase amplitude');
ylabel('Quadrature amplitude');


figure
scatter(time*10^9,rc_samples,'bx');
hold on
scatter(time*10^9,chromatically_dispersed_signal,'*');
xlabel('time (ns)');
ylabel('sample amplitude');
legend('samples','chromatically distorted samples');
hold off

figure
scatter(real(chromatically_dispersed_signal),imag(chromatically_dispersed_signal));
title('Signal with added chromatic dispersion');
xlabel('In-phase amplitude');
ylabel('Quadrature amplitude');
annotation('textbox',...
    [0.254307593307592 0.889523809523816 0.574714285714286 0.0552380952381099],...
    'String','sample rate: 200e9, z: 5000e3, D: 17e-6, \lambda: 1550e-9',...
    'LineStyle','none');

figure
scatter(real(cd_filtered_signal),imag(cd_filtered_signal));
title('Signal after simple filtering of chromatic dispersion');
xlabel('In-phase amplitude');
ylabel('Quadrature amplitude');
%text(-1,1.4,'sample rate: 200e9, z: 5e3, D: 17e-6, \lambda: 1550e-9'); 


figure
scatter(real(ls_filtered_signal),imag(ls_filtered_signal));
title('Signal after filtering with LS filter');
xlabel('In-phase amplitude');
ylabel('Quadrature amplitude');

figure
scatter(real(lms_filtered_signal),imag(lms_filtered_signal));
title('Signal after filtering with adaptive NLMS');
xlabel('In-phase amplitude');
ylabel('Quadrature amplitude');

%% Demodulation %%
x = reshape(cd_filtered_signal,8,[]);
matched_symbols1  = rx_ps_filter(reshape(cd_filtered_signal,no_of_samples/no_of_symbols,[]));
matched_symbols2 = rx_ps_filter(reshape(ls_filtered_signal,no_of_samples/no_of_symbols,[]));
matched_symbols3 = rx_ps_filter(reshape(lms_filtered_signal,no_of_samples/no_of_symbols,[]));

recovered1 = pskdemod(matched_symbols1,4,pi/4,'gray');
recovered2 = pskdemod(matched_samples2,4,pi/4,'gray');
recovered3 = pskdemod(matched_samples3,4,pi/4,'gray');


%% Deprecated %%
% function chromatically_dispersed_signal = chrom(analog_signal,K,sampling_rate)
% 
% 
%         %TO FIX
%         w = linspace(-(sampling_rate)/2, (sampling_rate)/2, length(analog_signal))*2*pi;
%         
%         %c = 299792458;
%         %wT = w*(1/fs);
%         
%         %K = (D*lambda^2*z)/(4*pi*c*(1/fs)^2); 
%         chrom_dispersion = exp(-1j*K.*(w.^2));
%         
%         %analog_signal = [analog_signal ,0];
%         chrom_dispersion = ifftshift(chrom_dispersion); %Center around origin
%         
%         AS = fft(analog_signal); %Freq domain
%         
%         Chrom_Dispersed_Spectrum = chrom_dispersion.*AS; %Freq domain
%        
%         chromatically_dispersed_signal = ifft(Chrom_Dispersed_Spectrum); %Time domain
% end


% function chromatic_dispersion_filter = chrom_filt(lambda,z,D,T)
%     
%     c = 299792458;
% 
%     filter_taps = 2*floor((abs(D)*lambda^2*z/(2*c*T^2))) + 1; %upper_bound
%      
%     k = linspace(-floor(filter_taps/2),floor(filter_taps/2),filter_taps); %arranged array from -L to L where L is half the taps
%     
%     %taps weight calculation
%     a_k_amplitude = sqrt((1j*c*T^2)/(D*lambda^2*z));  
%     a_k_exp = exp((-1j*pi*c*T^2).*(k.^2)/(D*lambda^2*z));
%    
%     chromatic_dispersion_filter = a_k_amplitude*a_k_exp;
%     
% end
    
% function h_hat = LSfilt(K,N_c)
% 
%          N = 2*floor(2*K*pi) + 1; %No of Filter Taps
%          
%          if N_c < N
%              N = N_c;
%          end
%          
%          n = linspace(-floor(N/2),floor(N/2),N);
%          
%          D = exp(-1j*((n.^2)/(4*K)+3*pi/4))/(4*sqrt(pi*K)).*(...
%                 erfi(exp(3j*pi/4)*(2*K*pi-n/(2*sqrt(K))))...
%                 +erfi(exp(3j*pi/4)*(2*K*pi+n/(2*sqrt(K))))...
%                 );
%          
%          h_hat = D';
%          
% end
