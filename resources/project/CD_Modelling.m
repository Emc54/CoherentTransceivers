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


%%% Delay is causing issues account with autocorrelation
%%plot(real(ifft(conj(fft(symbols)).*(fft(demoded_symbols)))))

no_of_symbols = 2048;
samples_per_symbol = 8;
no_of_samples = samples_per_symbol * no_of_symbols;
symbol_rate = 10e9; % Baud rate
fs = 8*symbol_rate; % Hz %Sampling Rate / Sampling Frequency
T = 1/fs; % s Sampling Period
z = 5000e3; % m
D = 34*10^-6; % s/m/m % Fiber dispersion in ps/nm/km (For non-dispersion-shifted fiber near 1550 nm this is typically 17.)
lambda = 1550*10^-9; % m
c = 299792458; % m/s
K = (D*lambda^2*z)/(4*pi*c*T^2); 
% constant defined in #Optimal Least-Squares FIR Digital Filters
% for Compensation of Chromatic Dispersion
% in Digital Coherent Optical Receivers#
Taps = 32;
N = 2*floor(2*K*pi) + 1;

%% Signal Vectors %%

symbols = pskmod(randi([0 3],1,no_of_symbols),4,pi/4,'gray');

samples = kron(symbols,ones(1,no_of_samples/no_of_symbols));

time = linspace(0, no_of_samples/fs, no_of_samples);
 
%Testing for pulse shaping
tx_ps_filter = comm.RaisedCosineTransmitFilter("FilterSpanInSymbols",32,"OutputSamplesPerSymbol",samples_per_symbol,"RolloffFactor",0.5,"Shape","Square root","Gain",1);
rx_ps_filter = comm.RaisedCosineReceiveFilter("FilterSpanInSymbols",32,"InputSamplesPerSymbol",samples_per_symbol,"RolloffFactor",0.5,"Shape","Square root","Gain",1);
specscope = dsp.SpectrumAnalyzer('SampleRate',fs);

rc_samples = tx_ps_filter(symbols')';

%eyediagram(rc_samples,2*samples_per_symbol);
%specscope(rc_samples');

%FFT frequency array should start at 0 go slightly below N/2 ((N-1)/2) then wrap to -N/2 and go to
%zero again
f_up = (1:(no_of_samples)/2);
f_down = (-(no_of_samples)/2+1:0);
f = [f_down f_up]*fs/no_of_samples;
w = 2*pi*f;
% wT = w*T;   %This should be from -pi to pi     
%% Modelling Chromatic Dispersion %%

chrom_dispersion_model = exp(-1j*K.*(w.^2)*T^2);
fvtool(chrom_dispersion_model);
%chrom_dispersion_model = ifftshift(chrom_dispersion_model); %Center around origin
        
analog_signal_spectrum = fftshift(fft(rc_samples)); %Freq domain and shift to -freq to 0 to +freq
        
chrom_dispersed_Spectrum = chrom_dispersion_model.*analog_signal_spectrum; %Freq domain
       
chromatically_dispersed_signal = (ifft(chrom_dispersed_Spectrum)); %Time domain

%%%ATTEMPT #2%%%
g_amp = sqrt(c/(1j*D*lambda^2*z));
g_phase = exp(1j*pi*c.*(time(1:N).^2)/(D*lambda^2*z));
g_impulse = g_amp*g_phase;

chrom_signal = conv(rc_samples,g_impulse,"same");

% figure
% scatter(real(chrom_signal),imag(chrom_signal),'x');
% grid on
% title('time domain attempt at chromatic dispersion');
% xlabel('In-phase amplitude');
% ylabel('Quadrature amplitude');


%% Modelling the Simple Chromatic Dispersion Filter %%

N = 2*floor(2*K*pi) + 1; % Number of Filter Taps (upper bound) 
     
n = linspace(-floor(N/2),floor(N/2),N); %arranged array from -L to L where L is N/2
    
%taps weight calculation
hn_amplitude = sqrt(1j/(4*K*pi));  
hn_exp = exp(-1j.*(n.^2)/(4*K));
chromatic_dispersion_filter = hn_amplitude*hn_exp;   %% This is h(n)

fvtool(chromatic_dispersion_filter);

%% Modelling the Least Squares Filter %%
N_c = N;
n2 = linspace(-floor(N_c/2),floor(N_c/2),N_c);
% Assumes full bandwidth from -pi to pi, can be reduced
         
D_hat = exp(-1j.*((n.^2)/(4*K)+3*pi/4))/(4*sqrt(pi*K)).*(...
       erfi(exp(3j*pi/4)*(2*K*pi-n/(2*sqrt(K))))...
       +erfi(exp(3j*pi/4)*(2*K*pi+n/(2*sqrt(K))))...
       );
   
D_n1 = exp(-1j*3*pi/4)*exp(-1j.*(n.^2)/(4*K))/(4*sqrt(pi*K));   
D_n2 = erfi((exp(1j*3*pi/4)*2*K*pi/(2*sqrt(K)))-(exp(1j*3*pi/4)/(2*sqrt(K))).*n);
D_n3 = erfi((exp(1j*3*pi/4)*2*K*pi/(2*sqrt(K)))+(exp(1j*3*pi/4)/(2*sqrt(K))).*n);
D_n = D_n1.*(D_n2+D_n3); 
%% Modelling the Adaptive NLMS Filter %%

lms = dsp.LMSFilter('Length',32,'StepSize',0.08,'Method','Normalized LMS');


%% Filtering the Signals %%

cd_filtered_signal = fftfilt(chromatic_dispersion_filter,chromatically_dispersed_signal);

ls_filtered_signal = fftfilt(D_n,chromatically_dispersed_signal);

[lms_filtered_signal,lms_error,lms_weights] = lms(chromatically_dispersed_signal',rc_samples');

%The last filter autocorrelates: Freq domain -> power spectrum -> fft back to time domain

%% Figures %%
% %eyediagram(cd_filtered_signal,2*samples_per_symbol);
% 
% figure
% scatter(real(symbols),imag(symbols),'x');
% hold on
% scatter(real(rc_samples),imag(rc_samples),'+');
% grid on
% title('Signal before and after pulse shaping');
% xlabel('In-phase amplitude');
% ylabel('Quadrature amplitude');
% legend('Modded Symbols','Pulse Shaped Symbols');
% 
% figure
% scatter(real(rc_samples),imag(rc_samples),'+');
% grid on
% title('Pulse Shaped Signal');
% xlabel('In-phase amplitude');
% ylabel('Quadrature amplitude');
% 
% 
% figure
% scatter(time*10^9,rc_samples,'bx');
% hold on
% scatter(time*10^9,chromatically_dispersed_signal,'*');
% title('Signal over time');
% xlabel('time (ns)');
% ylabel('sample amplitude');
% legend('samples','chromatically distorted samples');
% hold off
% 
% figure
% scatter(real(chromatically_dispersed_signal),imag(chromatically_dispersed_signal));
% title('Signal with added chromatic dispersion');
% xlabel('In-phase amplitude');
% ylabel('Quadrature amplitude');
% annotation('textbox',...
%     [0.254307593307592 0.889523809523816 0.574714285714286 0.0552380952381099],...
%     'String','sample rate: 200e9, z: 5000e3, D: 17e-6, \lambda: 1550e-9',...
%     'LineStyle','none');
% 
% figure
% scatter(real(cd_filtered_signal),imag(cd_filtered_signal));
% title('Signal after simple filtering of chromatic dispersion');
% xlabel('In-phase amplitude');
% ylabel('Quadrature amplitude');
% %text(-1,1.4,'sample rate: 200e9, z: 5e3, D: 17e-6, \lambda: 1550e-9'); 
% 
% 
% figure
% scatter(real(ls_filtered_signal),imag(ls_filtered_signal));
% title('Signal after filtering with LS filter');
% xlabel('In-phase amplitude');
% ylabel('Quadrature amplitude');
% 
% figure
% scatter(real(lms_filtered_signal),imag(lms_filtered_signal));
% title('Signal after filtering with adaptive NLMS');
% xlabel('In-phase amplitude');
% ylabel('Quadrature amplitude');

%% Demodulation %%

matched_symbols0 = rx_ps_filter(cd_filtered_signal')';
matched_symbols1 = rx_ps_filter(rc_samples')';
matched_symbols2 = rx_ps_filter(ls_filtered_signal')';
matched_symbols3 = rx_ps_filter(lms_filtered_signal)';

recovered0 = pskdemod(matched_symbols0,4,pi/4,'gray');
recovered1 = pskdemod(matched_symbols1,4,pi/4,'gray');
recovered2 = pskdemod(matched_symbols2,4,pi/4,'gray');
recovered3 = pskdemod(matched_symbols3,4,pi/4,'gray');

rec_symbols = pskdemod(symbols,4,pi/4,'gray');
rec_symbols2 = rec_symbols;
rec_symbols3 = rec_symbols;
%Shift by number of Taps, when taps known, also use autocorellation to get
%estimate when taps unknown or the channel has further shifted the signal
autocorr = real(ifft(conj(fft(rec_symbols)).*(fft(recovered1)))); %manual estimation

[correlation,lags] = xcorr(rec_symbols,recovered0); %Matlab defined function
[max_corr,idx_corr] = max(correlation);
lag = lags(idx_corr);

%plot(real(ifft(conj(fft(rec_symbols)).*(fft(recovered1)))));
plot(lags,correlation);

recovered0 = recovered0(1-lag:end);
recovered1 = recovered1((1+Taps):end);
recovered2 = recovered2((1+Taps):end);
rec_symbols = rec_symbols(1:end+lag);
rec_symbols3 = rec_symbols3(1:end-Taps);

[~,BER0] = biterr(rec_symbols,recovered0,2);
[~,BER1] = biterr(rec_symbols3,recovered1,2);

plot(rec_symbols3,recovered1);

[~,BER2] = biterr(rec_symbols3,recovered2,2);
[~,BER3] = biterr(rec_symbols2,recovered3,2);

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
