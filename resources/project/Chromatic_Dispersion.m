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
sampling_rate =  2*symbol_rate; % Hz %Also Sampling Frequency
T = 1/sampling_rate; % s Sampling Period
z = 5e3; % m
D = 17*10^-6; % s/m/m % Fiber dispersion in ps/nm/km (For non-dispersion-shifted fiber near 1550 nm this is typically 17.)
lambda = 1550*10^-9; % m
c = 299792458; % m/s


K = (D*lambda^2*z)/(4*pi*c*T^2); 
% constant defined in #Optimal Least-Squares FIR Digital Filters
% for Compensation of Chromatic Dispersion
% in Digital Coherent Optical Receivers#

%%/ Time and Freq Vectors /%%

symbols = pskmod(randi([0 3],1,no_of_symbols),4,pi/4,'gray');
samples = kron(symbols,ones(1,no_of_samples/no_of_symbols));

time = linspace(0, no_of_samples/sampling_rate, no_of_samples);
w = linspace(-sampling_rate/2, sampling_rate/2, no_of_samples);

cd_samples = chrom(samples,z,D,lambda,sampling_rate);

[cd_filt,trunc_cd_filt] = chrom_filt(lambda,z,D,T);

filtered_signal = conv2(cd_samples,cd_filt,'same');

figure
scatter(time*10^9,samples,'bx');
hold on
scatter(time*10^9,cd_samples,'*');
xlabel('time (ns)');
ylabel('sample amplitude');
legend('samples','chromatically distorted samples');
hold off

figure
scatter(real(cd_samples),imag(cd_samples));
title('Signal with added chromatic dispersion');
xlabel('In-phase amplitude');
ylabel('Quadrature amplitude');
annotation('textbox',...
    [0.254307593307592 0.889523809523816 0.574714285714286 0.0552380952381099],...
    'String','sample rate: 200e9, z: 5e3, D: 17e-6, \lambda: 1550e-9',...
    'LineStyle','none');

figure
scatter(real(filtered_signal),imag(filtered_signal));
title('Signal after filtering chromatic dispersion');
xlabel('In-phase amplitude');
ylabel('Quadrature amplitude');
text(-1,1.4,'sample rate: 200e9, z: 5e3, D: 17e-6, \lambda: 1550e-9'); 

%%/ LS CD Compensation Filter /%%

N_c = 1; %Has to be odd and less than N
E_min = 100;
Nc_best = 0;

for i = 1:length(cd_filt)/2+1

    least_squares_filter = LSfilt(samples,K,T,N_c)';

    E = sum(abs(cd_filt(1,(length(cd_filt)-N_c)/2+1:(length(cd_filt)+N_c)/2)-least_squares_filter).^2)/2/pi;

    if E < E_min
        E_min = E;
        Nc_best = N_c;
    end
    N_c = N_c + 2;
end
           

function chromatically_dispersed_signal = chrom(analog_signal,z,D,lambda,sampling_rate)

        c = 299792458;
        
        %sampling frequency F = 1/T
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
    
function h_hat = LSfilt(analog_signal,K,T,N_c)

%        N = 2*floor(2*K*pi) + 1; %No of Filter Taps
         
         n = linspace(-floor(N_c/2),floor(N_c/2),N_c);
         
         wT = linspace(-(1/T)/2, (1/T)/2, length(analog_signal))*2*pi;
         
%         h = sqrt(1j/(4*K*pi))*exp((-1j.*(n.^2))/(4*K));
%            
%          H = ones(1,length(n));
%          
%          for k = 1:length(n)
%             for i = 1:N_c
%                 H(k) = H(k) + h(i+(N-N_c)/2)*exp(-1j*n(i)*wT(k)); 
%                 % This line gets values in h from N-N_c/2 to N+N_c/2 using i  
%                 % It iterates through all of wT and n using k
%             end
%          end
         
         D = nan(1,length(n));
         
         for i = 1:length(n)
             D(i) = exp(-1j*((n(i)^2)/(4*K)+3*pi/4))/(4*sqrt(pi*K))...
                 *(...
                 erfi(exp(3j*pi/4)*(2*K*pi-n(i)/(2*sqrt(K))))...
                 +erfi(exp(3j*pi/4)*(2*K*pi+n(i)/(2*sqrt(K))))...
                 );
         end
         
         Q = eye(N_c);
%          
%          for i=1:N_c
%              if i == (N_c-1)/2
%                  Q(1,i) = (wT(1,end) - wT(1))/(2*pi);
%              else
%                  Q(1,i) = (exp(-1j*i*wT(1)) - exp(-1j*i*wT(end)))/(2*1j*pi*i);
%              end
%          end
%          
         for i=1:N_c
             if i == (N_c-1)/2
                 Q(1,i) = 1;
             else
                 Q(1,i) = (exp(-1j*i*(-pi)) - exp(-1j*i*pi))/(2*1j*pi*i);
             end
         end
         
         h_hat = Q\D';
         
end
    
    


