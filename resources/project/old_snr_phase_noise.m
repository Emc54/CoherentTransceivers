%%%SNR and BER for Random Walk Phase Noise

%This line initialises the arrays
%%%This section is not correct needs to be commented out

[phase_noise_SNR,noisePow,num_of_errs,phase_noise_BER] = deal(zeros(sets_num,1));

figure
hold on
title("SNR vs BER for Random Walk Noise")
symbol_amplitude = abs(modulated_symbols);
noise_amplitude = abs(phase_noise_sets);

for i=1:sets_num
    phase_noise_SNR(i,:) = snr(symbol_amplitude,noise_amplitude(i,:)); %All in dB

    [num_of_errs(i,:),phase_noise_BER(i,:)] = biterr(symbols(1,:),noisy_demod_symbols(i,:),2);
end
plot(phase_noise_SNR,phase_noise_BER,"bx--")
grid on
xlabel("SNR (dB)")
ylabel("BER")
axis([-1,1,0,1])
hold off
% 
% %Since the noise is the same magnitude and only the phase component
% %increases in each case, the SNR is continuously 0dB, while the BER
% %increases.