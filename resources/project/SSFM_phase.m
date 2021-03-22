%% Constants
num_symbols = 32*1024; 
symbol_rate = 32e9;
time_step = 1/symbol_rate;
upsampling = 16;
sample_rate = upsampling*symbol_rate; %optical sample rate
num_samples = upsampling*num_symbols; %optical samples
Ts = 1/sample_rate; %optical period
rng(1960);
span = 20;
upsampling = 16;
Rolloff = 0.1;
receiver_samples_per_symbol = 2; 
D = 17;                     % ps/nm/km
beta_2 = -D/0.784;           % ps^2/km
beta_2 = beta_2*10^-27;
D = D*10^-6;
lambda = 1550e-9;
z = 40e3;
f_up = (0:num_samples/2-1);
f_down = (-num_samples/2:-1);
fs = [f_up f_down]*symbol_rate/num_symbols;
w = 2*pi*fs.';


%%  


symbols_q = randi([0 3],1,num_symbols);
qpsk_mod = pskmod(symbols_q,4,pi/4,'gray');


tx_ps_filter = comm.RaisedCosineTransmitFilter("FilterSpanInSymbols",span,"OutputSamplesPerSymbol",upsampling,"RolloffFactor",Rolloff,"Shape","Square root","Gain",1);
rx_ps_filter = comm.RaisedCosineReceiveFilter("FilterSpanInSymbols",span,"InputSamplesPerSymbol",upsampling,"DecimationFactor",upsampling/receiver_samples_per_symbol,"RolloffFactor",Rolloff,"Shape","Square root","Gain",1);
    
upsampled_qpsk = tx_ps_filter(qpsk_mod.').';

% cd_signal = Chromatic_Dispersion(upsampled_qpsk,num_symbols,symbol_rate,D,z,lambda);

num_steps = z/dz;
x = upsampled_qpsk.';

dz=1000;
alpha = 0;
gamma = 0;

for i=1:num_steps

    eDh = exp(-(1j*beta_2.*w.*w + alpha)*dz/2);
    
    eNh = exp(-gamma*1j*dz.*(abs(x).^2));
    
    x = ifft(eDh.*fft(eNh.*x));
    
end
cd_signal = x.';

downsampled_qpsk = rx_ps_filter(cd_signal.').';

CD_Static_Filter = CD_Filter(symbol_rate,D,z,lambda);

downsampled_qpsk = fftfilt(CD_Static_Filter,downsampled_qpsk);
scatterplot(downsampled_qpsk,2);

if mean(abs(downsampled_qpsk(2:2:end)).^2)>mean(abs(downsampled_qpsk(1:2:end)).^2)
     moded_qpsk=downsampled_qpsk(2:2:end);
 else
     moded_qpsk=downsampled_qpsk(1:2:end);
end


chunk_length = 32768;
moded_qpsk = phase_correction(moded_qpsk.',chunk_length);
scatterplot(moded_qpsk(50:end));

 lag = finddelay(qpsk_mod,moded_qpsk);
 
 if lag > 0
    moded_qpsk = moded_qpsk(1+lag:end);
    symbols_q = symbols_q(1:end-lag);
else
    lag_q = abs(lag_q);
    moded_qpsk = moded_qpsk(1:end-lag);
    symbols_q = symbols_q(1+lag:end);
 end
  
qpsk = pskdemod(moded_qpsk,4,pi/4,'gray');

[errors,BER] = biterr(symbols_q,qpsk,log2(4))
    