function phase_corrected_signal = phase_correction(input_signal,chunk_length)
%PHASE_CORRECTION Tries to estimate and correct phase on qpsk
%   Detailed explanation goes here


if size(input_signal,2) > size(input_signal,1)
    input_signal = input_signal.';
end

%Attempt to find the mean angle in chunks
raised_phase_signal = input_signal.^4;
signal_length = length(input_signal);
chunk_size = [chunk_length,1];

raised_phase_signal = filtfilt(ones(100,1)/100,1,raised_phase_signal);


mean_filter_fun = @(theBlockStructure) mean2(theBlockStructure.data(:));
averages = blockproc(raised_phase_signal,chunk_size,mean_filter_fun);
phase_estimates = 1/4.*angle(averages);

multiple = signal_length - mod(signal_length,chunk_length);
reshaped = reshape(input_signal(1:multiple),chunk_length,[]);
leftover = input_signal(multiple+1:end);
leftover = padarray(leftover,[chunk_length-length(leftover) 0],0,'post');

combined = [reshaped leftover];

for iter=1:size(reshaped,2)
    %%% Check all four orientations
    combined(:,iter) = combined(:,iter).*exp(1j*pi/4-1j*phase_estimates(iter)+1j*3*pi/2);
end

phase_compensated = reshape(combined,[],1).';

phase_corrected_signal = phase_compensated(1:signal_length);
end

