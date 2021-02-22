
% filenames = ["BPSK CD","QPSK CD","16QAM CD","64QAM CD"];
% z_arr = (0:100:1100);
% 
%     
% for idx=1:length(filenames)
%     data = readmatrix(filenames(idx));
%     %A row containts all the data to be plotted for a specific distance
%     for row=1:size(data,1)
%         figure
%         txt=[filenames(idx),"with Z =",num2str(z_arr(row)),"km"];
%         scatter(real(data(row,:)),imag(data(row,:)));
%         legend(join(txt));
%     end
% end

z_arr = (0:10:110);

%%
for row=1:size(b_arr,1)
figure
txt=["Compensated BPSK with Z =",num2str(z_arr(row)),"km"];
scatter(real(b_arr(row,:)),imag(b_arr(row,:)));
xlabel("In-phase Amplitude");
ylabel("Quadrature Amplitude"); 
legend(join(txt));
end

figHandles = findall(0,'Type','figure');

for i=1:numel(figHandles)
saveas(figHandles(i),join(['Compensated BPSK CD ',num2str(numel(figHandles)-i+1)]),'bmp');
end
delete(findall(0));

%%
for row=1:size(q_arr,1)
figure
txt=["Compensated QPSK with Z =",num2str(z_arr(row)),"km"];
scatter(real(q_arr(row,:)),imag(q_arr(row,:))); 
xlabel("In-phase Amplitude");
ylabel("Quadrature Amplitude"); 
legend(join(txt));
end

figHandles = findall(0,'Type','figure');

for i=1:numel(figHandles)
saveas(figHandles(i),join(['Compensated QPSK CD ',num2str(numel(figHandles)-i+1)]),'bmp');
end
delete(findall(0));

%%
for row=1:size(q_arr,1)
figure
txt=["Compensated 16QAM with Z =",num2str(z_arr(row)),"km"];
scatter(real(qa_arr(row,:)),imag(qa_arr(row,:)));
xlabel("In-phase Amplitude");
ylabel("Quadrature Amplitude"); 
legend(join(txt));
end

figHandles = findall(0,'Type','figure');

for i=1:numel(figHandles)
saveas(figHandles(i),join(['Compensated 16QAM CD ',num2str(numel(figHandles)-i+1)]),'bmp');
end
delete(findall(0));

%%
for row=1:size(q_arr,1)
figure
txt=["Compensated 64QAM with Z =",num2str(z_arr(row)),"km"];
scatter(real(qam_arr(row,:)),imag(qam_arr(row,:)));
xlabel("In-phase Amplitude");
ylabel("Quadrature Amplitude"); 
legend(join(txt));
end

figHandles = findall(0,'Type','figure');

for i=1:numel(figHandles)
saveas(figHandles(i),join(['Compensated 64QAM CD ',num2str(numel(figHandles)-i+1)]),'bmp');
end
delete(findall(0));
%%
BERs = readmatrix('No CD BERs');
BPSK_BER = readmatrix('BPSK CD BER-SNR');
QPSK_BER = readmatrix('QPSK CD BER-SNR');
QAM16_BER = readmatrix('16QAM CD BER-SNR');
QAM64_BER = readmatrix('64QAM CD BER-SNR');

z = [10 100 1000]*1000;
SNR = [10 15 20 25 30 35 40];

figure
hold on
semilogy(SNR,BERs(:,1),'-bx');
semilogy(SNR,BERs(:,2),'-rx');
semilogy(SNR,BERs(:,3),'-x');
semilogy(SNR,BERs(:,4),'-gx');
legend("BPSK","QPSK","16QAM","64QAM");
title("SNR vs BER without Chromatic Dispersion");
ylabel("BER");
xlabel("SNR (dB)");

figure
hold on
semilogy(SNR,BPSK_BER(:,2),'-bx');
semilogy(SNR,QPSK_BER(:,2),'-rx');
semilogy(SNR,QAM16_BER(:,2),'-x');
semilogy(SNR,QAM64_BER(:,2),'-gx');
legend("BPSK","QPSK","16QAM","64QAM");
title("SNR vs BER with CD, z = 10km");
ylabel("BER");
xlabel("SNR (dB)");

figure
hold on
semilogy(SNR,BPSK_BER(:,3),'-bx');
semilogy(SNR,QPSK_BER(:,3),'-rx');
semilogy(SNR,QAM16_BER(:,3),'-x');
semilogy(SNR,QAM64_BER(:,3),'-gx');
legend("BPSK","QPSK","16QAM","64QAM");
title("SNR vs BER with CD, z = 100km");
ylabel("BER");
xlabel("SNR (dB)");

figure
hold on
semilogy(SNR,BPSK_BER(:,4),'-bx');
semilogy(SNR,QPSK_BER(:,4),'-rx');
semilogy(SNR,QAM16_BER(:,4),'-x');
semilogy(SNR,QAM64_BER(:,4),'-gx');
legend("BPSK","QPSK","16QAM","64QAM");
title("SNR vs BER with CD, z = 1000km");
ylabel("BER");
xlabel("SNR (dB)");
%%

BERvsSNR = readmatrix("BERvsSNR.txt");

SNR = [10 20 30 40];

BPSK_BER = BERvsSNR(1:4,2);
QPSK_BER = BERvsSNR(6:9,2);
QAM16_BER = BERvsSNR(11:14,2);
QAM64_BER = BERvsSNR(16:19,2);

figure
hold on
semilogy(SNR,BPSK_BER,'-bx');
semilogy(SNR,QPSK_BER,'-rx');
semilogy(SNR,QAM16_BER,'-x');
semilogy(SNR,QAM64_BER,'-gx');
legend("BPSK","QPSK","16QAM","64QAM");
title("SNR vs BER with CD, z = 10km");
ylabel("BER");
xlabel("SNR (dB)");

BPSK_BER = BERvsSNR(1:4,3);
QPSK_BER = BERvsSNR(6:9,3);
QAM16_BER = BERvsSNR(11:14,3);
QAM64_BER = BERvsSNR(16:19,3);

figure
hold on
semilogy(SNR,BPSK_BER,'-bx');
semilogy(SNR,QPSK_BER,'-rx');
semilogy(SNR,QAM16_BER,'-x');
semilogy(SNR,QAM64_BER,'-gx');
legend("BPSK","QPSK","16QAM","64QAM");
title("SNR vs BER with CD, z = 50km");
ylabel("BER");
xlabel("SNR (dB)");
BPSK_BER = BERvsSNR(1:4,4);
QPSK_BER = BERvsSNR(6:9,4);
QAM16_BER = BERvsSNR(11:14,4);
QAM64_BER = BERvsSNR(16:19,4);

figure
hold on
semilogy(SNR,BPSK_BER,'-bx');
semilogy(SNR,QPSK_BER,'-rx');
semilogy(SNR,QAM16_BER,'-x');
semilogy(SNR,QAM64_BER,'-gx');
legend("BPSK","QPSK","16QAM","64QAM");
title("SNR vs BER with CD, z = 100km");
ylabel("BER");
xlabel("SNR (dB)");

BPSK_BER = BERvsSNR(1:4,5);
QPSK_BER = BERvsSNR(6:9,5);
QAM16_BER = BERvsSNR(11:14,5);
QAM64_BER = BERvsSNR(16:19,5);

figure
hold on
semilogy(SNR,BPSK_BER,'-bx');
semilogy(SNR,QPSK_BER,'-rx');
semilogy(SNR,QAM16_BER,'-x');
semilogy(SNR,QAM64_BER,'-gx');
legend("BPSK","QPSK","16QAM","64QAM");
title("SNR vs BER with CD, z = 250km");
ylabel("BER");
xlabel("SNR (dB)");

BPSK_BER = BERvsSNR(1:4,6);
QPSK_BER = BERvsSNR(6:9,6);
QAM16_BER = BERvsSNR(11:14,6);
QAM64_BER = BERvsSNR(16:19,6);

figure
hold on
semilogy(SNR,BPSK_BER,'-bx');
semilogy(SNR,QPSK_BER,'-rx');
semilogy(SNR,QAM16_BER,'-x');
semilogy(SNR,QAM64_BER,'-gx');
legend("BPSK","QPSK","16QAM","64QAM");
title("SNR vs BER with CD, z = 500km");
ylabel("BER");
xlabel("SNR (dB)");
BPSK_BER = BERvsSNR(1:4,7);
QPSK_BER = BERvsSNR(6:9,7);
QAM16_BER = BERvsSNR(11:14,7);
QAM64_BER = BERvsSNR(16:19,7);

figure
hold on
semilogy(SNR,BPSK_BER,'-bx');
semilogy(SNR,QPSK_BER,'-rx');
semilogy(SNR,QAM16_BER,'-x');
semilogy(SNR,QAM64_BER,'-gx');
legend("BPSK","QPSK","16QAM","64QAM");
title("SNR vs BER with CD, z = 1000km");
ylabel("BER");
xlabel("SNR (dB)");
