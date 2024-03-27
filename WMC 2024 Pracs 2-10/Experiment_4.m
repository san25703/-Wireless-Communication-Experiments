clc;
clear all;
%%%%close all;

EsdB = -10:2:22; % Symbol energy vector in dB
Es = 10.^(EsdB/10); % Symbol energy in linear
No = 1; % Noise power in linear
SNR = Es/No; % SNR in linear
SNRdB = 10*log10(SNR); % SNR in dB
BER = zeros(length(Es),1); % initialize BER vector to store results
M=32; % No. of subcarriers
ITER = 10000; % No. of iterations/channel realizations

K_rice = 1; % K factor, the ratio of the power between the LOS path and all the NLOS paths, K = v^2/sigma_h^2

sigma_h = 1;
v_factor = sqrt((sigma_h^2)*K_rice); % get the amplitude of the LOS path

mod_ord = 2; % Modulation Order

for ite = 1:ITER % loop over all iterations
    ite

    TxBits = randi([0,1],[1, M*log2(mod_ord)]); % generate bits equal to M*log2(modulation order), to convert to M symbols
    shaped_bits = reshape(TxBits,log2(mod_ord),M); % reshape bit vector into matrix of size (log2(mod_ord) x M)
    transposed_bits = shaped_bits.'; % transpose the matrix
    decimal_syms = bi2de(transposed_bits); % apply binary to decimal conversion
    decimal_syms = decimal_syms.'; % transpose the decimal result
    Mod_syms = pskmod(decimal_syms,mod_ord); % perform psk modulation

%     h = 1; % to compare fading channel with simple AWGN channel, un-comment/comment this line and comment/un-comment the next ones
    mag_h = abs((sigma_h*sqrt(1/2)*(randn+1j*randn)) + v_factor); % generate the magnitude of a rician fading channel coefficient (Prob. Distribution ~Rice(v,sigma_h))
    theta_h = rand*2*pi; % generate the phase of the rician fading channel coefficient (Prob. Distribution ~U(0,2*pi))
    h = mag_h*(cos(theta_h) + 1j*sin(theta_h)); % generate a single-tap rician fading channel coefficient 

    Noise = sqrt(No/2)*(randn(1,M)+1j*randn(1,M)); % generate AWGN noise
    for ix = 1:length(Es) % loop over all SNR
        Loaded_syms = sqrt(Es(ix)).*Mod_syms; % scale info symbols by the available energy per symbol

        Txsamples = sqrt(M)*ifft(Loaded_syms); % perform ifft with power normalisation factor sqrt(M)*1/M = 1/sqrt(M)
        Rxsamples = h.*Txsamples + Noise; % multiply transmitted symbols by fading channel coefficient, and add AWGN noise
        RxsamplesFFT = 1/sqrt(M)*fft(Rxsamples); % perform ifft with power normalisation factor 1/sqrt(M)

        Eq_syms = RxsamplesFFT.*(h'/(abs(h)^2)); % equalize the received symbols on each subcarrier

        Rx_syms = pskdemod(Eq_syms,mod_ord); % perform psk demodulation to convert to integers
        Rx_syms = Rx_syms.'; % transpose the integer result
        Rx_decimal_syms = de2bi(Rx_syms,log2(mod_ord)); % apply decimal to binary conversion with number of bits equal to log2(mod_ord)
        Rx_transposed_bits = Rx_decimal_syms.'; % transpose the resultant bit matrix
        DecodedBits = reshape(Rx_transposed_bits,1, M*log2(mod_ord)); % reshape the bit matrix to a bit vector

        BER(ix) = BER(ix) + sum(DecodedBits ~= TxBits); % count all bit errors and sum them
    end
end
BER=BER/(M*log2(mod_ord)*ITER); % divide BER vector by total number of bits transmitted across all iterations

semilogy(SNRdB,BER,'b-s','linewidth',3.0,'MarkerFaceColor','b','MarkerSize',9.0);
grid on;
hold on;
axis tight;
title('OFDM BER vs SNR (dB)');