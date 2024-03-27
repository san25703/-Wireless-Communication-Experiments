EsdB = -10:2:10; % Symbol energy vector in dB, default -10:2:22
Es = 10.^(EsdB/10); % Symbol energy in linear
No = 1; % Noise power in linear
SNR = Es/No; % SNR in linear
SNRdB = 10*log10(SNR); % SNR in dB
BER = zeros(length(Es),1); % initialize BER vector to store results
M = 64; % No. of subcarriers in each transmitted OFDM symbol
L = 5; % No. of channel taps
Ncp = L - 1; % No. of cyclic prefix symbols between each transmitted block
PDP_dB = [0 -3 -6 -9 -12]; % Power Delay Profile in dB, length of this vector MUST be equal to L
PDP_lin = 10.^(PDP_dB./10); % Power Delay Profile in Linear scale
ITER = 5000; % No. of iterations/channel realizations

mod_ord = 16; % Modulation Order

for ite = 1:ITER % loop over all iterations
    ite

    TxBits = randi([0,1],[1, M*log2(mod_ord)]); % generate bits equal to M*log2(modulation order), to convert to M symbols
    shaped_bits = reshape(TxBits,log2(mod_ord),M); % reshape bit vector into matrix of size (log2(mod_ord) x M)
    transposed_bits = shaped_bits.'; % transpose the matrix
    decimal_syms = bi2de(transposed_bits); % apply binary to decimal conversion
    decimal_syms = decimal_syms.'; % transpose the decimal result
    Mod_syms = qammod(decimal_syms,mod_ord); % perform QAM modulation
    avg_power = mean(abs(Mod_syms).^2); % get average power of modulated symbol vector
    Mod_syms = Mod_syms./avg_power; % divide by average power to normalize power to 1

%     h_vec = [sqrt(PDP_lin(1))*(sqrt(1/2)*(randn+1j*randn)) zeros(1,L - 1)]; % single tap rayleigh fading channel
%     h_vec = [1 zeros(1,L-1)]; % AWGN channel
    h_vec = sqrt(PDP_lin).*(sqrt(1/2)*(randn(1,L)+1j*randn(1,L))); % generate an L-tap rayleigh fading channel w/ coefficients scaled by their PDP values

    Hfreq = fft(h_vec,M); % perform M-point FFT on the time domain channel vector, note that there is NO power normalisation
    Noise = sqrt(No/2)*(randn(1,M)+1j*randn(1,M)); % generate AWGN noise

    for ix = 1:length(Es) % loop over all SNR
        Loaded_syms = sqrt(Es(ix)).*Mod_syms; % scale info symbols by the available energy per symbol
        Txsamples = sqrt(M).*ifft(Loaded_syms,M); % IFFT with power normalisation
        TxsamplesCP = [Txsamples(M-Ncp+1:M) Txsamples]; % apply cyclic prefix to the transmit symbol vector

        RxsamplesCP = conv(TxsamplesCP,h_vec); % convolution in time domain by a multi-tap fading channel
        Rxsamples = RxsamplesCP(Ncp+1:Ncp+M) + Noise; % add AWGN noise, discard the cyclic prefix & other leftover symbols from the convolution

        RxsamplesFFT = 1/sqrt(M).*fft(Rxsamples,M); % FFT with power normalisation
        Eq_syms = RxsamplesFFT./Hfreq; % perform frequency domain equalization on each subcarrier

        Eq_syms = Eq_syms./sqrt(Es(ix)); % reverse the power scaling
        Eq_syms = Eq_syms.*avg_power; % reverse the power normalisation done at the transmitter side
        Rx_syms = qamdemod(Eq_syms,mod_ord); % perform QAM demodulation to convert to integers
        Rx_syms = Rx_syms.'; % transpose the integer result
        Rx_decimal_syms = de2bi(Rx_syms,log2(mod_ord)); % apply decimal to binary conversion with number of bits equal to log2(mod_ord)
        Rx_transposed_bits = Rx_decimal_syms.'; % transpose the resultant bit matrix
        DecodedBits = reshape(Rx_transposed_bits,1, M*log2(mod_ord)); % reshape the bit matrix to a bit vector

        BER(ix) = BER(ix) + sum(DecodedBits ~= TxBits); % count all bit errors and sum them
    end
end
BER=BER/(M*log2(mod_ord)*ITER); % divide BER vector by total number of bits transmitted across all iterations

semilogy(SNRdB,BER,'^--');
grid on;
hold on;
xlabel('SNR(dB)')
ylabel('BER')
axis tight;
title('OFDM BER vs SNR (dB) U21EC080');
legend('Location', 'best')