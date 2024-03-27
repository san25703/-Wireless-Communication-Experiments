clc;
clear all;
%%%%close all;

EsdB = -10:5:30; % Symbol energy vector in dB, default -10:5:30
Es = 10.^(EsdB/10); % Symbol energy in linear
No = 1; % Noise power in linear
SNR = Es/No; % SNR in linear
SNRdB = 10*log10(SNR); % SNR in dB
BER = zeros(length(Es),1); % initialize BER vector to store results
M = 32; % No. of subcarriers in each transmitted OFDM symbol (Default 32)
L = 5; % No. of channel taps (Default 5)
Ncp = L - 1; % No. of cyclic prefix symbols between each transmitted block
PDP_dB = [0 -3 -6 -9 -12]; % Power Delay Profile in dB, length of this vector MUST be equal to L
PDP_lin = 10.^(PDP_dB./10); % Power Delay Profile in Linear scale
ITER = 10000; % No. of iterations/channel realizations

K_rice = 1; % K factor, the ratio of the power between the LOS path and all the NLOS paths, K = v^2/sigma_h^2
sigma_h = sqrt(PDP_lin(1)); % the average amplitude of the NLOS paths is the first value in the PDP_lin vector
v_factor = sqrt((sigma_h^2)*K_rice); % get the amplitude of the LOS path

div_ord = 2; % diversity order sets the number of times an info-symbol is repeated per block, 1 = no repeat, 2 = repeated twice per block, etc.
%%%% set div_ord to one of these values [1, 2, 4, 8] if M = 32
%%%% M/div_ord *MUST* be an integer

mod_ord = 16; % Modulation Order

for ite = 1:ITER % loop over all iterations
    ite

    Meff = M/div_ord; % effective M is actual M divided by the diversity order div_ord, Meff is the no. of actual info symbols.
    TxBits = randi([0,1],[1, Meff*log2(mod_ord)]); % generate bits equal to M*log2(modulation order), to convert to M symbols
    shaped_bits = reshape(TxBits,log2(mod_ord),Meff); % reshape bit vector into matrix of size (log2(mod_ord) x M)
    transposed_bits = shaped_bits.'; % transpose the matrix
    decimal_syms = bi2de(transposed_bits); % apply binary to decimal conversion
    decimal_syms = decimal_syms.'; % transpose the decimal result
    Mod_syms = qammod(decimal_syms,mod_ord); % perform QAM modulation
    avg_power = mean(abs(Mod_syms).^2); % get average power of modulated symbol vector
    Mod_syms = Mod_syms./avg_power; % divide by average power to normalize power to 1
    Mod_syms = repmat(Mod_syms,1,div_ord); % COPY information symbols onto multiple subcarriers

%     h_vec = [sqrt(PDP_lin(1))*(sqrt(1/2)*(randn+1j*randn)) zeros(1,L - 1)]; % single tap rayleigh fading channel
%     h_vec = [1 zeros(1,L-1)]; % AWGN channel

    h_vec = sqrt(PDP_lin).*(sqrt(1/2)*(randn(1,L)+1j*randn(1,L))); % generate an L-tap rayleigh fading channel w/ coefficients scaled by their PDP values
    mag_h = abs((sigma_h*sqrt(1/2)*(randn+1j*randn)) + v_factor); % generate the magnitude of a rician fading channel coefficient (Prob. Distribution ~Rice(v,sigma_h))
    theta_h = rand*2*pi; % generate the phase of the rician fading channel coefficient (Prob. Distribution ~U(0,2*pi))
    h_rice = mag_h*(cos(theta_h) + 1j*sin(theta_h)); % generate a single-tap rician fading channel coefficient
    h_vec(1) = h_rice; % set the first channel tap in the L-tap channel vector to the rician fading channel coefficient

    Hfreq = fft(h_vec,M); % perform M-point FFT on the time domain channel vector, note that there is NO power normalisation
    Noise = sqrt(No/2)*(randn(1,M)+1j*randn(1,M)); % generate AWGN noise

    for ix = 1:length(Es) % loop over all SNR
        Loaded_syms = sqrt(Es(ix)).*Mod_syms; % scale info symbols by the available energy per symbol
        Txsamples = sqrt(M).*ifft(Loaded_syms,M); % IFFT with power normalisation
        TxsamplesCP = [Txsamples(M-Ncp+1:M) Txsamples]; % apply cyclic prefix to the transmit symbol vector

        RxsamplesCP = conv(TxsamplesCP,h_vec); % convolution in time domain by a multi-tap fading channel
        Rxsamples = RxsamplesCP(Ncp+1:Ncp+M) + Noise; % add AWGN noise, discard the cyclic prefix & other leftover symbols from the convolution

        RxsamplesFFT = 1/sqrt(M).*fft(Rxsamples,M); % FFT with power normalisation

        %%%% Receiver MAXIMAL RATIO COMBINING (MRC) section
        Eq_syms = zeros(1,Meff); % initialize equalized symbols vector
        Rx_Mat = reshape(RxsamplesFFT,Meff,div_ord).'; % gather all the copies of each unique info symbol (copies on rows, unique symbols on columns)
        Hfreq_Mat = reshape(Hfreq,Meff,div_ord).'; % gather all the M corresponding frequency domain channel coefficients accordingly, same structure as previous line
        for rx = 1:Meff % for each *unique* actual information symbol (i.e., each column of Rx_Mat)
            MRC_vec = conj(Hfreq_Mat(:,rx))./(sum(abs(Hfreq_Mat(:,rx)).^2)); % MRC vector of a particular unique symbol, size of vector is same as number of copies (div_ord)
            Eq_syms(rx) = sum(MRC_vec.*Rx_Mat(:,rx)); % dot product of MRC vector with the received vector of a particular unique symbol, size of vector same as number of copies
        end
        %%%%

        Eq_syms = Eq_syms./sqrt(Es(ix)); % reverse the power scaling
        Eq_syms = Eq_syms.*avg_power; % reverse the power normalisation done at the transmitter side
        Rx_syms = qamdemod(Eq_syms,mod_ord); % perform QAM demodulation to convert to integers
        Rx_syms = Rx_syms.'; % transpose the integer result
        Rx_decimal_syms = de2bi(Rx_syms,log2(mod_ord)); % apply decimal to binary conversion with number of bits equal to log2(mod_ord)
        Rx_transposed_bits = Rx_decimal_syms.'; % transpose the resultant bit matrix
        DecodedBits = reshape(Rx_transposed_bits,1, Meff*log2(mod_ord)); % reshape the bit matrix to a bit vector

        BER(ix) = BER(ix) + sum(DecodedBits ~= TxBits); % count all bit errors and sum them
    end
end
BER=BER/(Meff*log2(mod_ord)*ITER); % divide BER vector by total number of bits transmitted across all iterations

semilogy(SNRdB,BER,'b-s','linewidth',3.0,'MarkerFaceColor','b','MarkerSize',9.0);
grid on;
hold on;
axis tight;
title('OFDM BER vs SNR (dB)');