EsdB = -10:2:22; % Symbol energy vector in dB
Es = 10.^(EsdB/10); % Symbol energy in linear
No = 1; % Noise power in linear
SNR = Es/No; % SNR in linear
SNRdB = 10*log10(SNR); % SNR in dB
BER = zeros(length(Es),1); % initialize BER vector to store results
M=64; % No. of time domain symbols in each transmitted block
ITER = 5000; % No. of iterations/channel realizations

N_pil = 1; % No. of transmitted pilot symbols in one block, must be < M, but >= 1. Increase this to reduce BER (1, 4, 16).

sigma_h_squared = 1; % average channel gain/variance
sigma_h = sqrt(sigma_h_squared); % channel amplitude/standard deviation

mod_ord = 16; % Modulation Order

for ite = 1:ITER % loop over all iterations
    ite

    Meff = M - N_pil; % effective M is actual M subtracted by no. of pilot symbols transmitted per block
    TxBits = randi([0,1],[1, Meff*log2(mod_ord)]); % generate bits equal to Meff*log2(modulation order), to convert to Meff symbols
    shaped_bits = reshape(TxBits,log2(mod_ord),Meff); % reshape bit vector into matrix of size (log2(mod_ord) x Meff)
    transposed_bits = shaped_bits.'; % transpose the matrix
    decimal_syms = bi2de(transposed_bits); % apply binary to decimal conversion
    decimal_syms = decimal_syms.'; % transpose the decimal result
    Mod_syms = qammod(decimal_syms,mod_ord); % perform QAM modulation
    avg_power = mean(abs(Mod_syms).^2); % get average power of modulated symbol vector
    Mod_syms = Mod_syms./avg_power; % divide by average power to normalize power to 1

    h = sigma_h*sqrt(1/2)*(randn+1j*randn); % generate a single-tap rayleigh fading channel coefficient (Prob. Distribution ~CN(0,sigma_h^2))

    Noise = sqrt(No/2)*(randn(1,M)+1j*randn(1,M)); % generate AWGN noise
    for ix = 1:length(Es) % loop over all SNR
        Loaded_syms = sqrt(Es(ix)).*Mod_syms; % scale info symbols by the available energy per symbol
        pil_vec = sqrt(Es(ix)).*ones(1,N_pil); % generate pilot symbols and scale them by the available energy per symbol
        Txsamples = [pil_vec Loaded_syms]; % prepend pilot symbols to transmitted symbol vector

        Rxsamples = h.*Txsamples + Noise; % multiply transmitted symbols by fading channel coefficient, and add AWGN noise

        h_est = (conj(pil_vec)*(Rxsamples(1:N_pil).'))/(sum(abs(pil_vec).^2)); % estimate the channel coefficient using the pilot symbols

%         Eq_syms = Rxsamples(N_pil+1:M).*(h'/(abs(h)^2)); % equalize the received symbols using perfect channel estimate
        Eq_syms = Rxsamples(N_pil+1:M).*(h_est'/(abs(h_est)^2)); % equalize the received symbols using IMPERFECT channel estimate

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

semilogy(SNRdB,BER,'^--');
grid on;
hold on;
xlabel('SNR(dB)')
ylabel('BER')
axis tight;
title('BER vs SNR (dB) U21EC080');
legend('Location', 'best')