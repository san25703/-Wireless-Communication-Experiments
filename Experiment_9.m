EsdB = -10:2:10; % Symbol energy vector in dB
Es = 10.^(EsdB/10); % Symbol energy in linear
No = 1; % Noise power in linear
SNR = Es/No; % SNR in linear
SNRdB = 10*log10(SNR); % SNR in dB
BER = zeros(length(Es),1); % initialize BER vector to store results
M=64; % No. of time-domain symbols in one transmitted block
ITER = 5000; % No. of iterations/channel realizations
Nr = 4; % No. of Receiving Antennas, recommended values [1, 2, 3, 4]

K_rice = 1; % K factor, the ratio of the power between the LOS path and all the NLOS paths, K = v^2/sigma_h^2
sigma_h = 1;
v_factor = sqrt((sigma_h^2)*K_rice); % get the amplitude of the LOS path

mod_ord = 16; % Modulation Order

for ite = 1:ITER % loop over all iterations
    ite

    TxBits = randi([0,1],[1, M*log2(mod_ord)]); % generate bits equal to M*log2(modulation order), to convert to M symbols
    shaped_bits = reshape(TxBits,log2(mod_ord),M); % reshape bit vector into matrix of size (log2(mod_ord) x M)
    transposed_bits = shaped_bits.'; % transpose the matrix
    decimal_syms = bi2de(transposed_bits); % apply binary to decimal conversion
    decimal_syms = decimal_syms.'; % transpose the decimal result
    Mod_syms = pskmod(decimal_syms,mod_ord); % perform psk modulation

    h_vec = zeros(1,Nr); % initialize channel vector containing fading coefficients for each receiving antenna link
    for rx = 1:Nr
        mag_h = abs((sigma_h*sqrt(1/2)*(randn+1j*randn)) + v_factor); % generate the magnitude of a rician fading channel coefficient (Prob. Distribution ~Rice(v,sigma_h))
        theta_h = rand*2*pi; % generate the phase of the rician fading channel coefficient (Prob. Distribution ~U(0,2*pi))
        h_vec(rx) = mag_h*(cos(theta_h) + 1j*sin(theta_h)); % generate a single-tap rician fading channel coefficient and store into channel vector
    end

    Noise = sqrt(No/2)*(randn(Nr,M)+1j*randn(Nr,M)); % generate AWGN noise, different set of values for each antenna (Nr,M)

    for ix = 1:length(Es) % loop over all SNR
        Loaded_syms = sqrt(Es(ix)).*Mod_syms; % scale info symbols by the available energy per symbol
        Txsamples = Loaded_syms;

        Rxsamples = zeros(Nr,M);
        for kx = 1:Nr
            Rxsamples(kx,:) = h_vec(kx).*Txsamples + Noise(kx,:); % multiply transmitted symbols by a different fading channel coefficient for each antenna link, and add AWGN noise
        end

        %%%% SELECTION COMBINING (SC)
        [max_h_gain, index_max] = max(abs(h_vec).^2); % find the coefficient and the corresponding antenna with maximum gain
        Eq_syms = (conj(h_vec(index_max))/(abs(h_vec(index_max))^2)).*Rxsamples(index_max,:); % select only the received vector/antenna whose channel gain is maximum, and equalize the received symbols on that antenna.

        %%%% MAXIMAL RATIO COMBINING (MRC)
%         Eq_syms = (conj(h_vec)./(norm(h_vec)^2))*Rxsamples; % equalize the received symbols on each antenna per symbol slot using MRC

        Rx_syms = pskdemod(Eq_syms,mod_ord); % perform psk demodulation to convert to integers
        Rx_syms = Rx_syms.'; % transpose the integer result
        Rx_decimal_syms = de2bi(Rx_syms,log2(mod_ord)); % apply decimal to binary conversion with number of bits equal to log2(mod_ord)
        Rx_transposed_bits = Rx_decimal_syms.'; % transpose the resultant bit matrix
        DecodedBits = reshape(Rx_transposed_bits,1, M*log2(mod_ord)); % reshape the bit matrix to a bit vector

        BER(ix) = BER(ix) + sum(DecodedBits ~= TxBits); % count all bit errors and sum them
    end
end
BER=BER/(M*log2(mod_ord)*ITER); % divide BER vector by total number of bits transmitted across all iterations

semilogy(SNRdB,BER,'--^');
grid on;
hold on;
xlabel('SNR(dB)')
ylabel('BER')
axis tight;
title('SIMO BER vs SNR (dB) U21EC080');
legend('Location', 'best')