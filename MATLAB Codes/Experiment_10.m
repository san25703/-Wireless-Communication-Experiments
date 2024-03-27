clc
clear all
% close all

No = 1; % Noise power in linear scale
SNRdB = -20:5:20; % SNR in dB
SNR = 10.^(SNRdB./10); % SNR in linear scale
Es = SNR*No; % Symbol power is SNR*Noise Power
BER_zf = zeros(size(SNRdB)); % Initialize BER vector for ZF MIMO receiver
BER_mmse = zeros(size(SNRdB)); % Initialize BER vector for MMSE MIMO receiver
Nr = 4; % No. of Rx Antennas
Nt = 4; % No. of Tx Antennas
blklen = 1000; % Block length, No. of symbols sent in one block

ITER = 5000; % No. of Iterations

for ite = 1:ITER % iteration loop
    ite
    TxBits = randi([0,1],[Nt,blklen]); % Generate transmit bits for each Tx antenna and for each symbol in a block
    H = 1/sqrt(2).*(randn(Nr,Nt) + 1j*randn(Nr,Nt)); % Generate a MIMO channel matrix of size Nr X Nt, each coefficient is I.I.D. rayleigh fading.
    for ix = 1:length(SNRdB) % SNR loop
        Noise = sqrt(No/2).*(randn(Nr,blklen)) + 1j*randn(Nr,blklen); % Generate random noise matrix for each Rx antenna and for each symbol in a block
        TxSyms = TxBits.*2 - 1; % Convert Transmit bits to BPSK symbols directly (No use of pskmod() or qammod() functions)
        x = (sqrt(Es(ix))*1/sqrt(Nt)).*TxSyms; % Scale BPSK symbols by available power, but split that power across multiple antennas (1/sqrt(Nt) factor)
        y = H*x + Noise; % Received symbols on each Rx antenna are given by the MIMO system model y = H*x + noise
        r_zf = pinv(H)*y; % Zero forcing receiver, uses pseudo-inverse of channel matrix H to recover transmitted symbols
        r_mmse = (H'*H + (1/Es(ix)).*eye(Nt))\H'*y; % MMSE receiver, incorporates SNR information into the equation to improve BER performance
        EqBits_zf = (real(r_zf) >= 0); % decode the information bits from output of ZF receiver directly (without using pskdemod() or qamdemod())
        EqBits_mmse = (real(r_mmse) >= 0); % decode the information bits from output of MMSE receiver directly (without using pskdemod() or qamdemod())
        BER_zf(ix) = BER_zf(ix) + sum(sum(EqBits_zf ~= TxBits)); % compare bits received by ZF with original TxBits matrix and count all the errors
        BER_mmse(ix) = BER_mmse(ix) + sum(sum(EqBits_mmse ~= TxBits)); % compare bits received by MMSE with original TxBits matrix and count all the errors
    end
end
BER_zf = BER_zf./(ITER*Nt*blklen); % divide sum of all errors per SNR by total no. of bits transmitted across all iterations to get BER of ZF
BER_mmse = BER_mmse./(ITER*Nt*blklen); % divide sum of all errors per SNR by total no. of bits transmitted across all iterations to get BER of MMSE

semilogy(SNRdB,BER_zf,'b-s','linewidth',1.0,'MarkerFaceColor','b','MarkerSize',3.0);
grid on;
hold on;
axis tight;
title('MIMO BER v/s SNR(dB)');
xlabel('SNR(dB)');
ylabel('BER');
semilogy(SNRdB,BER_mmse,'r-s','linewidth',1.0,'MarkerFaceColor','r','MarkerSize',3.0);