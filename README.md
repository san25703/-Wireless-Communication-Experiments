# Wireless Communication Experiments

This repository contains MATLAB implementations for various wireless communication experiments. Each experiment focuses on different modulation schemes, channel conditions, and techniques used in wireless communication systems. The experiments are aimed at analyzing the performance of these systems under different scenarios.

## List of Experiments

1. **Various Modulation Schemes Simulation**: Implement and simulate different modulation schemes and analyze their performance.
2. **M-PSK Modulation with OFDM on AWGN Channel**: Implement and simulate M-PSK modulation technique with OFDM on an AWGN channel for M=2,4,16. Plot BER vs SNR.
3. **M-PSK Modulation with OFDM on Single-Tap Rayleigh Fading Channel**: Implement and simulate M-PSK modulation technique with OFDM on a single-tap Rayleigh fading channel with AWGN for M=2,4,16. Plot BER vs SNR.
4. **M-PSK Modulation with OFDM on Single-Tap Rician Fading Channel**: Implement and simulate M-PSK modulation technique with OFDM on a single-tap Rician fading channel with AWGN for M=2,4,16. Plot BER vs SNR for different values of the Rician K factor.
5. **Multi-Tap Rayleigh Fading Channel with ISI**: Simulate a multi-tap Rayleigh fading channel with inter-symbol interference (ISI). Perform time-domain block-wise zero-forcing equalization to eliminate ISI at the receiver and plot BER vs SNR for an M-QAM system.
6. **OFDM Communication System with Multi-Tap Fading Channel**: Implement and simulate an Orthogonal Frequency Division Multiplexing (OFDM) communication system with a multi-tap fading channel. Perform frequency domain equalization on each subcarrier at the receiver and plot BER vs SNR for M-QAM symbols.
7. **Single-Tap Rayleigh Fading Channel with Imperfect CSI**: Simulate a single-tap Rayleigh fading channel with imperfect channel state information (CSI) for an M-QAM communication system for M=4. Perform channel estimation at the receiver using pilot symbols and decode the information symbols using the channel estimates. Plot BER vs SNR for various numbers of pilot symbols (per block of transmitted symbols).
8. **OFDM System with Frequency Diversity on Multi-Tap Rician Fading Channel**: Implement an OFDM system with frequency diversity on a multi-tap Rician fading channel. Transmit the same M-QAM symbols on multiple subcarriers and perform maximal ratio combining (MRC) at the received frequency-domain symbols. Plot BER vs SNR and observe the effect of increasing the number of repeated symbols.
9. **Multi-Antenna Spatial Diversity Techniques on Single Tap Rician Fading Channel**: Implement and simulate multi-antenna Spatial Diversity techniques using a 2x1, 3x1, and 4x1 Single Input Multiple Output (SIMO) communication system on a single tap Rician fading channel. Plot BER vs SNR for Selection Combining and Maximal Ratio Combining (MRC) techniques.
10. **2x2 MIMO Communication System on Single Tap Rayleigh Fading Channel**: Implement and simulate a 2x2 Multiple Input Multiple Output (MIMO) communication system on a single tap Rayleigh fading channel. Plot BER vs SNR for Zero-Forcing (ZF) and Linear Minimum Mean Square Error (LMMSE) receivers.

## Common Terms used in this Repository

- **Modulation Schemes**: Techniques used to encode digital data onto an analog carrier signal for transmission. Examples include Amplitude Shift Keying (ASK), Frequency Shift Keying (FSK), Phase Shift Keying (PSK), and Quadrature Amplitude Modulation (QAM).
- **Orthogonal Frequency Division Multiplexing (OFDM)**: A multi-carrier modulation technique used in wireless communication systems. It divides the available spectrum into multiple subcarriers, which are modulated using QAM or PSK.
- **BER vs SNR**: Bit Error Rate (BER) vs Signal-to-Noise Ratio (SNR) is a plot showing the relationship between the bit error rate of a communication system and the signal-to-noise ratio of the received signal.
- **Rayleigh Fading Channel**: A fading channel model in which the magnitude of the channel response follows a Rayleigh distribution. It is often used to model wireless communication channels with no direct line of sight between transmitter and receiver.
- **Rician Fading Channel**: A fading channel model that combines both direct path (line of sight) and scattered path components. It is characterized by a Rician K factor, which represents the ratio of the power in the line-of-sight component to the power in the scattered component.
- **Inter-Symbol Interference (ISI)**: A distortion in a communication system caused by symbols interfering with each other due to channel characteristics.
- **Channel State Information (CSI)**: Information about the characteristics of the communication channel, such as fading, noise, and interference, which is used to optimize transmission and reception.
- **Channel Estimation**: The process of estimating the channel response at the receiver based on known training symbols or pilot signals.
- **Maximal Ratio Combining (MRC)**: A diversity combining technique used at the receiver to improve the reliability of communication in fading channels. It combines the signals received from multiple antennas or subcarriers with different weights to maximize the received signal power.
- **Frequency Diversity**: A technique used to combat fading in wireless communication by transmitting the same information over multiple frequency-selective paths.
- **Selection Combining**: A diversity combining technique where the receiver selects the signal with the highest instantaneous signal-to-noise ratio (SNR) among the available antennas.
- **SIMO (Single Input Multiple Output)**: A communication system with a single transmitter and multiple receivers.
- **MIMO (Multiple Input Multiple Output)**: A communication system with multiple antennas at both the transmitter and receiver.
- **Zero-Forcing (ZF)**: A receiver technique used in MIMO systems to nullify the interference caused by multiple transmit antennas.
- **Linear Minimum Mean Square Error (LMMSE)**: A receiver technique used in MIMO systems to minimize the mean square error between the transmitted and received signals while considering channel characteristics.
