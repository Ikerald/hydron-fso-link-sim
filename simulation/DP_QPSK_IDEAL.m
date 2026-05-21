%% 0. Environment Setup
clc; clear; close all;
rng(2025);                  % fixed seed

%% 1. System Parameters

% ESTOL 100G DP-QPSK baseline
Rs = 31.5e9;                % symbol rate

Ts = 1/Rs;

% Modulation --------------------------------------------------------------
% DP-QPSK:
% QPSK = 2 bits/symbol/polarization
% DP = 2 polarizations
num_symbols = 1e5;          % symbols per polarization ??
M = 4;                      % QPSK modulation order
k = log2(M);                % bits per symbol
Npol = 2;                   % dual polarization

Rb_raw = Rs * k * Npol;     % raw bit rate before FEC/overheads [bit/s]

% RRC Filter --------------------------------------------------------------
sps = 4;                        % samples per symbol
rolloff = 0.2;                  % RRC roll-off factor - OpenROAD MSA Spec
span = 16;                      % RRC span in symbols

Fs = Rs * sps;                  % sampling frequency [Hz]
BW_null = Rs * (1 + rolloff);   % theoretical null-to-null bandwidth [Hz]

% Visualization -----------------------------------------------------------
fprintf('\n=== DP-QPSK TX Parameters ===\n');
fprintf('Symbol rate per pol      = %.2f Gbaud\n', Rs/1e9);
fprintf('Raw bit rate             = %.2f Gbps\n', Rb_raw/1e9);
fprintf('Samples per symbol       = %d\n', sps);
fprintf('Sampling frequency       = %.2f GSa/s\n', Fs/1e9);
fprintf('RRC roll-off             = %.2f\n', rolloff);
fprintf('Expected null-null BW    = %.2f GHz\n', BW_null/1e9);
fprintf('=============================\n\n');

%% 2. TX - Data Generation (PRBS)

% Generate independent random binary data for X and Y polarizations
% randi([0 1], rows, columns) creates a column vector of 1s and 0s
bits_X = randi([0 1], num_symbols * k, 1);
bits_Y = randi([0 1], num_symbols * k, 1);

%% 3. Symbol Mapping (DP-QPSK)

% Convert the the 2-bit pairs (e.g., [1 0]) into integers (0 to 3). 
% MSB by default
ints_X = bit2int(bits_X, k);
ints_Y = bit2int(bits_Y, k);

% Modulate the symbols into a gray QPSK
symbols_X = pskmod(ints_X, M, pi/4, 'gray');
symbols_Y = pskmod(ints_Y, M, pi/4, 'gray');

% Check symbol energy
Es_X = mean(abs(symbols_X).^2);
Es_Y = mean(abs(symbols_Y).^2);

fprintf('Mean symbol energy X-pol = %.6f\n', Es_X);
fprintf('Mean symbol energy Y-pol = %.6f\n\n', Es_Y);

assert(abs(Es_X - 1) < 1e-6, 'X-pol symbol energy is not 1');
assert(abs(Es_Y - 1) < 1e-6, 'Y-pol symbol energy is not 1');

%% 4. TX - Pulse Shaping (Root-Raised Cosine)

% We create the Raised Cosine Filter
% https://es.mathworks.com/help/signal/ref/rcosdesign.html
rrc = rcosdesign(rolloff, span, sps, 'sqrt');
%impz(rrc)

% Filter length
% N_h = span * sps + 1;
fprintf('RRC filter length        = %d taps\n', length(rrc));
fprintf('RRC filter energy        = %.6f\n\n', sum(abs(rrc).^2));

% Upsample and filter the data for pulse shaping
% https://www.mathworks.com/help/signal/ref/upfirdn.html
% sps - upsampling
tx_X = upfirdn(symbols_X, rrc, sps);
tx_Y = upfirdn(symbols_Y, rrc, sps);

% Average waveform sample power
Ptx_X_norm = mean(abs(tx_X).^2);
Ptx_Y_norm = mean(abs(tx_Y).^2);

fprintf('Mean sample power X-pol  = %.6f\n', Ptx_X_norm);
fprintf('Mean sample power Y-pol  = %.6f\n', Ptx_Y_norm);
fprintf('Expected approx.         = %.6f (= 1/sps)\n\n', 1/sps);

%% 5. TX - Visualization

% 5.1 Ideal QPSK constellation before pulse shaping -----------------------
figure('Name','Tx symbols');

subplot(1,2,1);
plot(real(symbols_X), imag(symbols_X), 'bo', 'MarkerFaceColor','b');
grid on; axis square;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
title('X-Polarization (Ideal)'); 
xlim([-1.5 1.5]); ylim([-1.5 1.5]);

subplot(1,2,2);
plot(real(symbols_Y), imag(symbols_Y), 'ro', 'MarkerFaceColor','r');
grid on; axis square;
xlabel('In-Phase (I)'); 
ylabel('Quadrature (Q)');
title('Y-Polarization'); 
xlim([-1.5 1.5]); ylim([-1.5 1.5]);

% 5.2 RRC impulse response ----------------------------------------------
figure('Name','RRC Impulse Response');

t_taps = (-(length(rrc)-1)/2 : (length(rrc)-1)/2) / sps;

stem(t_taps, rrc, 'filled');
grid on;
xlabel('Time [symbols]');
ylabel('Amplitude');
title(sprintf('RRC impulse response: rolloff = %.2f, span = %d, sps = %d', ...
    rolloff, span, sps));

% 5.3 RRC frequency response ----------------------------------------------
figure('Name','RRC Frequency Response');

% H = freq response
% f = frequencies in Hz
[H, f] = freqz(rrc, 1, 4096, Fs);

plot(f/1e9, 20*log10(abs(H)/max(abs(H))), 'LineWidth', 1.2);
grid on;
xlabel('Frequency [GHz]');
ylabel('Magnitude [dB]');
title('RRC frequency response');
ylim([-80 5]);

% 5.4 Transmitted spectrum ------------------------------------------------
figure('Name','Tx Spectrum');

nfft = 4096;
win = hann(2048);
noverlap = 1024;

% https://es.mathworks.com/help/signal/ref/pwelch.html
[Pxx, f_psd] = pwelch(tx_X, win, noverlap, nfft, Fs, 'centered');
[Pyy, ~]     = pwelch(tx_Y, win, noverlap, nfft, Fs, 'centered');

plot(f_psd/1e9, 10*log10(Pxx/max(Pxx)), 'b', 'LineWidth', 1.2);
hold on;
plot(f_psd/1e9, 10*log10(Pyy/max(Pyy)), 'r--', 'LineWidth', 1.2);
grid on;

xlabel('Frequency [GHz]');
ylabel('Normalized PSD [dB]');
title(sprintf('DP-QPSK Tx spectrum, expected null-null BW = %.2f GHz', ...
    BW_null/1e9));
legend('X-pol', 'Y-pol', 'Location', 'best');
ylim([-80 5]);

% 5.5 Eye diagrams of pulse-shaped waveform -------------------------------
eyediagram(real(tx_X(span*sps+1 : span*sps+4000)), 2*sps);
title('X-pol I-component tx eye diagram');

eyediagram(real(tx_Y(span*sps+1 : span*sps+4000)), 2*sps);
title('Y-pol I-component tx eye diagram');

%% 6. Ideal Channel: no noise, no fading, no impairments

% Physically this is one optical DP-QPSK signal.
% In baseband simulation we keep the two polarizations separate.
tx_DP = [tx_X, tx_Y];

% Ideal channel: receiver gets exactly what transmitter sends
rx_X = tx_DP(:,1);
rx_Y = tx_DP(:,2);

%% 7. RX - Matched Filter

% Receiver matched filter uses the same RRC filter
% No up/downsampling
rxmf_X = upfirdn(rx_X, rrc);
rxmf_Y = upfirdn(rx_Y, rrc);

% Total delay:

% Tx RRC group delay = ((span*sps+1) - 1) / 2
% Rx RRC group delay = ((span*sps+1) - 1) / 2
% Total delay  = (((span*sps+1) - 1) / 2) + (((span*sps+1) - 1) / 2)
total_delay = span * sps;

% Downsample at symbol instants
rx_symbols_X = rxmf_X(total_delay + 1 : sps : total_delay + num_symbols*sps);
rx_symbols_Y = rxmf_Y(total_delay + 1 : sps : total_delay + num_symbols*sps);

%% 8. RX - QPSK Demodulation

% Demodulate the gray QPSK into symbols
ints_hat_X = pskdemod(rx_symbols_X, M, pi/4, 'gray');
ints_hat_Y = pskdemod(rx_symbols_Y, M, pi/4, 'gray');

% Convert the the integers (0 to 3) into 2-bit pairs (e.g., [1 0]). 
% MSB by default
bits_hat_X = int2bit(ints_hat_X, k);
bits_hat_Y = int2bit(ints_hat_Y, k);

%% 9. BER Calculation

BER_X = mean(bits_X ~= bits_hat_X);
BER_Y = mean(bits_Y ~= bits_hat_Y);

BER_total = mean([bits_X ~= bits_hat_X; bits_Y ~= bits_hat_Y]);

fprintf('\n=== Ideal Receiver Check ===\n');
fprintf('BER X-pol   = %.3e\n', BER_X);
fprintf('BER Y-pol   = %.3e\n', BER_Y);
fprintf('BER total   = %.3e\n', BER_total);
fprintf('============================\n\n');

%% 10. RX - Visualization

% 10.1 Eye diagrams after matched filter ----------------------------------
eyediagram(real(rx_X(span*sps+1 : span*sps+4000)), 2*sps);
title('X-pol I-component rx eye diagram after filter');

eyediagram(real(rx_Y(span*sps+1 : span*sps+4000)), 2*sps);
title('Y-pol I-component rx eye diagram after filter');

% 10.2 Received constellation after matched filter ------------------------
figure('Name','Received symbols after ideal matched filter');

subplot(1,2,1);
plot(real(rx_symbols_X(1:2000)), imag(rx_symbols_X(1:2000)), 'bo', ...
    'MarkerFaceColor','b');
grid on; axis square;
grid on; axis square;
xlabel('In-Phase (I)'); 
ylabel('Quadrature (Q)');
title('Recovered X-Pol. symbols'); 
xlim([-1.5 1.5]); ylim([-1.5 1.5]);

subplot(1,2,2);
plot(real(rx_symbols_Y(1:2000)), imag(rx_symbols_Y(1:2000)), 'ro', ...
    'MarkerFaceColor','r');
grid on; axis square;
xlabel('In-Phase (I)'); 
ylabel('Quadrature (Q)');
title('Recovered Y-Pol. symbols'); 
xlim([-1.5 1.5]); ylim([-1.5 1.5]);