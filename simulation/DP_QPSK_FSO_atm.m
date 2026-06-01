%% 0. Environment Setup
clc; clear; close all;
rng(2025);                  % fixed seed

%% 1. System Parameters

% 1.1 ESTOL 100G DP-QPSK baseline modulation config -----------------------
% DP-QPSK:
% QPSK = 2 bits/symbol/polarization
% DP = 2 polarizations
mod.Rs = 31.5e9;                % symbol rate [bauds]
mod.num_sym = 1e6;              % symbols per polarization
mod.M = 4;                      % QPSK modulation order
mod.Npol = 2;                   % dual polarization
% RRC Filter
mod.sps = 4;                    % samples per symbol
mod.rolloff = 0.2;              % RRC roll-off factor - OpenROAD MSA Spec
mod.span = 16;                  % RRC mod.span in symbols

% Calculated variables
Ts = 1/mod.Rs;
k = log2(mod.M);                        % bits per symbol
Rb_raw = mod.Rs * k * mod.Npol;         % raw bit rate before FEC/overheads [bit/s]
Fs = mod.Rs * mod.sps;                  % sampling frequency [Hz]
BW_null = mod.Rs * (1 + mod.rolloff);   % theoretical null-to-null bandwidth [Hz]

% 1.2 FSO scenario inputs -------------------------------------------------
% Worst-case design point: el = 15 deg, 1550 nm band, bad conditions.
% Check Link Budget for more info
scen.el_deg = 15;               % design elevation [deg]
scen.wl = 1554.13e-9;           % wavelength [m], ESTOL L2
scen.t_z= 0.891;                % zenith transmittance bad, 1550 nm
% Future additions (do not delete; comment out until used):
% scen.Cn2_model    = 'HV57';
% scen.sigma_jit_rad = ...;
% scen.theta_tx_rad  = 380e-6;
% scen.h_orbit       = 530e3;
% scen.h_ogs         = 578;

% 1.3 Visualization -------------------------------------------------------
fprintf('\n=== DP-QPSK TX Parameters ===\n');
fprintf('Symbol rate per pol      = %.2f Gbaud\n', mod.Rs/1e9);
fprintf('Raw bit rate             = %.2f Gbps\n', Rb_raw/1e9);
fprintf('Samples per symbol       = %d\n', mod.sps);
fprintf('Sampling frequency       = %.2f GSa/s\n', Fs/1e9);
fprintf('RRC roll-off             = %.2f\n', mod.rolloff);
fprintf('Expected null-null BW    = %.2f GHz\n', BW_null/1e9);
fprintf('=============================\n');

fprintf('\n=== FSO Scenario ===\n');
fprintf('Design elevation         = %.1f deg\n', scen.el_deg);
fprintf('Wavelength               = %.2f nm\n',  scen.wl*1e9);
fprintf('Zenith transmittance     = %.3f\n',     scen.t_z);
fprintf('====================\n\n');

%% 2. TX - Data Generation (PRBS)

% Not truly PRBS, in 2020_Nazir_32GBaud_DP-QPSK, they use 2^17 - 1 PRBS
% Generate independent random binary data for X and Y polarizations
% randi([0 1], rows, columns) creates a column vector of 1s and 0s
bits_X = randi([0 1], mod.num_sym * k, 1);
bits_Y = randi([0 1], mod.num_sym * k, 1);

%% 3. Symbol Mapping (DP-QPSK)

% Convert the the 2-bit pairs (e.g., [1 0]) into integers (0 to 3). 
% MSB by default
ints_X = bit2int(bits_X, k);
ints_Y = bit2int(bits_Y, k);

% Modulate the symbols into a gray QPSK
symbols_X = pskmod(ints_X, mod.M, pi/4, 'gray');
symbols_Y = pskmod(ints_Y, mod.M, pi/4, 'gray');
% Based on the OpenROAD standard the modulation are different
    %bX = reshape(bits_X, k, []).';     % rows: [bI bQ]
    %symbols_X = ((2*bX(:,1)-1) + 1j*(2*bX(:,2)-1)) / sqrt(2);
%bY = reshape(bits_Y, k, []).';
%symbols_Y = ((2*bY(:,1)-1) + 1j*(2*bY(:,2)-1)) / sqrt(2);

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
rrc = rcosdesign(mod.rolloff, mod.span, mod.sps, 'sqrt');
%impz(rrc)

% Filter length
% N_h = mod.span * mod.sps + 1;
fprintf('RRC filter length        = %d taps\n', length(rrc));
fprintf('RRC filter energy        = %.6f\n\n', sum(abs(rrc).^2));

% Upsample and filter the data for pulse shaping
% https://www.mathworks.com/help/signal/ref/upfirdn.html
% mod.sps - upsampling
tx_X = upfirdn(symbols_X, rrc, mod.sps);
tx_Y = upfirdn(symbols_Y, rrc, mod.sps);

% Average waveform sample power
Ptx_X_norm = mean(abs(tx_X).^2);
Ptx_Y_norm = mean(abs(tx_Y).^2);

fprintf('Mean sample power X-pol  = %.6f\n', Ptx_X_norm);
fprintf('Mean sample power Y-pol  = %.6f\n', Ptx_Y_norm);
fprintf('Expected approx.         = %.6f (= 1/sps)\n\n', 1/mod.sps);

%% 5. TX - Visualization

% 5.1 Ideal QPSK constellation before pulse shaping -----------------------
figure('Name','Tx Symbols');

subplot(1,2,1);
plot(real(symbols_X), imag(symbols_X), 'bo', 'MarkerFaceColor','b');
grid on; axis square;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
title('X-Polarization'); 
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

t_taps = (-(length(rrc)-1)/2 : (length(rrc)-1)/2) / mod.sps;

stem(t_taps, rrc, 'filled');
grid on;
xlabel('Time [symbols]');
ylabel('Amplitude');
title(sprintf('RRC impulse response: rolloff = %.2f, span = %d, sps = %d', ...
    mod.rolloff, mod.span, mod.sps));

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
eye_tx1 = real(tx_X(mod.span*mod.sps+1 : mod.span*mod.sps+4000));
eye_tx2 = real(tx_Y(mod.span*mod.sps+1 : mod.span*mod.sps+4000));
eye_tx = [eye_tx1, eye_tx2];

eyediagram(eye_tx, 2*mod.sps);

% Grab the axes handles from the current figure
ax = findobj(gcf, 'Type', 'axes');

% Note: findobj grabs handles in reverse order of creation!
% ax(1) is the RIGHT plot, ax(2) is the LEFT plot.
title(ax(2), 'X-pol I-component eye diagram');
title(ax(1), 'Y-pol I-component eye diagram');

%% 6. Free Space Optical Channel

% Beer-Lambert with plane-parallel airmass m = 1/sin(el).
% Valid for el >~ 10 deg; below that, use Kasten-Young airmass.
% Refer to Giggenbach paper
el_rad = deg2rad(scen.el_deg);
m_air = 1/sin(el_rad);                      % airmass
h_atm = scen.t_z ^ m_air;                   % linear power transmittance
a_atm_dB = 10*log10(h_atm);     % [dB], negative = loss

% Physically this is one optical DP-QPSK signal.
tx_DP = [tx_X, tx_Y];

% AWGN channel: receiver gets a noised signal
EbN0_dB_vec = 0:1:12;               % Sweep of values to plot the BER
EbN0_lin = 10.^(EbN0_dB_vec/10);

BER_X = zeros(size(EbN0_dB_vec));
BER_Y = zeros(size(EbN0_dB_vec));
BER_total = zeros(size(EbN0_dB_vec));

% Exact uncoded Gray-QPSK theory
BER_theory = qfunc(sqrt(2*EbN0_lin));

% Filters group delay
% Tx RRC group delay = ((mod.span*mod.sps+1) - 1) / 2
% Rx RRC group delay = ((mod.span*mod.sps+1) - 1) / 2
% Total delay  = (((mod.span*mod.sps+1) - 1) / 2) + (((mod.span*mod.sps+1) - 1) / 2)
total_delay = mod.span * mod.sps;

for i = 1:length(EbN0_dB_vec)
    EbN0 = 10^(EbN0_dB_vec(i)/10);

    % QPSK: Es/N0 = k*Eb/N0, with Es = 1
    EsN0 = k * EbN0;
    N0 = 1 / EsN0;
    sigma = sqrt(N0/2);

    % nI_X + nQ_X = sigma * randn(size(tx_DP(:,1))) + sigma * randn(size(tx_DP(:,1)))
    noise_X = sigma * (randn(size(tx_DP(:,1))) ...
        + 1j*randn(size(tx_DP(:,1))));
    noise_Y = sigma * (randn(size(tx_DP(:,2))) ...
        + 1j*randn(size(tx_DP(:,2))));

    % Apply FSO channel (amplitude scaling = sqrt of power gain)
    % h_ch = h_atm + h_tur and h_poin
    % TO DO
    h_ch = h_atm;

    rx_X = sqrt(h_ch) * tx_DP(:,1) + noise_X;
    rx_Y = sqrt(h_ch) * tx_DP(:,2) + noise_Y;

    % Receiver matched filter uses the same RRC filter
    % No up/downsampling
    rxmf_X = upfirdn(rx_X, rrc);
    rxmf_Y = upfirdn(rx_Y, rrc);

    % Downsample at symbol instants
    rx_symbols_X = rxmf_X(total_delay + 1 : mod.sps : total_delay + mod.num_sym*mod.sps);
    rx_symbols_Y = rxmf_Y(total_delay + 1 : mod.sps : total_delay + mod.num_sym*mod.sps);

    % Demodulate the gray QPSK into symbols
    ints_hat_X = pskdemod(rx_symbols_X, mod.M, pi/4, 'gray');
    ints_hat_Y = pskdemod(rx_symbols_Y, mod.M, pi/4, 'gray');

    % Convert the the integers (0 to 3) into 2-bit pairs (e.g., [1 0]). 
    % MSB by default
    bits_hat_X = int2bit(ints_hat_X, k);
    bits_hat_Y = int2bit(ints_hat_Y, k);

    % BER
    BER_X(i) = mean(bits_X ~= bits_hat_X);
    BER_Y(i) = mean(bits_Y ~= bits_hat_Y);
    BER_total(i) = mean([bits_X ~= bits_hat_X; bits_Y ~= bits_hat_Y]);

end

% 9.1 Visualization -------------------------------------------------------
fprintf('=== FSO Channel: atmosphere ===\n');
fprintf('Airmass                  = %.3f\n',      m_air);
fprintf('h_atm (linear power)     = %.4f\n',      h_atm);
fprintf('a_atm                    = %.2f dB\n',   a_atm_dB);
fprintf('===============================\n\n');

%% 10. RX - Visualization

% 10.1 Received spectrum --------------------------------------------------
figure('Name','Rx Spectrum after Matched Filter');

% https://es.mathworks.com/help/signal/ref/pwelch.html
[Pxx_r, f_psd_r] = pwelch(rxmf_X, win, noverlap, nfft, Fs, 'centered');
[Pyy_r, ~]     = pwelch(rxmf_Y, win, noverlap, nfft, Fs, 'centered');

plot(f_psd_r/1e9, 10*log10(Pxx_r/max(Pxx_r)), 'b', 'LineWidth', 1.2);
hold on;
plot(f_psd_r/1e9, 10*log10(Pyy_r/max(Pyy_r)), 'r--', 'LineWidth', 1.2);
grid on;

xlabel('Frequency [GHz]');
ylabel('Normalized PSD [dB]');
title(sprintf('DP-QPSK Rx spectrum after matched filter, expected null-null BW = %.2f GHz', ...
    BW_null/1e9));
legend('X-pol', 'Y-pol', 'Location', 'best');
ylim([-80 5]);

% 10.2 Eye diagrams after channel -----------------------------------------
eye_rx1 = real(rxmf_X(mod.span*mod.sps+1 : mod.span*mod.sps+4000));
eye_rx2 = real(rxmf_Y(mod.span*mod.sps+1 : mod.span*mod.sps+4000));
eye_rx = [eye_rx1, eye_rx2];

eyediagram(eye_rx, 2*mod.sps);

% Grab the axes handles from the current figure
ax2 = findobj(gcf, 'Type', 'axes');

title(ax2(2), 'X-pol I-component eye diagram after matched filter');
title(ax2(1), 'Y-pol I-component eye diagram after matched filter');

% 10.2 Received constellation after matched filter ------------------------
figure('Name','Rx Symbols after Matched Filter');

subplot(1,2,1);
plot(real(rx_symbols_X(1:2000)), imag(rx_symbols_X(1:2000)), 'bo', ...
    'MarkerFaceColor','b');
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

% 10.4 BER table ----------------------------------------------------------

% OpenROADM / ESTOL oFEC pre-FEC BER threshold
BER_preFEC = 2.0e-2;

% Theoretical Eb/N0 at which Gray-QPSK reaches the pre-FEC threshold:
EbN0_thr_lin = erfcinv(2*BER_preFEC)^2;
EbN0_thr_dB  = 10*log10(EbN0_thr_lin);

% Ber for applied channel
BER_theory_ch = qfunc(sqrt(2*EbN0_lin*h_ch));

% Print BER table
BER_table = table( ...
    EbN0_dB_vec(:), ...
    BER_theory(:), ...
    BER_theory_ch(:), ...
    BER_X(:), ...
    BER_Y(:), ...
    BER_total(:), ...
    'VariableNames', {'EbN0_dB','BER_theory','BER_theory_ch','BER_X',...
    'BER_Y','BER_total'} );

disp(' ');
disp('=== BER vs Eb/N0 Table ===');
disp(BER_table);

fprintf('\nPre-FEC BER threshold = %.2e\n', BER_preFEC);
fprintf('Theoretical Eb/N0 at pre-FEC threshold = %.3f dB\n\n', EbN0_thr_dB);

% 10.4 BER plot -----------------------------------------------------------

figure('Name','BER vs EbN0_lin');

semilogy(EbN0_dB_vec, BER_theory, 'k-', 'LineWidth', 1.5);
hold on;
semilogy(EbN0_dB_vec, BER_theory_ch, 'g--', 'LineWidth', 1.5);
semilogy(EbN0_dB_vec, BER_X, 'bo-');
semilogy(EbN0_dB_vec, BER_Y, 'Rs-');
semilogy(EbN0_dB_vec, BER_total, 'md-');
% Horizontal pre-FEC threshold
yline(BER_preFEC, 'k--', 'Pre-FEC BER = 2.0\times10^{-2}', ...
    'LineWidth', 1.2, ...
    'LabelHorizontalAlignment','left', ...
    'LabelVerticalAlignment','bottom');
% Vertical theoretical threshold crossing
xline(EbN0_thr_dB, 'k:', sprintf('%.2f dB', EbN0_thr_dB), ...
    'LineWidth', 1.2, ...
    'LabelOrientation','horizontal', ...
    'LabelHorizontalAlignment','left', ...
    'LabelVerticalAlignment','bottom');
grid on;
grid minor;
xlabel('Tx. E_b/N_0 before atmospheric attenuation [dB]');
ylabel('BER');
title('Gray-coded DP-QPSK with deterministic atmospheric attenuation + AWGN');
legend('Theory AWGN QPSK', 'Theory Channel QPSK', 'X-Pol.', 'Y-Pol.', ...
    'Total', 'Pre-FEC threshold', 'Threshold crossing');
ylim([1e-6 1]);
xlim([min(EbN0_dB_vec) max(EbN0_dB_vec)]);
