%% 0. Environment Setup
clc; clear; close all;
rng(2025);                  % fixed seed

% Latex setup
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
% Typography
set(groot, 'defaultAxesFontSize', 12);                 % tick labels
set(groot, 'defaultAxesLabelFontSizeMultiplier', 13/12); % axis labels ≈ 13 pt
set(groot, 'defaultAxesTitleFontSizeMultiplier', 14/12); % title ≈ 14 pt
set(groot, 'defaultAxesTitleFontWeight', 'normal');
% Lines and axes
set(groot, 'defaultLineLineWidth', 1.3);
set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');


%% 1. System Parameters

% 1.1 ESTOL 100G DP-QPSK baseline modulation config -----------------------
% DP-QPSK:
% QPSK = 2 bits/symbol/polarization
% DP = 2 polarizations
sig.Rs = 31.5e9;                % symbol rate [bauds]
sig.num_sym = 1e6;              % symbols per polarization
sig.M = 4;                      % QPSK modulation order
sig.Npol = 2;                   % dual polarization
% RRC Filter
sig.sps = 4;                    % samples per symbol
sig.rolloff = 0.2;              % RRC roll-off factor - OpenROAD MSA Spec
sig.span = 16;                  % RRC sig.span in symbols
% Calculated variables
sig.Ts = 1/sig.Rs;
sig.k = log2(sig.M);                        % bits per symbol
sig.Rb_raw = sig.Rs * sig.k * sig.Npol;         % raw bit rate before FEC/overheads [bit/s]
sig.Fs = sig.Rs * sig.sps;                  % sampling frequency [Hz]
sig.BW_null = sig.Rs * (1 + sig.rolloff);   % theoretical null-to-null bandwidth [Hz]

% 1.2 FSO scenario inputs -------------------------------------------------
% Worst-case design point: el = 15 deg, 1550 nm band, bad conditions.
% Check Link Budget for more info
% Atmospheric effects
scen.el_deg = 15;           % design elevation [deg]
scen.el_rad = deg2rad(scen.el_deg);
scen.wl = 1554.13e-9;       % wavelength [m], ESTOL L2
scen.t_z= 0.891;            % zenith transmittance bad, 1550 nm
% Turbulence effects
scen.h_orbit = 530e3;       % satellite altitude [m]
scen.h_site = 578;          % site altitude above MSL [m]
scen.h_site_ogs = 598;      % OGS altitude above MSL [m]
scen.v_wind = 21;           % HV-5/7 RMS high-altitude wind [m/s]
scen.A0 = 1.7e-14;          % ground-level Cn^2 [m^-2/3], HV-5/7 default
% Calculated variables
scen.h_ogs = scen.h_site_ogs - scen.h_site;

% Future additions (do not delete; comment out until used):
% scen.sigma_jit_rad = ...;
% scen.theta_tx_rad  = 380e-6;


% 1.3 FSO channel inputs --------------------------------------------------
% Block Fading
% Physical fading block = T_c * R_s
%ch.nsym_fade = 5e3;                        % Manually increased to notice effect
ch.nsym_fade = round(0.2e-3 * sig.Rs);      % symbols per turbulence block


% 1.4 Visualization -------------------------------------------------------
fprintf('\n=== DP-QPSK TX Parameters ===\n');
fprintf('Symbol rate per pol      = %.2f Gbaud\n', sig.Rs/1e9);
fprintf('Raw bit rate             = %.2f Gbps\n', sig.Rb_raw/1e9);
fprintf('Samples per symbol       = %d\n', sig.sps);
fprintf('Sampling frequency       = %.2f GSa/s\n', sig.Fs/1e9);
fprintf('RRC roll-off             = %.2f\n', sig.rolloff);
fprintf('Expected null-null BW    = %.2f GHz\n', sig.BW_null/1e9);
fprintf('=============================\n');

fprintf('\n=== FSO Scenario ===\n');
fprintf('Site altitude            = %.1f m\n',   scen.h_site);
fprintf('OGS height               = %.1f m\n',   scen.h_ogs);
fprintf('Site OGS altitude        = %.1f m\n',   scen.h_site_ogs);
fprintf('Design elevation         = %.1f deg\n', scen.el_deg);
fprintf('Wavelength               = %.2f nm\n',  scen.wl*1e9);
fprintf('Zenith transmittance     = %.3f\n',     scen.t_z);
fprintf('====================\n\n');


%% 2. TX - Data Generation (PRBS)

% Not truly PRBS, in 2020_Nazir_32GBaud_DP-QPSK, they use 2^17 - 1 PRBS
% Generate independent random binary data for X and Y polarizations
% randi([0 1], rows, columns) creates a column vector of 1s and 0s
bits_X = randi([0 1], sig.num_sym * sig.k, 1);
bits_Y = randi([0 1], sig.num_sym * sig.k, 1);


%% 3. TX - Symbol Mapping (DP-QPSK)

% Convert the the 2-bit pairs (e.g., [1 0]) into integers (0 to 3). 
% MSB by default
ints_X = bit2int(bits_X, sig.k);
ints_Y = bit2int(bits_Y, sig.k);

% Modulate the symbols into a gray QPSK
symbols_X = pskmod(ints_X, sig.M, pi/4, 'gray');
symbols_Y = pskmod(ints_Y, sig.M, pi/4, 'gray');
% Based on the OpenROAD standard the modulation are different
    %bX = reshape(bits_X, sig.k, []).';     % rows: [bI bQ]
    %symbols_X = ((2*bX(:,1)-1) + 1j*(2*bX(:,2)-1)) / sqrt(2);
%bY = reshape(bits_Y, sig.k, []).';
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
rrc = rcosdesign(sig.rolloff, sig.span, sig.sps, 'sqrt');
%impz(rrc)

% Filter length
% N_h = sig.span * sig.sps + 1;
fprintf('RRC filter length        = %d taps\n', length(rrc));
fprintf('RRC filter energy        = %.6f\n\n', sum(abs(rrc).^2));

% Upsample and filter the data for pulse shaping
% https://www.mathworks.com/help/signal/ref/upfirdn.html
% sig.sps - upsampling
tx_X = upfirdn(symbols_X, rrc, sig.sps);
tx_Y = upfirdn(symbols_Y, rrc, sig.sps);

% Average waveform sample power
Ptx_X_norm = mean(abs(tx_X).^2);
Ptx_Y_norm = mean(abs(tx_Y).^2);

fprintf('Mean sample power X-pol  = %.6f\n', Ptx_X_norm);
fprintf('Mean sample power Y-pol  = %.6f\n', Ptx_Y_norm);
fprintf('Expected approx.         = %.6f (= 1/sps)\n\n', 1/sig.sps);


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

t_taps = (-(length(rrc)-1)/2 : (length(rrc)-1)/2) / sig.sps;

stem(t_taps, rrc, 'filled');
grid on;
xlabel('Time [symbols]');
ylabel('Amplitude');
title(sprintf('RRC impulse response: rolloff = %.2f, span = %d, sps = %d', ...
    sig.rolloff, sig.span, sig.sps));

% 5.3 RRC frequency response ----------------------------------------------
figure('Name','RRC Frequency Response');

% H = freq response
% f = frequencies in Hz
[H, f] = freqz(rrc, 1, 4096, sig.Fs);

plot(f/1e9, 20*log10(abs(H)/max(abs(H))));
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
[Pxx, f_psd] = pwelch(tx_X, win, noverlap, nfft, sig.Fs, 'centered');
[Pyy, ~]     = pwelch(tx_Y, win, noverlap, nfft, sig.Fs, 'centered');

plot(f_psd/1e9, 10*log10(Pxx/max(Pxx)), 'b');
hold on;
plot(f_psd/1e9, 10*log10(Pyy/max(Pyy)), 'r--');
grid on;

xlabel('Frequency [GHz]');
ylabel('Normalized PSD [dB]');
title(sprintf('DP-QPSK Tx spectrum, expected null-null BW = %.2f GHz', ...
    sig.BW_null/1e9));
legend('X-pol', 'Y-pol', 'Location', 'best');
ylim([-80 5]);

% 5.5 Eye diagrams of pulse-shaped waveform -------------------------------
eye_tx1 = real(tx_X(sig.span*sig.sps+1 : sig.span*sig.sps+4000));
eye_tx2 = real(tx_Y(sig.span*sig.sps+1 : sig.span*sig.sps+4000));
eye_tx = [eye_tx1, eye_tx2];

eyediagram(eye_tx, 2*sig.sps);

% Grab the axes handles from the current figure
ax = findobj(gcf, 'Type', 'axes');

% Note: findobj grabs handles in reverse order of creation!
% ax(1) is the RIGHT plot, ax(2) is the LEFT plot.
title(ax(2), 'X-pol I-component eye diagram');
title(ax(1), 'Y-pol I-component eye diagram');


%% 6. Free Space Optical Channel - Setup

% 6.1 Atmospheric effects -------------------------------------------------
% Beer-Lambert with plane-parallel airmass m = 1/sin(el).
% Valid for el >~ 10 deg; below that, use Kasten-Young airmass.
% Refer to Giggenbach paper
m_air = 1/sin(scen.el_rad);         % airmass
h_atm = scen.t_z ^ m_air;           % linear power transmittance
a_atm_dB = 10*log10(h_atm);         % [dB], negative = loss

fprintf('=== FSO Channel: atmosphere ===\n');
fprintf('Airmass                  = %.3f\n',      m_air);
fprintf('h_atm (linear power)     = %.4f\n',      h_atm);
fprintf('a_atm                    = %.2f dB\n',   a_atm_dB);
fprintf('===============================\n\n');

% 6.2 Turbulence effects --------------------------------------------------
% h starts with the ogs height
h_low  = linspace(scen.h_site_ogs, scen.h_site_ogs + 2e3, 500);     % OGS-2 km
h_mid  = linspace(scen.h_site_ogs + 2e3, 20e3, 500);                % 2–20 km
h_high = linspace(20e3, scen.h_orbit, 500);                         % 20 km-sat. orbit
h_vec = unique([h_low, h_mid, h_high]);                             % MSL height vector [m]
h_agl_vec = h_vec - scen.h_site;                                    % AGL height vector [m]

% Optical / geometry constants
k_opt = 2*pi/scen.wl;                   % optical wavenumber [rad/m]
zeta_rad = pi/2 - scen.el_rad;          % zenith angle [rad]            

% Different models
% HV-5/7. Refer to Carrasco papaer
% 2020_Carrasco_Free-space_optical_links_for_space_communication_n.pdf
prof(1) = struct('name','HV night', 'type','HV', 'A0',1.7e-14, 'v',21, ...
                 'col',[0 0.45 0.74], 'style','-' );
prof(2) = struct('name','HV day', 'type','HV', 'A0',1.7e-13, 'v',30, ...
                 'col',[0.85 0.33 0.10], 'style','-' );
% 2021_Knoedler_Atmospheric Turbulence Statistics and Profile Modeling
prof(3) = struct('name','MHV night', 'type','MHV', 'A0',4.64e-15, 'v',21, ...
                 'col',[0 0.45 0.74], 'style','--');
prof(4) = struct('name','MHV day', 'type','MHV', 'A0',1.7e-13, 'v',30, ...
                 'col',[0.85 0.33 0.10], 'style','--');

% Compute Cn2s
Cn2_all = cell(1, numel(prof));
sigma_R2_all = zeros(1, numel(prof));

for p = 1:numel(prof)
    switch prof(p).type
        case 'HV'               % Carrasco Eq. (8.22), h in AGL
            h_arg = h_agl_vec;
        case 'MHV'              % Knoedler Eq. (1), h + h_s = MSL
            h_arg = h_vec;
    end

    Cn2_all{p} = 0.00594 * (prof(p).v/27)^2 .* (1e-5 .* h_arg).^10 .* exp(-h_arg/1000) ...
               + 2.7e-16 .* exp(-h_arg/1500) ...
               + prof(p).A0 .* exp(-h_arg/100);

    integrand = Cn2_all{p} .* (h_agl_vec - scen.h_ogs).^(5/6);
    sigma_R2_all(p) = 2.25 * k_opt^(7/6) * sec(zeta_rad)^(11/6) * trapz(h_agl_vec, integrand);
end

% Reference profile (MHV night) for downstream plots
Cn2 = Cn2_all{3};
sigma_R2_MHV = sigma_R2_all(3);
sigma_R12_5_MHV = sigma_R2_MHV^(6/5);

% Point-receiver scintillation index (Carrasco / Andrews-Phillips)
% sigma_R2_MHV < 1 => weak turbulence
% Giggenbach 2008 for weak turbulence. sigma_R2 aprox= sigma_I2. CORRECT
sigma_I2_dl = exp( 0.49*sigma_R2_MHV / (1 + 1.11*sigma_R12_5_MHV)^(7/6) ...
                + 0.51*sigma_R2_MHV / (1 + 0.69*sigma_R12_5_MHV)^(5/6) ) - 1;

% Aperture averaging to be introduced
sigma_2_logn = log(sigma_I2_dl + 1);
mu_logn = -0.5 * sigma_2_logn;

% MonteCarlo sampler
sample_logn = @(N) exp(mu_logn + sqrt(sigma_2_logn)*randn(N,1));

% Test values for pdf plot
N_test = 1e5;
I_test = sample_logn(N_test);

% Turbulence channel simulated by lognormal distribution at symbol rate
% We used independent fading. Did not make sense Tc >> Ts.
% We cannot assume so much fading
% The practical BER was better than the theoretical
%h_turb = sample_logn(sig.num_sym);
% Repeat each symbol fading value over the samples-per-symbol
%h_turb_samp = repelem(h_turb, sig.sps);

% Block fading - Slow  log-normal fading
Nblk = ceil(sig.num_sym / ch.nsym_fade);    % number of fading blocks needed

h_turb_blk = sample_logn(Nblk);             % one power gain per block

h_turb = repelem(h_turb_blk, ch.nsym_fade);
h_turb = h_turb(1:sig.num_sym);             % trim to exactly num_sym
h_turb = h_turb(:);                         % forced to column vector to withstand different ch.nsym_fade

h_turb_samp = repelem(h_turb, sig.sps);

% tx_X is longer do to the rrc filter
if length(h_turb_samp) < length(tx_X)
    h_turb_samp = [h_turb_samp; ...
        h_turb_samp(end)*ones(length(tx_X) - length(h_turb_samp), 1)];
else
    h_turb_samp = h_turb_samp(1:length(tx_X));
end

fprintf('=== FSO Channel: turbulence (HV-5/7) ===\n');
fprintf('Wind v (HV-5/7)              = %.1f m/s\n',        scen.v_wind);
fprintf('A0 (ground Cn^2)             = %.2e m^(-2/3)\n',   scen.A0);
fprintf('Zenith angle                 = %.2f deg\n',        rad2deg(zeta_rad));
fprintf('Rytov variance (sigma_R^2)   = %.2f \n',           sigma_R2_MHV);
fprintf('Scint. index (sigma_I^2)     = %.2f \n',           sigma_I2_dl);
fprintf('Monte Carlo lognormal sampler check (N = %.0e):\n', N_test);
fprintf('  mean(I)                    = %.4f (target 1.0000)\n', mean(I_test));
fprintf('  var(I)                     = %.4f (target %.4f)\n', ...
    var(I_test), sigma_I2_dl);
fprintf('========================================\n\n');

% 6.3 Pointing effects ----------------------------------------------------







%% 7. Free Space Optical Channel - Setup Visualization 

% 7.1 Cn2 profiles and Rytov, A0 and wind values --------------------------
figure('Name','Cn2 profile models: HV vs MHV');

hold on; grid on;

set(gca, 'XScale', 'log', 'YScale', 'log');

for p = 1:numel(prof)

    loglog(Cn2_all{p}, max(h_agl_vec,1), ...
        'Color', prof(p).col, ...
        'LineStyle', prof(p).style, ...
        'DisplayName', sprintf(['%s, $A_0 = %.2e\\,\\mathrm{m}^{-2/3}$, ' ...
                      '$v = %.0f\\,\\mathrm{m/s}$'], ...
                      prof(p).name, prof(p).A0, prof(p).v));
end

xlabel('$C_n^2(h)$ [$\mathrm{m}^{-2/3}$]');
ylabel('$h - h_0$ [$\mathrm{m}$]');
title(sprintf('HV / MHV turbulence profiles, $\\epsilon = %.0f^\\circ$', scen.el_deg));
legend('Location','northeast');

xlim([1e-19 1e-12]);
%ylim([1e-3 30]);

% 7.2 Lognormal probability distribution ----------------------------------
figure('Name','Log-normal intensity PDF');

[counts, edges] = histcounts(I_test, 120, 'Normalization','pdf');
centers = 0.5*(edges(1:end-1) + edges(2:end));

I_grid = linspace(1e-3, max(I_test), 2000);
% Theoretical PDF of the lognormal distribution
p_I = 1 ./ (I_grid*sqrt(2*pi*sigma_2_logn)) .* ...
      exp( -(log(I_grid) - mu_logn).^2 ./ (2*sigma_2_logn) );

stairs(centers, counts, 'LineWidth', 1.0); 
hold on; grid on;
semilogy(I_grid, p_I);

set(gca,'YScale','log');

xlabel('$I_n = I/\langle I\rangle$');
ylabel('$p_{I_n}(i)$');
title(sprintf(['Log-normal normalized intensity PDF, ' ...
    '$\\sigma_P^2 = %.3f$, $\\sigma_{\\ln I}^2 = %.4f$'], ...
    sigma_I2_dl, sigma_2_logn));

legend({'MC histogram','Analytical PDF'}, 'Location','northeast');
xlim([0, 1 + 5*sqrt(sigma_I2_dl)]);
ylim([1e-4 2]);


%% 8. Free Space Optical Channel - Simulation

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
% Tx RRC group delay = ((sig.span*sig.sps+1) - 1) / 2
% Rx RRC group delay = ((sig.span*sig.sps+1) - 1) / 2
% Total delay  = (((sig.span*sig.sps+1) - 1) / 2) + (((sig.span*sig.sps+1) - 1) / 2)
total_delay = sig.span * sig.sps;

for i = 1:length(EbN0_dB_vec)
    EbN0 = 10^(EbN0_dB_vec(i)/10);

    % QPSK: Es/N0 = sig.k*Eb/N0, with Es = 1
    EsN0 = sig.k * EbN0;
    N0 = 1 / EsN0;
    sigma = sqrt(N0/2);

    % nI_X + nQ_X = sigma * randn(size(tx_DP(:,1))) + sigma * randn(size(tx_DP(:,1)))
    noise_X = sigma * (randn(size(tx_DP(:,1))) ...
        + 1j*randn(size(tx_DP(:,1))));
    noise_Y = sigma * (randn(size(tx_DP(:,2))) ...
        + 1j*randn(size(tx_DP(:,2))));

    % Apply FSO channel (amplitude scaling = sqrt of power gain)
    % h_ch = h_atm * h_tur * h_poin
    % TO DO
    h_ch_samp = h_atm .* h_turb_samp;

    rx_X_AWGN = tx_DP(:,1) + noise_X;
    rx_Y_AWGN = tx_DP(:,2) + noise_Y;

    rx_X_atm = sqrt(h_atm) .* tx_DP(:,1) + noise_X;
    rx_Y_atm = sqrt(h_atm) .* tx_DP(:,2) + noise_Y;

    rx_X_turb = sqrt(h_atm .* h_turb_samp) .* tx_DP(:,1) + noise_X;
    rx_Y_turb = sqrt(h_atm .* h_turb_samp) .* tx_DP(:,2) + noise_Y;

    rx_X = sqrt(h_ch_samp) .* tx_DP(:,1) + noise_X;
    rx_Y = sqrt(h_ch_samp) .* tx_DP(:,2) + noise_Y;

    % Receiver matched filter uses the same RRC filter
    % No up/downsampling
    rxmf_X_AWGN = upfirdn(rx_X_AWGN, rrc);
    rxmf_Y_AWGN = upfirdn(rx_Y_AWGN, rrc);

    rxmf_X_atm = upfirdn(rx_X_atm, rrc);
    rxmf_Y_atm = upfirdn(rx_Y_atm, rrc);

    rxmf_X_turb = upfirdn(rx_X_turb, rrc);
    rxmf_Y_turb = upfirdn(rx_Y_turb, rrc);

    rxmf_X = upfirdn(rx_X, rrc);
    rxmf_Y = upfirdn(rx_Y, rrc);

    % Downsample at symbol instants
    rx_symbols_X_AWGN = rxmf_X_AWGN(total_delay + 1 : sig.sps : total_delay + sig.num_sym*sig.sps);
    rx_symbols_Y_AWGN = rxmf_Y_AWGN(total_delay + 1 : sig.sps : total_delay + sig.num_sym*sig.sps);

    rx_symbols_X_atm = rxmf_X_atm(total_delay + 1 : sig.sps : total_delay + sig.num_sym*sig.sps);
    rx_symbols_Y_atm = rxmf_Y_atm(total_delay + 1 : sig.sps : total_delay + sig.num_sym*sig.sps);

    rx_symbols_X_turb = rxmf_X_turb(total_delay + 1 : sig.sps : total_delay + sig.num_sym*sig.sps);
    rx_symbols_Y_turb = rxmf_Y_turb(total_delay + 1 : sig.sps : total_delay + sig.num_sym*sig.sps);

    rx_symbols_X = rxmf_X(total_delay + 1 : sig.sps : total_delay + sig.num_sym*sig.sps);
    rx_symbols_Y = rxmf_Y(total_delay + 1 : sig.sps : total_delay + sig.num_sym*sig.sps);

    % Demodulate the gray QPSK into symbols
    ints_hat_X_AWGN = pskdemod(rx_symbols_X_AWGN, sig.M, pi/4, 'gray');
    ints_hat_Y_AWGN = pskdemod(rx_symbols_Y_AWGN, sig.M, pi/4, 'gray');

    ints_hat_X_atm = pskdemod(rx_symbols_X_atm, sig.M, pi/4, 'gray');
    ints_hat_Y_atm = pskdemod(rx_symbols_Y_atm, sig.M, pi/4, 'gray');

    ints_hat_X_turb = pskdemod(rx_symbols_X_turb, sig.M, pi/4, 'gray');
    ints_hat_Y_turb = pskdemod(rx_symbols_Y_turb, sig.M, pi/4, 'gray');

    ints_hat_X = pskdemod(rx_symbols_X, sig.M, pi/4, 'gray');
    ints_hat_Y = pskdemod(rx_symbols_Y, sig.M, pi/4, 'gray');

    % Convert the the integers (0 to 3) into 2-bit pairs (e.g., [1 0]). 
    % MSB by default
    bits_hat_X = int2bit(ints_hat_X, sig.k);
    bits_hat_Y = int2bit(ints_hat_Y, sig.k);

    % Optional, we already have the right dimensions
    bits_hat_X = bits_hat_X(:);
    bits_hat_Y = bits_hat_Y(:);

    % BER
    BER_X(i) = mean(bits_X ~= bits_hat_X);
    BER_Y(i) = mean(bits_Y ~= bits_hat_Y);
    BER_total(i) = mean([bits_X ~= bits_hat_X; bits_Y ~= bits_hat_Y]);
end


%% 9. RX - Visualization

% 9.1 Received spectrum ---------------------------------------------------
figure('Name','Rx Spectrum after Matched Filter - AWGN, atmosphere and turbulence');

% https://es.mathworks.com/help/signal/ref/pwelch.html
[Pxx_r, f_psd_r] = pwelch(rxmf_X, win, noverlap, nfft, sig.Fs, 'centered');
[Pyy_r, ~]     = pwelch(rxmf_Y, win, noverlap, nfft, sig.Fs, 'centered');

plot(f_psd_r/1e9, 10*log10(Pxx_r/max(Pxx_r)), 'b');
hold on;
plot(f_psd_r/1e9, 10*log10(Pyy_r/max(Pyy_r)), 'r--');
grid on;

xlabel('Frequency [GHz]');
ylabel('Normalized PSD [dB]');
title(sprintf('DP-QPSK Rx spectrum after matched filter, expected null-null BW = %.2f GHz', ...
    sig.BW_null/1e9));
legend('X-pol', 'Y-pol', 'Location', 'best');
ylim([-80 5]);

% 9.2 Eye diagrams after channel ------------------------------------------
eye_rx1 = real(rxmf_X(sig.span*sig.sps+1 : sig.span*sig.sps+4000));
eye_rx2 = real(rxmf_Y(sig.span*sig.sps+1 : sig.span*sig.sps+4000));
eye_rx = [eye_rx1, eye_rx2];

eyediagram(eye_rx, 2*sig.sps);

% Grab the axes handles from the current figure
ax2 = findobj(gcf, 'Type', 'axes');

title(ax2(2), 'X-pol I-component eye diagram after matched filter - AWGN, atmosphere and turbulence');
title(ax2(1), 'Y-pol I-component eye diagram after matched filter - AWGN, atmosphere and turbulence');

% 9.3 Received constellations after matched filter -------------------------
figure('Name','Rx Symbols after Matched Filter - AWGN');

subplot(1,2,1);
plot(real(rx_symbols_X_AWGN(1:2000)), imag(rx_symbols_X_AWGN(1:2000)), 'bo', ...
    'MarkerFaceColor','b');
grid on; axis square;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
title('Recovered X-Pol. symbols'); 
xlim([-1.5 1.5]); ylim([-1.5 1.5]);

subplot(1,2,2);
plot(real(rx_symbols_Y_AWGN(1:2000)), imag(rx_symbols_Y_AWGN(1:2000)), 'ro', ...
    'MarkerFaceColor','r');
grid on; axis square;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
title('Recovered Y-Pol. symbols'); 
xlim([-1.5 1.5]); ylim([-1.5 1.5]);

% Constellation AWGN + atmospheric effects
figure('Name','Rx Symbols after Matched Filter - AWGN and atmospheric effects');

subplot(1,2,1);
plot(real(rx_symbols_X_atm(1:2000)), imag(rx_symbols_X_atm(1:2000)), 'bo', ...
    'MarkerFaceColor','b');
grid on; axis square;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
title('Recovered X-Pol. symbols'); 
xlim([-1.5 1.5]); ylim([-1.5 1.5]);

subplot(1,2,2);
plot(real(rx_symbols_Y_atm(1:2000)), imag(rx_symbols_Y_atm(1:2000)), 'ro', ...
    'MarkerFaceColor','r');
grid on; axis square;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
title('Recovered Y-Pol. symbols'); 
xlim([-1.5 1.5]); ylim([-1.5 1.5]);

% Constellation AWGN + atmospheric effects + turbulence effects
figure('Name','Rx Symbols after Matched Filter - AWGN, atmospheric and turbulence effects');

subplot(1,2,1);
plot(real(rx_symbols_X_turb(1:2000)), imag(rx_symbols_X_turb(1:2000)), 'bo', ...
    'MarkerFaceColor','b');
grid on; axis square;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
title('Recovered X-Pol. symbols'); 
xlim([-1.5 1.5]); ylim([-1.5 1.5]);

subplot(1,2,2);
plot(real(rx_symbols_Y_turb(1:2000)), imag(rx_symbols_Y_turb(1:2000)), 'ro', ...
    'MarkerFaceColor','r');
grid on; axis square;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
title('Recovered Y-Pol. symbols'); 
xlim([-1.5 1.5]); ylim([-1.5 1.5]);

% Constellation AWGN + atmospheric effects + turbulence effects
figure('Name','Rx Symbols after Matched Filter - AWGN, atmospheric and turbulence effect');

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

% 9.4 BER table -----------------------------------------------------------
% OpenROADM / ESTOL oFEC pre-FEC BER threshold
BER_preFEC = 2.0e-2;

% Theoretical Eb/N0 at which Gray-QPSK reaches the pre-FEC threshold:
EbN0_thr_lin = erfcinv(2*BER_preFEC)^2;
EbN0_thr_dB  = 10*log10(EbN0_thr_lin);

BER_theory_atm = qfunc(sqrt(2*EbN0_lin.*h_atm));
BER_theory_atm_turb = zeros(size(EbN0_lin));

for i = 1:numel(EbN0_lin)
    BER_theory_atm_turb(i) = mean(qfunc( ...
        sqrt(2*EbN0_lin(i).*h_atm.*h_turb) ));
end

% Print BER table
BER_table = table( ...
    EbN0_dB_vec(:), ...
    BER_theory(:), ...
    BER_theory_atm(:), ...
    BER_theory_atm_turb(:), ...
    BER_X(:), ...
    BER_Y(:), ...
    BER_total(:), ...
    'VariableNames', {'EbN0_dB','BER_AWGN','BER_atm','BER_atm_turb_theory', ...
    'BER_X','BER_Y','BER_total'} );

disp(' ');
disp('=== BER vs Eb/N0 Table ===');
disp(BER_table);

fprintf('\nPre-FEC BER threshold = %.2e\n', BER_preFEC);
fprintf('Theoretical Eb/N0 at pre-FEC threshold = %.3f dB\n\n', EbN0_thr_dB);

% 9.5 BER plot ------------------------------------------------------------
figure('Name','BER vs EbN0_lin');

semilogy(EbN0_dB_vec, BER_theory, 'k--');
hold on;
semilogy(EbN0_dB_vec, BER_theory_atm, 'g--');
semilogy(EbN0_dB_vec, BER_theory_atm_turb, 'c-.');
semilogy(EbN0_dB_vec, BER_X, 'bo-');
semilogy(EbN0_dB_vec, BER_Y, 'rs-');
semilogy(EbN0_dB_vec, BER_total, 'md-');
% Horizontal pre-FEC threshold
yline(BER_preFEC, 'k--', 'Pre-FEC BER $=2.0\times10^{-2}$', ...
    'LineWidth', 1.0, ...
    'LabelHorizontalAlignment','left', ...
    'LabelVerticalAlignment','bottom', ...
    'FontSize', 12, ...
    'Interpreter','latex');
% Vertical theoretical threshold crossing
xline(EbN0_thr_dB, 'k:', ...
    sprintf('$%.2f\\,\\mathrm{dB}$', EbN0_thr_dB), ...
    'LineWidth', 1.0, ...
    'LabelOrientation','horizontal', ...
    'LabelHorizontalAlignment','left', ...
    'LabelVerticalAlignment','bottom', ...
    'FontSize', 12, ...
    'Interpreter','latex');
grid on;
grid minor;
xlabel('Transmitted $E_b/N_0$ before atmospheric attenuation [$\mathrm{dB}$]');
ylabel('BER');
title('Gray-coded DP-QPSK with deterministic atmospheric attenuation + AWGN');
legend('Theory AWGN QPSK', ...
       'Theory atmosphere only', ...
       'Theory atmosphere + lognormal turb.', ...
       'X-Pol.', 'Y-Pol.', 'Total', ...
       'Location','southwest');
ylim([1e-6 1]);
xlim([min(EbN0_dB_vec) max(EbN0_dB_vec)]);
