
% =========================================
% Experiment 1: DFT vs DTFT
% =========================================

% -----------------------------------------
% Part 1(a): Create a MATLAB Function for Computing the DFT
% -----------------------------------------

% Function definition for DFT
function Xk = calculate_dft(x)
    N = length(x);
    n = 0:N-1;
    k = n;  % Frequency indices
    % Create the DFT matrix
    W = exp(-1j * 2 * pi / N * (n') * k);
    % Compute the DFT
    Xk = W' * x(:);
end

% Note: The DTFT function calculates the sum over n for each frequency w.
% For the DFT, w = 2*pi*k/N for k = 0:N-1.

% -----------------------------------------
% Part 1(b): Compute DFTs of Rectangular Signals with Different Lengths
% -----------------------------------------

% Define lengths N
N_values = [16, 32, 64, 128, 256];

% Initialize cell arrays to store signals and their DFTs
signals = cell(1, length(N_values));
DFTs = cell(1, length(N_values));

% Loop over each N
for idx = 1:length(N_values)
    N = N_values(idx);
    % Create the rectangular signal
    x = zeros(1, N);
    x(1:16) = 1;  % Rectangular pulse of length 16
    signals{idx} = x;
    % Compute the DFT
    Xk = calculate_dft(x);
    DFTs{idx} = Xk;
    % Plot the magnitude of the DFT
    figure(1);
    subplot(length(N_values), 1, idx);
    stem(0:N-1, abs(Xk));
    title(['Magnitude of DFT for N = ', num2str(N)]);
    xlabel('Frequency Index k');
    ylabel('|X[k]|');
end

f = gcf;
exportgraphics(f,'1d_homemade.png','Resolution',1200)

% -----------------------------------------
% Part 1(d): Compute and Plot the DTFT of the Same Input Set
% -----------------------------------------

% Function for DTFT (from last lab)
function Xw = calculate_dtft(x, w)
    n = 0:length(x)-1;
    Xw = arrayfun(@(omega) sum(x .* exp(-1j * omega * n)), w);
end

% Compute and plot DTFTs
for idx = 1:length(N_values)
    N = N_values(idx);
    x = signals{idx};
    % Frequency vector w matching the DFT frequencies
    w = (2 * pi / N) * (0:N-1);
    % Compute DTFT
    Xw = calculate_dtft(x, w);
    % Plot the magnitude of the DTFT
    figure(2);
    subplot(length(N_values), 1, idx);
    plot(w, abs(Xw));
    title(['Magnitude of DTFT for N = ', num2str(N)]);
    xlabel('Frequency \omega (rad/sample)');
    ylabel('|X(e^{j\omega})|');
end
f = gcf;
exportgraphics(f,'1d_builtin.png','Resolution',1200)


% -----------------------------------------
% Part 1(e): Use fft() to Compute the DFT and Compare Results
% -----------------------------------------

% Compute and compare fft() results
for idx = 1:length(N_values)
    N = N_values(idx);
    x = signals{idx};
    % Compute DFT using fft()
    Xk_fft = fft(x);
    % Compare with custom DFT function
    Xk_custom = DFTs{idx};
    % Plot comparison
    figure(3);
    subplot(length(N_values), 1, idx);
    stem(0:N-1, abs(Xk_fft), 'b', 'DisplayName', 'fft()');
    hold on;
    stem(0:N-1, abs(Xk_custom), 'r--', 'DisplayName', 'Custom DFT');
    title(['Comparison of DFT Methods for N = ', num2str(N)]);
    xlabel('Frequency Index k');
    ylabel('|X[k]|');
    legend;
    hold off;
end
f = gcf;
exportgraphics(f,'1e_comparison.png','Resolution',1200)
%% 

% =========================================
% Experiment 2: Analyzing a Synthesized Vowel Sound
% =========================================

% -----------------------------------------
% Part 2(a): Load the .wav File and Determine Sampling Frequency
% -----------------------------------------

% Load the audio file 'aaa.wav'
[y, Fs] = audioread('aaa.wav');

% Display the sampling frequency
disp(['Sampling Frequency Fs = ', num2str(Fs), ' Hz']);

% -----------------------------------------
% Part 2(b): Plot the Waveform of the First 300 Samples
% -----------------------------------------

% Extract the first 300 samples
y_segment = y(1:300);

% Plot the waveform
figure(4);
plot(0:299, y_segment);
title('Waveform of the First 300 Samples');
xlabel('Sample Index n');
ylabel('Amplitude');
f = gcf;
exportgraphics(f,'2b.png','Resolution',1200)

% period is about 80
% can also use autocorrelation to find the period

% -----------------------------------------
% Part 2(c): Compare DFT Magnitudes of One and Two Periods
% -----------------------------------------

T0 = 80;  % period determined from 2b
one_period = y(1:T0);
two_periods = y(1:2*T0);

% Compute DFTs
X1 = fft(one_period);
X2 = fft(two_periods);

% Plot the magnitudes
figure(5);
subplot(2,1,1);
stem(0:length(X1)-1, abs(X1));
title('Magnitude of DFT of One Period');
xlabel('Frequency Index k');
ylabel('|X1[k]|');

subplot(2,1,2);
stem(0:length(X2)-1, abs(X2));
title('Magnitude of DFT of Two Periods');
xlabel('Frequency Index k');
ylabel('|X2[k]|');
f = gcf;
exportgraphics(f,'2c.png','Resolution',1200)

% -----------------------------------------
% Part 2(d): Zero-Pad the Signal Pieces to 1024 Points and Plot DFTs
% -----------------------------------------

% Zero-pad to 1024 points
N_fft = 1024;
one_period_zp = [one_period; zeros(N_fft - length(one_period), 1)];
two_periods_zp = [two_periods; zeros(N_fft - length(two_periods), 1)];

% Compute DFTs
X1_zp = fft(one_period_zp);
X2_zp = fft(two_periods_zp);

% Compute frequency vector
f = (0:N_fft-1) * Fs / N_fft;

% Plot the magnitudes
figure(6);
subplot(2,1,1);
stem(f, abs(X1_zp));
title('Magnitude of DFT of One Period (Zero-Padded to 1024 Points)');
xlabel('Frequency (Hz)');
ylabel('|X1\_zp(f)|');

subplot(2,1,2);
stem(f, abs(X2_zp));
title('Magnitude of DFT of Two Periods (Zero-Padded to 1024 Points)');
xlabel('Frequency (Hz)');
ylabel('|X2\_zp(f)|');
f = gcf;
exportgraphics(f,'2d.png','Resolution',1200)

% -----------------------------------------
% Part 2(e): Compute and Plot DFTs for 1 to 5 Periods
% -----------------------------------------

% Compute DFTs for 1 to 5 periods 
figure(7);
for k = 1:5
    % Extract k periods of the signal
    signal = y(1:k*T0);
    N_signal = length(signal);
    
    % Compute DFT 
    Xk = fft(signal);
    
    % Compute frequency vector
    f = (0:N_signal-1) * Fs / N_signal;
    
    % Plot magnitude of the DFT
    subplot(5,1,k);
    stem(f, abs(Xk));
    title(['Magnitude of DFT of ', num2str(k), ' Period(s)']);
    xlabel('Frequency (Hz)');
    ylabel(['|X', num2str(k), '(f)|']);
    grid on;
end
f = gcf;
exportgraphics(f,'2e.png','Resolution',1200)
