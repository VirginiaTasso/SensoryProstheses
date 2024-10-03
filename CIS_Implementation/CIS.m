%% Read the mp3 file
filename = 'Sam_excerpt.wav'; % or change with Sam_excerpt.wav
[audio_sgn, fs] = audioread(filename); % this is a stereo audio -> 2 channels
audio_sgn = audio_sgn(:, 1); % select one channel
%% if audio signal select just 20 seconds
if strcmp(filename, 'time.mp3')
 % The audio is long, just select 10 seconds, from min 2:14
start_time = 134; %desired starting second
end_time = 154; %  desired ending second

% need to convert seconds in sample idxs
start_sample = round(start_time * fs);
end_sample = round(end_time * fs);

sgn = audio_sgn(start_sample:end_sample, 1); % select just one channel
len = length(sgn);
else 
    sgn = audio_sgn;
    len = length(sgn);
end


%% Plot signal

figure;
plot(sgn, 'linewidth', 1.5);
title('Plot of the original signal', 'FontSize', 18)
xlabel('Time [s]', 'FontSize', 18);
ylabel('Signal', 'FontSize', 18)
grid on

%%
any(isnan(audio_sgn)) 
any(isinf(audio_sgn)) 

%% ========= Compute Spectrum ========= %%
spect = fft(sgn);
f = (0:length(spect)-1)*fs/length(spect);
% Plot Spectrum
figure;
plot(f(1:length(f)/2), abs(spect(1:length(f)/2)), 'linewidth', 1.5)
title('Spectrum of the signal','FontSize', 20)
xlabel('Frequency [Hz]', 'FontSize', 18)
ylabel('Magnitude','FontSize',18)
grid on

%% ========= Select different frequency bands ========= %%
f_bands = [100 400; 400 1000; 1000 3000; 3000 8000]; % four music
%f_band = [20 200; 200 500; 500 2000; 2000 5000]; % for speech
f_nyq = fs /2; % Nyquist frequency -> needed because the cut off frequency must be normalized
filt_sgn = {}; % initialize cell array to store the filtered signal


for i=1:length(f_bands)
    f_band = f_bands(i,:) ./ f_nyq;
    [b, a] = butter(4, f_band, "bandpass");
    filt_sgn{i} = filtfilt(b, a, sgn);
end

% Plot the filtered signal
figure;
for i=1:length(f_bands)
    subplot(2, 2, i)
    plot(filt_sgn{i})
    title(['Frequency Band ', num2str(f_bands(i, 1)), '-', num2str(f_bands(i, 2)), ' Hz',], 'FontSize', 18);
    xlabel('Time [s]', 'FontSize', 18)
    ylabel('Amplitude', 'FontSize', 18)
    grid on
end
sgtitle('Filtered Signal Divided into Different Frequency Bands', 'FontSize', 18);

%% ========= Envelope Extraction ========= %%
% rectification of the filtered signal -> extraction of the abosulut value
rectified = {}; % initialization of the cell array to store the rectified signal
for i=1:length(filt_sgn)
    rectified{i} = abs(filt_sgn{i});
end

% Now extract the envelop with a LP filter
envelope = {}; % initialize cell array to store the envelope
lp_cutoff = 400 / f_nyq;
[b, a] = butter(4, lp_cutoff, "low");
% Extract the envelop from each band
for i=1:length(f_bands)
    envelope{i} = filtfilt(b, a, rectified{i});
end

% Plot the different envelopes

figure;
for i=1:length(f_bands)
    subplot(2, 2, i)
    plot(envelope{i})
    title(['Frequency Band ', num2str(f_bands(i, 1)), '-', num2str(f_bands(i, 2)), ' Hz']);
    xlabel('Time [s]', 'FontSize', 18)
    ylabel('Amplitude', 'FontSize',18)
    grid on
end
sgtitle('Envelopes of the different bands');

%% ========= Logarithmic Compression ========= %%
compressed_sgn = {}; % initialize cell array 
% Compress each band
for i=1:length(f_bands)
    compressed_sgn{i} = log_compression(envelope{i}, 1);
end

figure;
for i=1:length(f_bands)
    subplot(2,2,i);
    plot(compressed_sgn{i}, 'LineWidth',1.5)
    title(['Frequency Band ', num2str(f_bands(i, 1)), '-', num2str(f_bands(i, 2))])

    xlabel('Time [s]', 'FontSize',18);
    ylabel('Amplitude (dB)', 'FontSize',18);
    grid on
end
sgtitle('Compressed Signals', 'FontSize', 18)

%% ========= Biphasic pulses ========= %%
pulse_rate = 800; % pps
time_btw_pulses = 1 / pulse_rate;
samples_btw_pulses = round(time_btw_pulses * fs);
modulated_signals = cell(length(f_bands), 1); % initialize cell array to store pulse trains for each band

for i = 1:length(f_bands)
    curr_sgn = compressed_sgn{i};
    count = 1;
    modulated_signals{i} = zeros(size(curr_sgn));
    for j=1:samples_btw_pulses:length(curr_sgn)
        if mod(count, 2) == 0
            modulated_signals{i}(j) = curr_sgn(j);
        else
            modulated_signals{i}(j) = -curr_sgn(j);
        end
        count = count + 1;
    end
end
%Combine Pulse Trains to Recreate Audio this is what our goal is 
final_pulses = zeros(size(sgn));  %(Initialize the final pulse signal)
for i = 1:length(f_bands)
    final_pulses = final_pulses + modulated_signals{i}; % ( Sum the pulse trains from all bands)
end

%% Plot Biphasic Signals for Each Band
time = (0:length(sgn)-1) / fs;
colors = ['r', 'g', 'b', 'm'];
figure;

% Plot separated bands
for i = 1:length(f_bands)
    subplot(length(f_bands)+1, 1, i); %this is to create a subplot for each band
    plot(time, modulated_signals{i}, 'color', colors(i)); %to plot a portion of the signal
    title(['Biphasic Pulses for Band ' num2str(i) ' (' num2str(f_bands(i,1)) '-' num2str(f_bands(i,2)) ' Hz)'], 'FontSize',18);
    xlabel('Time [s]', 'FontSize',18);
    ylabel('Amplitude', 'FontSize',18);
    grid on
end

% plot the final combined signal
subplot(length(f_bands)+1, 1, length(f_bands)+1);
plot(final_pulses);
title('Reconstruced Signal', 'FontSize',18);
xlabel('Time [s]', 'FontSize',18);
ylabel('Amplitude', 'FontSize', 18);
grid on
%

figure;
% Plot separated bands overlapped
seconds_to_plot = 5;
for i = 1:length(f_bands)
    hold on;
    plot(time(1:seconds_to_plot*fs), modulated_signals{i}(1:seconds_to_plot*fs), 'color', colors(i)); %to plot a portion of the signal
    title('Overlapped Pulses');
    xlabel('Time [s]', 'FontSize', 18);
    ylabel('Amplitude', 'FontSize',18);
    legend(arrayfun(@(b) ['Band ' num2str(b) ' (' num2str(f_bands(b,1)) '-' num2str(f_bands(b,2)) ' Hz)'], 1:length(f_bands), 'UniformOutput', false));
    grid on
end

figure;
plot(time, final_pulses);
title('Reconstructed Signal', 'FontSize', 18);
xlabel('Time [s]', 'FontSize',18);
ylabel('Amplitude', 'FontSize',18);
grid on

%% Plot in the frequency
% Calculate the frequency axis
N = length(time); % Number of points in the signal
frequencies = (0:N-1)*(fs/N); % Frequency axis (in Hz)
colors = ['r', 'g', 'b', 'm'];

figure;

% Plot spectra for each band
for i = 1:length(f_bands)
    % Perform FFT
    spectrum = abs(fft(modulated_signals{i}));
    
    % Plot spectrum
    subplot(length(f_bands)+1, 1, i);
    plot(frequencies, spectrum, 'LineWidth', 1.5, 'color', colors(i));
    xlim([0, max(f_bands(:))*1.5]); % Limit to the relevant frequency range
    title(['Spectrum of Band ' num2str(i) ' (' num2str(f_bands(i,1)) '-' num2str(f_bands(i,2)) ' Hz)']);
    xlabel('Frequency [Hz]', 'FontSize',18);
    ylabel('Amplitude', 'FontSize',18);
    ylim([0 1500])
    grid on;
end

% Plot spectrum of the final combined signal
spectrum_combined = abs(fft(final_pulses));
subplot(length(f_bands)+1, 1, length(f_bands)+1);
plot(frequencies, spectrum_combined, 'LineWidth', 1.5);
xlim([0, max(f_bands(:))*1.5]); % Limit to the relevant frequency range
title('Spectrum of Reconstructed Signal');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
grid on;

figure;

% Plot overlapped spectra
for i = 1:length(f_bands)
    hold on;
    spectrum = abs(fft(modulated_signals{i}));
    plot(frequencies, spectrum, 'color', colors(i)); % Plot the spectrum for each band
end
xlim([0, max(f_bands(:))*1.5]); % Limit to the relevant frequency range
title('Overlapped Spectra','FontSize',18);
xlabel('Frequency [Hz]', 'FontSize',18);
ylabel('Amplitude', 'FontSize',18);
legend(arrayfun(@(b) ['Band ' num2str(b) ' (' num2str(f_bands(b,1)) '-' num2str(f_bands(b,2)) ' Hz)'], 1:length(f_bands), 'UniformOutput', false));
grid on;

% Plot spectrum of the final combined signal
figure;
plot(frequencies, spectrum_combined);
xlim([0, max(f_bands(:))*1.5]); % Limit to the relevant frequency range
title('Spectrum of Reconstructed Signal', 'FontSize',18);
xlabel('Frequency [Hz]', 'FontSize',18);
ylabel('Amplitude', 'FontSize',18);
grid on;

%% Save signal

%audiowrite('recreated_sound_speech.wav', final_pulses, fs); % change with "recreated_sound_speech"

