clear all;
close all;
clc;
warning('off', 'all');

%% =============================================================================================
%% START Parameter Definition

SF = 8; % Spreading factor
M = 2^SF; % no. of samples in one symbol
BW = 125e3; % Bandwidth
Fs = 10e6;  % sampling freq
Ts = 2^SF / BW;   % Symbol period
fc = 915e6;  % carrier center frequency
Power = 14;  % Tx power 14 dB
noise_sigma = 1;
numPreambleSymbols = 8; % Standard LoRa preamble symbol count
Fc = 921.5e6; % spectrum center frequency
cfo = 0;
df = Fc - fc;   % if CFO exists

% RTL-SDR Setup
radio = comm.SDRRTLReceiver('CenterFrequency', fc, ...
    'SampleRate', Fs, ...
    'EnableTunerAGC', true, ...
    'OutputDataType', 'double', ...
    'SamplesPerFrame', 370000);

phy = LoRaPHY(fc, SF, BW, Fs);
phy.has_header = 1;                         % explicit header mode
phy.cr = 1;                                 % code rate = 4/8 (1:4/5 2:4/6 3:4/7 4:4/8)
phy.crc = 1;                                % enable payload CRC checksum
phy.preamble_len = numPreambleSymbols;      % preamble: 8 basic upchirps

%% END Parameter Definition
%% =============================================================================================
%% START Message Reception Loop

while true
    % Receive and process signal
    fprintf('Receiving LoRa signal...\n');
    rxSig = step(radio); % Receive RF signal from RTL-SDR
    received_signal = rxSig;

    % spectrogram for one LoRa packet
    figure
    spectrogram(received_signal, 1000, 0, 1000, Fs, 'yaxis', 'centered')
    ylim([-BW / 2000000, BW / 2000000])
    saveas(gcf,'LoRa_Spectrogram.png')

    i = 1;
    j = 1;
    shift = [];
    final_shift = [];
    decoded_messages = [];
    [shift(i), chirp_len, k] = cross_corr(SF, BW, Fs, cfo, received_signal, phy);
    final_shift(j) = shift(i);
    [decoded_message, ~] = phy.demodulate(received_signal(final_shift(j):end));
    decoded_messages = {decoded_message(:, 1).'};

    % Loop over received signal, run through AlignTrack algorithm, then use
    % output to advance to next peak corresponding to next window alignment
    % point
    while (shift(i) + chirp_len - 1) <= length(received_signal)
        % received signal after windowing and FFT
        [f, received_fft] = LoRa_demod_1(received_signal(shift(i):shift(i) + chirp_len - 1), fc, SF, BW, Fs, cfo, "hann");
        AA = AlignTrack(received_fft);
        
        if length(AA) > 1
            [B, IX] = sort(AA); %order the amplitudes
            f2 = f(IX(end - 1)); %frequency of second peak
        else
            break
        end

        i = i + 1;
        shift(i) = shift(i - 1) + abs(Fs * ((f2 - BW) / k));
        [decoded_message, ~] = phy.demodulate(received_signal(shift(i):end));

        % Append the new row if it's not a duplicate
        if (~isempty(decoded_message)) && (~any(cellfun(@(row) isequal(row, decoded_message(:, 1).'), decoded_messages)))
            decoded_messages{end + 1} = decoded_message(:, 1).';
            final_shift(end + 1) = shift(i);
        end
    end

    for i = 1:length(decoded_messages)
        fprintf('Row %d: ', i);
        disp(decoded_messages{i}); % Display
        [decoded_data, ~] = phy.decode(decoded_messages{i}.');
        fprintf('Received Message: %s\n', char(decoded_data.'));
    end
end

%% END Message Reception Loop
%% =============================================================================================