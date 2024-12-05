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
noise_sigma = 0;
numPreambleSymbols = 8; % Standard LoRa preamble symbol count
message1 = "Hello World!";
message2 = "Hello AlignTrack!";
message3 = "Hello Goodbye!";
Fc = 921.5e6; % spectrum center frequency
cfo = 0;
df = Fc - fc;   % if CFO exists

phy = LoRaPHY(fc, SF, BW, Fs);
phy.has_header = 1;                         % explicit header mode
phy.cr = 1;                                 % code rate = 4/8 (1:4/5 2:4/6 3:4/7 4:4/8)
phy.crc = 1;                                % enable payload CRC checksum
phy.preamble_len = numPreambleSymbols;      % preamble: 8 basic upchirps

%% END Parameter Definition
%% =============================================================================================
%% START Symbol Generation

% Encode payload 1
fprintf("[encode] message 1:\n");
disp(double(char(message1)).')
symbols1 = phy.encode(double(char(message1)).');
fprintf("[encode] symbols 1:\n");
disp(symbols1);

% Encode payload 2
fprintf("[encode] message 2:\n");
disp(double(char(message2)).')
symbols2 = phy.encode(double(char(message2)).');
fprintf("[encode] symbols 2:\n");
disp(symbols2);

% Encode payload 3
fprintf("[encode] message 3:\n");
disp(double(char(message3)).')
symbols3 = phy.encode(double(char(message3)).');
fprintf("[encode] symbols 3:\n");
disp(symbols3);

offset_time = round(min([(length(phy.modulate(symbols1)) / Fs) * 1000, (length(phy.modulate(symbols2)) / Fs) * 1000, (length(phy.modulate(symbols3)) / Fs) * 1000])); % ms

windows = {'rectangular', 'welch', 'sine', 'hann', 'hamming'};
tEnds = [];
SER = [];
cum_SER = [];
new_SER = [];
trials = 10;

for trial = 1:1:trials
    trial
    % signal recieved at LoRa gateway
    received_signal = send_receive_signal(offset_time, symbols1, symbols2, symbols3, Fs, Ts, phy, noise_sigma);

    % spectrogram for one LoRa packet
    if trial == trials
        figure
        spectrogram(received_signal, 1000, 0, 1000, Fs, 'yaxis', 'centered')
        ylim([-BW / 2000000, BW / 2000000])
        saveas(gcf,'LoRa_Spectrogram.png')
    end

%% END Symbol Generation
%% =============================================================================================
%% START AlignTrack Trials

    for window = 1:1:length(windows)
        tStart = tic;

        i = 1;
        j = 1;
        shift = [];
        final_shift = [];
        decoded_messages = [];
        [shift(i), chirp_len, k] = cross_corr(SF, BW, Fs, cfo, received_signal, phy);
        final_shift(j) = shift(i);
        while true
            try
                [decoded_message, ~] = phy.demodulate(received_signal(final_shift(j):end));
                decoded_messages = {decoded_message(:, 1).'};
                break
            catch            
                % signal recieved at LoRa gateway
                received_signal = send_receive_signal(offset_time, symbols1, symbols2, symbols3, Fs, Ts, phy, noise_sigma);
            
                % spectrogram for one LoRa packet
                if trial == trials
                    figure
                    spectrogram(received_signal, 1000, 0, 1000, Fs, 'yaxis', 'centered')
                    ylim([-BW / 2000000, BW / 2000000])
                    saveas(gcf,'LoRa_Spectrogram.png')
                end

                tStart = tic;

                i = 1;
                j = 1;
                [shift(i), chirp_len, k] = cross_corr(SF, BW, Fs, cfo, received_signal, phy);
                final_shift(j) = shift(i);
            end
        end

        % Loop over received signal, run through AlignTrack algorithm, then use
        % output to advance to next peak corresponding to next window alignment
        % point
        while (shift(i) + chirp_len - 1) <= length(received_signal)
            % received signal after windowing and FFT
            [f, received_fft] = LoRa_demod_1(received_signal(shift(i):shift(i) + chirp_len - 1), fc, SF, BW, Fs, cfo, string(windows(window)));
            AA = AlignTrack(received_fft);
            
            if length(AA) > 1
                [B, IX] = sort(AA); %order the amplitudes
                % A1 = B(end); %amplitude of second peak
                % A2 = B(end - 1); %amplitude of second peak
                % f1 = f(IX(end)); %frequency of first peak
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

        final_shift

        % Loop through results, cross-correlate with the inputs to determine
        % closest match, then compare to determine SER
        symbol_corr = [];
        for i = 1:length(decoded_messages)
            fprintf('Row %d: ', i);
            disp(decoded_messages{i}); % Display the contents of the cell
            [decoded_data, ~] = phy.decode(decoded_messages{i}.');
            fprintf('Received Message: %s\n', char(decoded_data.'));
            max_num = 0;

            % Loop through each output row
            for j = 1:3
                switch j
                    case 1
                        symbols = symbols1;
                    case 2
                        symbols = symbols2;
                    case 3
                        symbols = symbols3;
                end
                correlation = abs(xcorr(symbols, decoded_messages{i}));

                if max(correlation) > max_num
                    max_num = max(correlation);
                    symbol_corr(i) = j;
                end
            end

            switch symbol_corr(i)
                case 1
                    symbols = symbols1;
                case 2
                    symbols = symbols2;
                case 3
                    symbols = symbols3;
            end

            % Calculate the number of symbol errors
            uniqueToArray1 = setdiff(symbols.', decoded_messages{i});
            uniqueToArray2 = setdiff(decoded_messages{i}, symbols.');
            numErrors = length(uniqueToArray1) + length(uniqueToArray2);

            % Calculate the symbol error rate (SER)
            SER(window, trial, i) = numErrors / length(symbols.');

            % Display results
            fprintf('Correlated Input: %d\n', symbol_corr(i));
            fprintf('Number of symbol errors: %d\n', numErrors);
            fprintf('Symbol Error Rate (SER): %.4f\n', SER(window, trial, i));

        end
        % Calculate the cumulative SER
        if trial == 1
            cum_SER(window, trial) = sum(SER(window, trial, :));
        else
            cum_SER(window, trial) = cum_SER(window, trial - 1) + sum(SER(window, trial, :));
        end
        tEnds(window, trial) = toc(tStart);
    end
end

%% END AlignTrack Trials
%% =============================================================================================
%% START Plot/Display Results

windows
mean(tEnds, 2).'

% Plot time for each window by trial
figure
hold on
title("Time vs. Window by Trial")
for i = 1:1:size(tEnds, 1)
    x = 1:1:trials;
    plot(x, tEnds(i, :))
    legend(windows)
end
saveas(gcf,'Time_vs_Window_by_Trial.png')
hold off

% Plot SER for each window by trial
figure
hold on
title("SER vs. Window by Trial")
for i = 1:1:size(SER, 1)
    x = 1:1:trials;
    plot(x, sum(SER(i, :, :), 3))
    legend(windows)
end
saveas(gcf,'SER_vs_Window_by_Trial.png')
hold off

% Plot cumulative SER for each window
figure
hold on
title("Cumulative SER vs. Window")
for i = 1:1:size(cum_SER, 1)
    x = 1:1:trials;
    plot(x, cum_SER(i, :))
    legend(windows)
end
saveas(gcf,'Cumulative_SER_vs_Window.png')
hold off

%% END Plot/Display Results
%% =============================================================================================