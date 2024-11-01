clear all;
close all;
clc;

SDR = 0;
AlignTrack = 1;
collision = 1;

% LoRa Parameters
SF = 8; % Spreading factor
M=2^SF; % no. of samples in one symbol
BW = 125e3 ; % Bandwidth
Fs = 1e6;  % sampling freq
fc = 915e6;  % carrier center frequency
noise_sigma=0;
message = "Hello world!";
message2 = "Hello AlignTrack!";
numPreambleSymbols = 8; % Standard LoRa preamble symbol count
symbolRate = BW / (2^SF); % Symbol rate
shift = 0;
payload_offset = 5000;

if SDR == 1
    % RTL-SDR Setup
    radio = comm.SDRRTLReceiver('CenterFrequency', fc, ...
        'SampleRate', Fs, ...
        'EnableTunerAGC', true, ...
        'OutputDataType', 'double', ...
        'SamplesPerFrame', 370000);
end

% LoRa Setup
phy = LoRaPHY(fc, SF, BW, Fs);
phy.has_header = 1;                         % explicit header mode
phy.cr = 1;                                 % code rate = 4/8 (1:4/5 2:4/6 3:4/7 4:4/8)
phy.crc = 1;                                % enable payload CRC checksum
phy.preamble_len = numPreambleSymbols;      % preamble: 8 basic upchirps

while true
    if SDR == 1
        % Receive and process signal
        fprintf('Receiving LoRa signal...\n');
        rxSig = step(radio); % Receive RF signal from RTL-SDR
        received_signal = rxSig;
    else
        if collision == 1
            % Encode payload
            fprintf("[encode] message:\n");
            disp(double(char(message)).')
            symbols = phy.encode(double(char(message)).');
            fprintf("[encode] symbols:\n");
            disp(symbols);

            % Encode payload2
            fprintf("[encode] message2:\n");
            disp(double(char(message2)).')
            symbols2 = phy.encode(double(char(message2)).');
            fprintf("[encode] symbols2:\n");
            disp(symbols2);

            if find([size(symbols), size(symbols2)]==max(max([size(symbols), size(symbols2)]))) == 1
                % Baseband Modulation
                signalIQ1 = [phy.modulate(symbols); zeros(payload_offset, 1)];

                % Baseband Modulation2
                signalIQ2 = [zeros(numel(signalIQ1) - numel(phy.modulate(symbols2)), 1); phy.modulate(symbols2)];
            else
                % Baseband Modulation2
                signalIQ2 = [zeros(payload_offset, 1); phy.modulate(symbols2)];

                % Baseband Modulation
                signalIQ1 = [phy.modulate(symbols); zeros(numel(signalIQ2) - numel(phy.modulate(symbols)), 1)];
            end

            signalIQ3 = signalIQ1 + signalIQ2;
            figure(1)
            spectrogram(signalIQ3,1000,0,1000,Fs,'yaxis','centered')
            noise=noise_sigma*randn(length(signalIQ3),1);
            received_signal=signalIQ3+noise;
        else
            % Encode payload
            fprintf("[encode] message:\n");
            disp(double(char(message)).')
            symbols = phy.encode(double(char(message)).');
            fprintf("[encode] symbols:\n");
            disp(symbols);

            % Baseband Modulation
            signalIQ1 = phy.modulate(symbols);
            figure(1)
            spectrogram(signalIQ1,1000,0,1000,Fs,'yaxis','centered')
            noise=noise_sigma*randn(length(signalIQ1),1);
            received_signal=signalIQ1+noise;
        end
    end

    if AlignTrack == 0
        new_received_signal = received_signal;
        figure(2)
        spectrogram(new_received_signal,1000,0,1000,Fs,'yaxis','centered')

        % Demodulation
        length(new_received_signal)
        [symbols_d, cfo] = phy.demodulate(new_received_signal);
        fprintf("[demodulate] symbols:\n");
        disp(symbols_d);

        % Decoding
        [data, checksum] = phy.decode(symbols_d);
        fprintf("[decode] data:\n");
        disp(data(1:end-2));
        fprintf("[decode] checksum:\n");
        disp(checksum);
    else
        %% AlignTrack Implementation
        figure(2)
        spectrogram(received_signal,1000,0,1000,Fs,'yaxis','centered')
        [symbols_d, cfo] = phy.demodulate(received_signal);
        [received_fft] = fft(symbols_d);
        new_received_fft=received_fft;

        % peak extraction algorithm with AlignTrack decoding for complete packet
        [m n]=size(received_fft);

        % Message is unable to be decoded i.e. size(received_fft) == 0, try
        % AlignTrack deconfliction
        if m == 0 || n == 0
            fprintf("Deconfliction Necessary\n")
            [received_fft] = phy.LoRa_demod_1(payload_offset);
            new_received_fft=received_fft;
            [m n]=size(received_fft);
        end

        for j1=0:1:m-1
            for i1=1:1:n
                [val(j1+1),idx(j1+1)]=max(new_received_fft(j1+1,:));
                r(j1+1)=mean(new_received_fft(j1+1,:))+6*std(new_received_fft(j1+1,:));
                if val(j1+1)>=r(j1+1)
                    I(j1+1,i1)=idx(j1+1);
                    before=(idx(j1+1)-1);
                    after=(idx(j1+1)+1);
                    % finding local min before index value
                    ii1=1;
                    while ii1<n
                        if before==0
                            before=idx(j1+1);
                        end
                        if before-1<=0
                            before=idx(j1+1)+1;
                        end
                        if received_fft(before+1)>received_fft(before) && received_fft(before)<received_fft(before-1)
                            b=before;
                            before=before-1;
                            break;
                        else
                            before=before-1;
                        end
                        ii1=ii1+1;
                    end
                    % finding local min after index value
                    ii2=1;
                    while ii2<n
                        if after>=n
                            after=idx(j1+1);
                        end
                        if after+1>=n
                            after=idx(j1+1)-1;
                        end
                        if received_fft(after)<received_fft(after-1) && received_fft(after+1)>received_fft(after)
                            a=after;
                            after=after+1;
                            break;
                        else
                            after=after+1;
                        end
                        ii2=ii2+1;
                    end
                    %  remove local min points before and after index value
                    for y=before:after
                        new_received_fft(j1+1,y+1)=0; % changed this line from "y" to "y + 1"
                    end
                else
                    break;
                end
            end
        end
        [m1 n1]=size(I);
        % sort the index value for further processing
        I=sort(I,2);
        % side lobe elimination
        for ee=1:1:m1
            for kk=1:1:length(I(ee,:))
                if I(ee,kk)==0
                    continue;
                else
                    d1=kk;
                    break;
                end
            end
            A.(sprintf('RandomVariable_%d', ee))=I(ee,[d1:end]);
            AA=A.(sprintf('RandomVariable_%d', ee));
            issidelobe=zeros(1,length(AA));
            for i=1:1:length(AA)
                for j=i+1:1:length(AA)
                    for k=1:1:length(AA)
                        if AA(j)-AA(i)==AA(i)-AA(k)
                            if k~=j && received_fft(AA(k))==received_fft(AA(j))
                                issidelobe(k)=1;
                                issidelobe(j)=1;
                            end
                        end
                    end
                end
            end
            AA=AA(issidelobe~=1);
        end
        plot(received_fft)
        symbols_d = ifft(received_fft(:,1));

        for payload_num=1:1:size(received_fft,2)
            symbols_d = ifft(received_fft(:,payload_num));

            % Decoding
            [data, checksum] = phy.decode(symbols_d);
            fprintf("[decode] data:\n");
            disp(data(1:end-2));
            fprintf("[decode] checksum:\n");
            disp(checksum);
        end
    end
end