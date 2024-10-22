clear all;
close all;
clc;
SF = 8; % Spreading factor
M=2^SF; % no. of samples in one symbol
BW = 125e3 ; % Bandwidth
Fs = 10e6 ;  % sampling freq
fc = 915e6 ;  % carrier center frequency
noise_sigma=1;
message = "Hello world!" ;
Fc = 921.5e6; % spectrum center frequency
df=Fc-fc;   % if CFO exists
shift = 0;

phy = LoRaPHY(fc, SF, BW, Fs);
phy.has_header = 1;         % explicit header mode
phy.cr = 1;                 % code rate = 4/8 (1:4/5 2:4/6 3:4/7 4:4/8)
phy.crc = 1;                % enable payload CRC checksum
phy.preamble_len = 8;       % preamble: 8 basic upchirps

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

%% AlignTrack Implementation
%[received_fft] = fft(received_signal);
[received_fft] = LoRa_demod_1(received_signal,SF,BW,Fs,shift); % changed this line from "LoRa_demod_1_new" to "LoRa_demod_1"
new_received_fft=received_fft;

% peak extraction algorithm with AlignTrack decoding for complete packet
[m n]=size(received_fft);
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

new_received_signal = ifft(new_received_fft(:,1));
figure(3)
spectrogram(new_received_signal,1000,0,1000,Fs,'yaxis','centered')

% Demodulation
[symbols_d, cfo] = phy.demodulate(new_received_signal);
fprintf("[demodulate] symbols:\n");
disp(symbols_d);

% Decoding
[data, checksum] = phy.decode(symbols_d);
fprintf("[decode] data:\n");
disp(data);
fprintf("[decode] checksum:\n");
disp(checksum);

str = char(data(1:end-2));
disp(str.');