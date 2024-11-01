clear all;
close all;
clc;
SF = 8; % Spreading factor
M=2^SF; % no. of samples in one symbol
BW = 125e3 ; % Bandwidth
Fs = 10e6 ;  % sampling freq
Ts=2^SF/BW;   % Symbol period
fc = 915e6 ;  % carrier center frequency
Power = 14 ;  % Tx power 14 dB
noise_sigma=1;
message = "Hello akanksha!" ;
message2 = "Hello AlignTrack!" ;
Fc = 921.5e6; % spectrum center frequency
df=Fc-fc;   % if CFO exists
shift=0;
payload_offset = 10000;
%% CSS modulated Signal
signalIQ1 = [LoRa_Tx_new(message,BW,SF,Power,Fs); zeros(payload_offset + length(LoRa_Tx_new(message2,BW,SF,Power,Fs)) - length(LoRa_Tx_new(message,BW,SF,Power,Fs)), 1)];
signalIQ2 = [zeros(payload_offset, 1); LoRa_Tx_new(message,BW,SF,Power,Fs)];
signalIQ3 = signalIQ1 + signalIQ2;

% spectrogram for one LoRa packet
figure(1)
spectrogram(signalIQ3,1000,0,1000,Fs,'yaxis','centered')
%AWGN noise
noise=noise_sigma*randn(length(signalIQ3),1);
% signal recieved at  LoRa gateway
recieved_signal=signalIQ3+noise;
% recieved signal after windowing and FFT
[recieved_fft] = LoRa_demod_1(recieved_signal,SF,BW,Fs,0);
new_recieved_fft(1,:)=recieved_fft.';
[recieved_fft] = LoRa_demod_1(recieved_signal,SF,BW,Fs,payload_offset);
new_recieved_fft(2,:)=recieved_fft.';
% peak extraction algorith with align track decoding for complete packet
[m n]=size(new_recieved_fft);
for j1=0:1:m-1
    for i1=1:1:n
        [val(j1+1),idx(j1+1)]=max(new_recieved_fft(j1+1,:));
        r(j1+1)=mean(new_recieved_fft(j1+1,:))+6*std(new_recieved_fft(j1+1,:));
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
                if recieved_fft(before+1)>recieved_fft(before) && recieved_fft(before)<recieved_fft(before-1)
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
                if recieved_fft(after)<recieved_fft(after-1) && recieved_fft(after+1)>recieved_fft(after)
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
                new_recieved_fft(j1+1,y+1)=0;
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
                    if k~=j && recieved_fft(AA(k))==recieved_fft(AA(j))
                        issidelobe(k)=1;
                        issidelobe(j)=1;
                    end
                end
            end
        end
        i
    end
    AA=AA(issidelobe~=1);
end