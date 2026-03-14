function received_signal = send_receive_signal(offset_time, symbols1, symbols2, symbols3, Fs, Ts, phy, noise_sigma)
    while true
        try
            payload_offset1 = round(Fs * (randi([round(offset_time / 2), offset_time]) / 1000)); % samples
            payload_offset2 = round(Fs * (randi([round(offset_time / 2), offset_time]) / 1000)); % offset from second signal, not first
    
            switch find([numel(symbols1), numel(symbols2), numel(symbols3)] == max(max([numel(symbols1), numel(symbols2), numel(symbols3)])), 1, 'first')
                case 1
                    diff = (max(max([numel(symbols1), numel(symbols2), numel(symbols3)])) - min(min([numel(symbols1), numel(symbols2), numel(symbols3)]))) * Ts * Fs;
                    signalIQ1 = zeros(numel(phy.modulate(symbols1)) + payload_offset1 + payload_offset2 - diff, 1);
                    signalIQ2 = signalIQ1;
                    signalIQ3 = signalIQ1;
                case 2
                    diff = (max(max([numel(symbols1), numel(symbols2), numel(symbols3)])) - min(min([numel(symbols1), numel(symbols2), numel(symbols3)]))) * Ts * Fs;
                    signalIQ2 = zeros(numel(phy.modulate(symbols2)) + payload_offset1 + payload_offset2 - diff, 1);
                    signalIQ1 = signalIQ2;
                    signalIQ3 = signalIQ2;
                case 3
                    diff = (max(max([numel(symbols1), numel(symbols2), numel(symbols3)])) - min(min([numel(symbols1), numel(symbols2), numel(symbols3)]))) * Ts * Fs;
                    signalIQ3 = zeros(numel(phy.modulate(symbols3)) + payload_offset1 + payload_offset2 - diff, 1);
                    signalIQ1 = signalIQ3;
                    signalIQ2 = signalIQ3;
                otherwise
                    error("Can't determine max signal length");
            end
    
            signalIQ1(1:numel(phy.modulate(symbols1))) = phy.modulate(symbols1);
            signalIQ2(payload_offset1 + 1:numel(phy.modulate(symbols2)) + payload_offset1) = phy.modulate(symbols2);
            signalIQ3(payload_offset1 + payload_offset2 + 1:numel(phy.modulate(symbols3)) + payload_offset1 + payload_offset2) = phy.modulate(symbols3);
            signalIQtotal = signalIQ1 + signalIQ2 + signalIQ3;
            break
        end
    end
    
    % Add zero leader
    initial_offset = Fs * (randi([1, 10]) / 1000);
    signalIQtotal = [zeros(initial_offset, 1); signalIQtotal];
    
    % AWGN noise
    noise = noise_sigma * randn(length(signalIQtotal), 1);
    
    % signal recieved at LoRa gateway
    received_signal = signalIQtotal + noise;
end