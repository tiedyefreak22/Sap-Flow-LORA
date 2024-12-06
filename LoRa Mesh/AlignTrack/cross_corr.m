function [shift, chirp_len, k] = cross_corr(SF, BW, Fs, cfo, received_signal, phy)
    % Cross-Correlation for Preamble Detection
    ref_upchirp = phy.chirp(true, SF, BW, Fs, 0, cfo, 0);
    chirp_len = length(ref_upchirp);
    k = BW / (chirp_len / Fs); % Chirp rate
    correlation = abs(xcorr(received_signal, ref_upchirp));
    [pks, locs] = findpeaks(correlation, 'MinPeakHeight', 0.5e4);
    shift = locs(1) - length(received_signal);
end