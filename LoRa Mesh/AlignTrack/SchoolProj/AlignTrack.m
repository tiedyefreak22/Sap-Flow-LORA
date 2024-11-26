function AA = AlignTrack(received_fft)
    new_received_fft = received_fft.';
    % peak extraction algorithm with AlignTrack decoding for complete packet
    [row_fft, col_fft] = size(new_received_fft);
    k = 6;
    for fft_row_idx = 0:1:row_fft - 1
        for fft_col_idx = 1:1:col_fft
            [val(fft_row_idx + 1), idx(fft_row_idx + 1)] = max(new_received_fft(fft_row_idx + 1, :));
            r(fft_row_idx + 1) = mean(new_received_fft(fft_row_idx + 1, :)) + k * std(new_received_fft(fft_row_idx + 1, :));
            if val(fft_row_idx + 1) >= r(fft_row_idx + 1) % test against dynamic peak extraction threshold
                I(fft_row_idx + 1, fft_col_idx) = idx(fft_row_idx + 1);
                before = (idx(fft_row_idx + 1) - 1);
                after = (idx(fft_row_idx + 1) + 1);
                % finding local min before index value
                ii1 = 1;
                while ii1 < col_fft
                    if before == 0
                        before = idx(fft_row_idx + 1);
                    end
                    if before-1 <= 0
                        before = idx(fft_row_idx + 1) + 1;
                    end
                    if received_fft(before + 1) > received_fft(before) && received_fft(before) < received_fft(before - 1)
                        b = before;
                        before = before-1;
                        break;
                    else
                        before = before - 1;
                    end
                    ii1 = ii1 + 1;
                end
                % finding local min after index value
                ii2 = 1;
                while ii2 < col_fft
                    if after >= col_fft
                        after = idx(fft_row_idx + 1);
                    end
                    if after + 1 >= col_fft
                        after = idx(fft_row_idx + 1) - 1;
                    end
                    if received_fft(after) < received_fft(after - 1) && received_fft(after + 1) > received_fft(after)
                        a = after;
                        after = after + 1;
                        break;
                    else
                        after = after + 1;
                    end
                    ii2 = ii2 + 1;
                end
                %  remove local min points before and after index value
                for y = before:after
                    new_received_fft(fft_row_idx + 1, y) = 0;
                end
            else
                break;
            end
        end
    end
    [m1, n1] = size(I);
    % sort the index value for further processing
    I = sort(I, 2);
    % side lobe elimination
    for ee = 1:1:m1
        for kk = 1:1:length(I(ee, :))
            if I(ee, kk) == 0
                continue;
            else
                d1 = kk;
                break;
            end
        end
        A.(sprintf('RandomVariable_%d', ee)) = I(ee, d1:end);
        AA = A.(sprintf('RandomVariable_%d', ee));
        issidelobe = zeros(1, length(AA));
        for i = 1:1:length(AA)
            for j = i + 1:1:length(AA)
                for k = 1:1:length(AA)
                    if AA(j) - AA(i) == AA(i) - AA(k)
                        if k ~= j && received_fft(AA(k)) == received_fft(AA(j))
                            issidelobe(k) = 1;
                            issidelobe(j) = 1;
                        end
                    end
                end
            end
        end
        AA = AA(issidelobe ~= 1);
    end
end