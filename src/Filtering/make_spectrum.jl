################## Check these functions because they might be deprecated #####################################
function fft_spectrum(data::Experiment)
    #FFTW filtering
    t = data.t
    dt = t[2] - t[1]
    freqs = FFTW.fftfreq(length(t), 1.0 / dt) |> fftshift
    over_0 = findall(freqs .> 0)
    n_sweep, n_data, n_ch = size(data)
    fft_data = zeros(Complex, n_sweep, n_data, n_ch)
    for swp in 1:n_sweep
        for ch in 1:n_ch
            fft_data[swp, :, ch] = fft(data.data_array[swp, :, ch]) |> fftshift
        end
    end
    return freqs[over_0], fft_data[:, over_0, :]
end
