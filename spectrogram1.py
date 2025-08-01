import numpy as np
import scipy.io.wavfile as wav
import math
import cmath

cut_freq = 5000
dspratio = 4
windowsamples = 1024
step = 32
peakwinlen = 4


def fft(input_data):
    N = len(input_data)

    # Convert to complex numbers (imaginary part = 0)
    complex_array = [complex(val, 0) for val in input_data]

    return recursive_fft(complex_array)

def recursive_fft(x):
    N = len(x)
    if N <= 1:
        return x

    even = recursive_fft(x[0::2])
    odd = recursive_fft(x[1::2])

    result = [0] * N
    for k in range(N // 2):
        t = cmath.exp(-2j * math.pi * k / N) * odd[k]
        result[k] = even[k] + t
        result[k + N // 2] = even[k] - t

    return result

def stereo_to_mono(audio, seconds=None):
    if audio.ndim == 2:
        audio = audio.mean(axis=1)

    if seconds is not None:
        max_samples = int(seconds * 44100)
        audio = audio[:max_samples]

    return audio


def low_pass_filter(signal, sample_rate, cutoff_freq):
    rc = 1.0 / (2 * math.pi * cutoff_freq)
    dt = 1.0 / sample_rate
    alpha = dt / (rc + dt)

    filtered_signal = np.zeros_like(signal)
    prev_output = 0.0

    for i in range(len(signal)):
        x = signal[i]
        if i == 0:
            filtered_signal[i] = x * alpha
        else:
            filtered_signal[i] = alpha * x + (1 - alpha) * prev_output
        prev_output = filtered_signal[i]

    return filtered_signal

def downsample(signal, dspratio): #factor 4 means 44.1kHz -> 11.025kHz sample rate
    dsampled = []
    for i in range(len(signal) // dspratio):
        dsampled.append(np.mean(signal[dspratio * i : dspratio * i + dspratio]))
    return dsampled

def slide(postdsp):
    numwindows = math.ceil((len(postdsp) - windowsamples) // step)  # number of windows that fit sliding by step
    magnitudes = np.zeros((numwindows, 513))      # 2D array to store magnitudes per window

    for i in range(numwindows):
        start = i * step
        end = start + windowsamples
        window_data = postdsp[start:end]

        #apply hamming function
        for j in range(windowsamples):
            window_data[j] *= 0.54 - 0.46*math.cos(2*math.pi*j/(windowsamples-1))

        fft_result = fft(window_data)  # your recursive FFT returns complex list
        # convert to magnitude (abs)
        mag = np.array([abs(c) for c in fft_result])

        mag = mag[:513] #mirroring, only first 513 are used

        magnitudes[i, :] = mag


    return magnitudes

import numpy as np
import matplotlib.pyplot as plt

def plot(array2d):
    # Avoid log(0) by adding a small constant
    array2d = np.log1p(array2d)

    plt.figure(figsize=(10, 6))
    plt.imshow(
        array2d.T,
        aspect='auto',
        origin='lower',
        cmap='inferno'
    )
    plt.xlabel('Time Frame')
    plt.ylabel('Frequency Bin')
    plt.title('Spectrogram (log-scaled)')
    plt.colorbar(label='log(1 + Magnitude)')
    plt.show()


def process_wav(input_file, cutoff_freq):
    sample_rate, data = wav.read(input_file)
    data = data.astype(np.float64)
    audiolength = data.shape[0] / sample_rate

    # Step 1: Stereo to Mono
    mono = stereo_to_mono(data, 2)

    # Step 2: Low-pass Filter
    filtered = low_pass_filter(mono, sample_rate, cutoff_freq)

    # Step 3: Downsample by 4x
    downsampled = downsample(filtered, dspratio)

    # Step 4: perform sliding window 1024-sample FFT, step size = 32 and return
    slid = slide(downsampled)

    # Step 5: Plot
    plot(slid)
    

    #optional plot spectrogram
    #plot(slid, 44100, 32, 4)

    # Step 5: choosing peaks

    

# Example usage
process_wav("test.wav", 5000)
