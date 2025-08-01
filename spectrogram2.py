import numpy as np
from scipy.io import wavfile

sample_rate = 44100
win_size = 2048
overlap = 50
max_mag = 0


def read_wav(input_file):
    # Read the WAV file
    sr, data = wavfile.read(input_file)
    sample_rate = sr
    return data


def mono(audio):

    global max_mag

    # Ensure input is int16
    if audio.dtype != np.int16:
        raise ValueError("Input must be int16")

    # Compute average safely in int32 to prevent overflow
    avg = ((audio[:, 0].astype(np.int32) + audio[:, 1].astype(np.int32)) // 2)

    # Note max magnitude
    max_mag = np.max(np.abs(avg))

    # Convert to float in range [-1, 1]
    mono_signal = avg.astype(np.float32) / 32768.0  # 32768 is max abs value for int16

    return mono_signal  # 1D float32 array


def slide(signal):
    """
    Performs Hamming windowing and STFT on a 1D signal.

    Parameters:
        signal (np.ndarray): Input 1D audio signal.

    Returns:
        np.ndarray: 2D complex array of STFT results (frames x frequency bins)
    """

    hop_size = int(win_size * (1 - overlap / 100))  # step size between windows
    window = np.hamming(win_size)                  # create Hamming window

    # Number of frames
    num_frames = 1 + (len(signal) - win_size) // hop_size

    # STFT matrix
    stft_matrix = np.zeros((num_frames, win_size), dtype=np.complex64)

    for i in range(num_frames):
        start = i * hop_size
        frame = signal[start:start + win_size]

        # Apply window
        windowed_frame = frame * window

        # FFT
        stft_matrix[i, :] = np.fft.fft(windowed_frame)

    return stft_matrix


def inv_slide(stft_matrix):
    """
    Performs inverse STFT using overlap-add, compensating for the Hamming window.

    Parameters:
        stft_matrix (np.ndarray): 2D STFT matrix (frames x frequency bins)

    Returns:
        np.ndarray: Reconstructed 1D signal
    """
    hop_size = int(win_size * (1 - overlap / 100))  # step size
    window = np.hamming(win_size)                  # Hamming window

    num_frames = stft_matrix.shape[0]
    signal_length = (num_frames - 1) * hop_size + win_size

    # Output buffers
    output = np.zeros(signal_length)
    window_sum = np.zeros(signal_length)  # for normalization

    for i in range(num_frames):
        start = i * hop_size

        # Inverse FFT to get windowed frame
        frame = np.fft.ifft(stft_matrix[i]).real

        # Apply window again
        windowed_frame = frame * window

        # Overlap-add
        output[start:start + win_size] += windowed_frame
        window_sum[start:start + win_size] += window**2  # track overlap energy

    # Normalize to account for window overlap
    nonzero_indices = window_sum > 1e-8
    output[nonzero_indices] /= window_sum[nonzero_indices]

    return output


def inv_mono(mono_signal):
    """
    Convert a normalized mono float signal back to a stereo int16 signal.

    Parameters:
        mono_signal (np.ndarray): 1D float signal in range [-1, 1]
        max_mag (int): Maximum magnitude (from original mono() function)

    Returns:
        np.ndarray: Stereo int16 array of shape (N, 2)
    """
    if not np.issubdtype(mono_signal.dtype, np.floating):
        raise ValueError("Input must be a float array")

    # Clip to [-1, 1] to avoid overflow
    mono_signal = np.clip(mono_signal, -1.0, 1.0)

    # Scale to match the original max magnitude (in int32)
    scaled_signal = mono_signal * max_mag

    # Convert to int16
    int_signal = scaled_signal.astype(np.int16)

    # Duplicate to stereo (N, 2) where L = R
    stereo_signal = np.column_stack((int_signal, int_signal))

    return stereo_signal


def export_wav(signal, filename):
    wavfile.write(filename, sample_rate, signal)



# flow

a = read_wav("les1sec.wav")
b = mono(a)
c = slide(b)
d = inv_slide(c)
e = inv_mono(d)
export_wav(e, "test.wav")