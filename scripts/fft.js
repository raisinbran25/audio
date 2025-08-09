export function fft(inputData) {
    // Convert input array to complex (as {real, imag})
    const complexArray = inputData.map(val => ({ real: val, imag: 0 }));
    return recursiveFFT(complexArray);
}

export function recursiveFFT(x) {
    const N = x.length;
    if (N <= 1) return x;

    // Split even and odd
    const even = recursiveFFT(x.filter((_, i) => i % 2 === 0));
    const odd = recursiveFFT(x.filter((_, i) => i % 2 === 1));

    const result = new Array(N);
    for (let k = 0; k < N / 2; k++) {
        const angle = (-2 * Math.PI * k) / N;

        // twiddle factor e^(-j2πk/N)
        const twiddle = {
            real: Math.cos(angle),
            imag: Math.sin(angle)
        };

        // t = twiddle * odd[k]
        const t = {
            real: twiddle.real * odd[k].real - twiddle.imag * odd[k].imag,
            imag: twiddle.real * odd[k].imag + twiddle.imag * odd[k].real
        };

        // result[k] = even[k] + t
        result[k] = {
            real: even[k].real + t.real,
            imag: even[k].imag + t.imag
        };

        // result[k + N/2] = even[k] - t
        result[k + N / 2] = {
            real: even[k].real - t.real,
            imag: even[k].imag - t.imag
        };
    }

    return result;
}

export function ifft(inputData) {
    const N = inputData.length;

    // Conjugate
    const conjugated = inputData.map(c => ({ real: c.real, imag: -c.imag }));

    // Forward FFT on conjugated input
    const fftResult = recursiveFFT(conjugated);

    // Conjugate again and scale
    return fftResult.map(c => ({
        real: c.real / N,
        imag: -c.imag / N
    }));
}