const windowSize = 2048;
const overlap = 50; // percentage value

var sampleRate;

const fileInput = document.getElementById("fileInput");

fileInput.addEventListener("change", (event) => {

    const file = event.target.files[0]; // get file

    if (file.name.split('.').pop() != "mp3" && file.name.split('.').pop() != "wav") {

        console.log("invalid file type");
        return;

    }
    else {

        console.log("valid file type!")
        flow(file)

    };
});

async function flow(file) {

    arr = await fileToMono(file);

    console.log(slide(arr))

    // playSignal(arr);

    console.log("Done");

}

async function fileToMono(file) {

    const audioContext = new AudioContext();
    await audioContext.resume();

    // Read file into ArrayBuffer
    const arrayBuffer = await file.arrayBuffer();

    // Decode MP3/WAV into AudioBuffer
    const audioBuffer = await audioContext.decodeAudioData(arrayBuffer);

    const numChannels = audioBuffer.numberOfChannels;
    const length = audioBuffer.length;
    sampleRate = audioBuffer.sampleRate;

    console.log("Channels:", numChannels);
    console.log("Length:", length);
    console.log("Sample rate:", sampleRate);

    // Extract each channel into an array
    const signal = [];
    for (let ch = 0; ch < numChannels; ch++) {
        signal.push(audioBuffer.getChannelData(ch)); // Float32Array [-1, 1]
    }

    // average signal to mono
    const mono = [];
    for (let i = 0; i < length; i++) {
        mono.push((signal[0][i] + signal[1][i]) / 2);
    };

    return mono; // return the mono signal
}

import FFT from 'fft.js'; // npm install fft.js

function slide(signal) {
    /**
     * Performs Hamming windowing and STFT on a 1D signal.
     *
     * @param {Float32Array|Array} signal - Input 1D audio signal.
     * @param {number} windowSize - Window size (number of samples per frame).
     * @param {number} overlap - Overlap percentage (0–100).
     * @returns {Array<Array<{real:number, imag:number}>>} 2D STFT result
     */

    const hopSize = Math.floor(windowSize * (1 - overlap / 100)); // step size
    const numFrames = 1 + Math.floor((signal.length - windowSize) / hopSize);

    // Hamming window
    const window = new Array(windowSize).fill(0).map((_, n) =>
        0.54 - 0.46 * Math.cos((2 * Math.PI * n) / (windowSize - 1))
    );

    const stftMatrix = [];

    for (let i = 0; i < numFrames; i++) {
        const start = i * hopSize;
        const frame = signal.slice(start, start + windowSize);

        // Apply window
        const windowedFrame = frame.map((val, idx) => val * window[idx]);

        // FFT
        const fft = new FFT(windowSize);
        const input = fft.createComplexArray();
        const output = fft.createComplexArray();

        // Real input, imag part 0
        for (let j = 0; j < windowSize; j++) {
            input[2 * j] = windowedFrame[j];
            input[2 * j + 1] = 0;
        }

        fft.transform(output, input);

        // Store FFT result as array of complex numbers
        const frameFFT = [];
        for (let j = 0; j < windowSize; j++) {
            frameFFT.push({
                real: output[2 * j],
                imag: output[2 * j + 1],
            });
        }

        stftMatrix.push(frameFFT);
    }

    return stftMatrix;
}


function playSignal(signal) {
    const audioCtx = new window.AudioContext();
    floatSignal = new Float32Array(signal)

    // Use global sample_rate variable
    const buffer = audioCtx.createBuffer(2, floatSignal.length, sampleRate);
    buffer.copyToChannel(floatSignal, 0);
    buffer.copyToChannel(floatSignal, 1);


    const source = audioCtx.createBufferSource();
    source.buffer = buffer;
    source.connect(audioCtx.destination);
    source.start();
}



