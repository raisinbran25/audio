import { fft, ifft } from "./fft.js"

const winSize = 1024;
const overlap = 256; // number of windows
const dspRatio = 3;

var sampleRate;
var maxVal;
var minVal;
let liveMatrix = [];
let zeroBin = [];

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

    let arr = await fileToMono(file);

    // low pass filter

    let low = lowPass(arr);

    // downsample audio

    let dsp = downSample(low);

    let slid = slide(dsp);

    liveMatrix = truncate(slid);

    let heat = magMatrix(slid)

    paintCanvas(heat);

}

function truncate(fftOutput) {
    let newMatrix = [];

    for (let i=0;i<fftOutput.length;i++) {
        let inner = [];

        for (let j=1;j<513;j++) {
            inner.push(fftOutput[i][j]);
        }
        newMatrix.push(inner);

        zeroBin.push(fftOutput[i][0])
    }

    return newMatrix;
}

function magMatrix(fftOutput) {
    let newMatrix = [];

    for (let i=0;i<fftOutput.length;i++) {
        let inner = [];

        for (let j=1;j<513;j++) {
            inner.push(Math.sqrt(Math.pow(fftOutput[i][j].real, 2) + Math.pow(fftOutput[i][j].imag, 2)));
        }
        newMatrix.push(inner);
    }

    return newMatrix;
}

function paintCanvas(heatMatrix) {
    const canvas = document.getElementById("spectrogramCanvas");
    const ctx = canvas.getContext("2d");

    const width = heatMatrix.length;     // number of time frames (columns)
    const height = 511;                   // 

    canvas.width = width;

    const trimmedMatrix = heatMatrix.map(col => col.slice(0, 511));

    // Flatten and find min/max for normalization
    const flat = trimmedMatrix.flat();
    maxVal = flat.reduce((a, b) => Math.max(a, b), -Infinity);
    minVal = flat.reduce((a, b) => Math.min(a, b), Infinity);

    const imgData = ctx.createImageData(width, height);

    for (let x = 0; x < width; x++) {
        for (let y = 0; y < height; y++) {
            // Normalize and log-scale for sensitivity
            const normalized = (trimmedMatrix[x][y] - minVal) / (maxVal - minVal);
            const mag = Math.log10(1 + 9 * normalized); // log compression for better contrast
            const { r, g, b } = heatmapColor(mag);

            // Bottom pixel = bin 0
            const pixelIndex = ((height - 1 - y) * width + x) * 4;

            imgData.data[pixelIndex] = r;
            imgData.data[pixelIndex + 1] = g;
            imgData.data[pixelIndex + 2] = b;
            imgData.data[pixelIndex + 3] = 255;
        }
    }

    ctx.putImageData(imgData, 0, 0);
}

function heatmapColor(value) {
    const v = Math.max(0, Math.min(1, value));
    const biased = Math.pow(v, 0.5); // sensitivity boost

    let r = 0, g = 0, b = 0;
    if (biased < 0.25) {        // Blue → Cyan
        r = 0;
        g = biased * 4 * 255;
        b = 255;
    } else if (biased < 0.5) {  // Cyan → Green
        r = 0;
        g = 255;
        b = (1 - (biased - 0.25) * 4) * 255;
    } else if (biased < 0.75) { // Green → Yellow
        r = (biased - 0.5) * 4 * 255;
        g = 255;
        b = 0;
    } else {                    // Yellow → Red
        r = 255;
        g = (1 - (biased - 0.75) * 4) * 255;
        b = 0;
    }

    return { r: Math.round(r), g: Math.round(g), b: Math.round(b) };
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

function lowPass(signal) {
    const cutFreq = (sampleRate * 0.9) / (2 * dspRatio); // cutoff frequency
    const rc = 1.0 / (2 * Math.PI * cutFreq);
    const dt = 1.0 / sampleRate;
    const alpha = dt / (rc + dt);

    const filteredSignal = new Array(signal.length).fill(0);
    let prevOutput = 0.0;

    for (let i = 0; i < signal.length; i++) {
        const x = signal[i];
        if (i === 0) {
            filteredSignal[i] = x * alpha;
        } else {
            filteredSignal[i] = alpha * x + (1 - alpha) * prevOutput;
        }
        prevOutput = filteredSignal[i];
    }

    return filteredSignal;
}

function downSample(signal) {
    const dsampled = [];
    for (let i = 0; i < Math.floor(signal.length / dspRatio); i++) {
        dsampled.push(signal[dspRatio * i]); // pick sample instead of averaging
    }
    return dsampled;
}

function slide(signal) {
    const hopSize = overlap;
    const numFrames = 1 + Math.floor((signal.length - winSize) / hopSize);

    const hammingWindow = new Array(winSize).fill(0).map((_, n) =>
        0.54 - 0.46 * Math.cos((2 * Math.PI * n) / (winSize - 1))
    );

    const stftMatrix = [];

    for (let i = 0; i < numFrames; i++) {
        const start = i * hopSize;
        const frame = signal.slice(start, start + winSize);

        const windowedFrame = frame.map((val, idx) => val * hammingWindow[idx]);

        const frameFFT = fft(windowedFrame);
        stftMatrix.push(frameFFT);
    }

    return stftMatrix;
}

function invSlide(stftMatrix) {
    const hopSize = overlap;
    const numFrames = stftMatrix.length;
    const signalLength = (numFrames - 1) * hopSize + winSize;

    const hammingWindow = new Array(winSize).fill(0).map((_, n) =>
        0.54 - 0.46 * Math.cos((2 * Math.PI * n) / (winSize - 1))
    );

    const output = new Array(signalLength).fill(0);
    const windowSum = new Array(signalLength).fill(0);

    for (let i = 0; i < numFrames; i++) {
        const start = i * hopSize;

        // Inverse FFT → real frame
        const frameIFFT = ifft(stftMatrix[i]).map(c => c.real);

        // Apply window
        const windowedFrame = frameIFFT.map((val, idx) => val * hammingWindow[idx]);

        // Overlap-add
        for (let j = 0; j < winSize; j++) {
            output[start + j] += windowedFrame[j];
            windowSum[start + j] += hammingWindow[j] ** 2;
        }
    }

    // Normalize
    for (let i = 0; i < signalLength; i++) {
        if (windowSum[i] > 1e-8) {
            output[i] /= windowSum[i];
        }
    }

    return output;
}

function playSignal(signal) {
    const audioCtx = new window.AudioContext();
    const floatSignal = new Float32Array(signal);
    const newRate = sampleRate / dspRatio;

    // Use new sample rate
    const buffer = audioCtx.createBuffer(2, floatSignal.length, newRate);
    buffer.copyToChannel(floatSignal, 0);
    buffer.copyToChannel(floatSignal, 1);


    const source = audioCtx.createBufferSource();
    source.buffer = buffer;
    source.connect(audioCtx.destination);
    source.start();
}







const brushSize = 2;
const brushStrength = 0.9;

const canvas = document.getElementById("spectrogramCanvas");

canvas.addEventListener("mousedown", (e) => {
    paintBrush(e);
    canvas.addEventListener("mousemove", paintBrush);
});

canvas.addEventListener("mouseup", () => {
    canvas.removeEventListener("mousemove", paintBrush);
});

function paintBrush(e) {
    const rect = canvas.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;

    const frameX = Math.floor(mouseX); // time frame index
    const binY = 511 - Math.floor(mouseY); // freq bin (invert Y-axis)

    for (let dx = -brushSize; dx <= brushSize; dx++) {
        for (let dy = -brushSize; dy <= brushSize; dy++) {
            const x = frameX + dx;
            const y = binY + dy;

            if (x >= 0 && x < liveMatrix.length && y >= 0 && y < liveMatrix[0].length) {
                const distance = Math.sqrt(dx * dx + dy * dy);

                if (distance <= brushSize) {
                    // Calculate falloff multiplier
                    const falloff = 1 - (distance / brushSize); 
                    const multiplier = brushStrength * falloff; 

                    // Update real & imaginary parts in liveMatrix
                    liveMatrix[x][y].real *= multiplier;
                    liveMatrix[x][y].imag *= multiplier;

                    // Update pixel color based on new magnitude
                    const canvas = document.getElementById("spectrogramCanvas");
                    const ctx = canvas.getContext("2d");

                    function getMag(num1, num2) {
                        return Math.sqrt(Math.pow(num1, 2) + Math.pow(num2, 2))
                    }

                    // Normalize and log-scale for sensitivity
                    const normalized = ((getMag(liveMatrix[x][y].real, liveMatrix[x][y].imag)) - minVal) / (maxVal - minVal);
                    const mag = Math.log10(1 + 9 * normalized); // log compression for better contrast
                    const { r, g, b } = heatmapColor(mag);

                    setPixel(ctx, x, (511-y), r, g, b);

                }
            }
        }
    }
}








const playButton = document.getElementById("play");
playButton.addEventListener("click", (event) => {

    let unfolded = unfold(liveMatrix);
    console.log(unfolded)

    let inv = invSlide(unfolded);
    console.log("inv slide")

    playSignal(inv);
    console.log("played")

    console.log("done")
})

function unfold(truncatedMatrix) {
    let outMatrix = [];
    const length = truncatedMatrix.length;
    for (let i=0;i<length;i++) {
        let inner = [];
        inner.push(zeroBin[i]);
        for (let j=0;j<512;j++) {
            inner.push(truncatedMatrix[i][j]);
        }
        for (let j=510;j>=0;j--) {
            inner.push(truncatedMatrix[i][j])
            inner[inner.length - 1].real *= -1;
            inner[inner.length - 1].imag *= -1;
        }
        outMatrix.push(inner);
    }
    return outMatrix;
}

function setPixel(ctx, x, y, r, g, b) {
    const imgData = ctx.getImageData(x, y, 1, 1); // read that pixel
    imgData.data[0] = r;
    imgData.data[1] = g;
    imgData.data[2] = b;
    imgData.data[3] = 255;
    ctx.putImageData(imgData, x, y); // write it back
}
