#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <cassert>

namespace FourierTransforms {

    using Complex = std::complex<double>;
    using ComplexVector = std::vector<Complex>;
    
    const double PI = 3.141592653589793;
    const Complex I(0.0, 1.0); // Imaginary unit
    
    // ===== UTILITY FUNCTIONS =====
    
    // Check if n is a power of 2
    bool isPowerOfTwo(size_t n) {
        return n > 0 && (n & (n - 1)) == 0;
    }
    
    // Find next power of 2 >= n
    size_t nextPowerOfTwo(size_t n) {
        if (n <= 1) return 1;
        size_t result = 1;
        while (result < n) {
            result <<= 1;
        }
        return result;
    }
    
    // Bit-reverse permutation for FFT
    void bitReversePermutation(ComplexVector& data) {
        size_t n = data.size();
        size_t j = 0;
        
        for (size_t i = 1; i < n - 1; ++i) {
            size_t bit = n >> 1;
            while (j & bit) {
                j ^= bit;
                bit >>= 1;
            }
            j ^= bit;
            
            if (i < j) {
                std::swap(data[i], data[j]);
            }
        }
    }
    
    // ===== CORE FFT IMPLEMENTATIONS =====
    
    // Cooley-Tukey FFT algorithm (radix-2, decimation-in-time)
    void fftRadix2(ComplexVector& data, bool inverse = false) {
        size_t n = data.size();
        
        if (n <= 1) return;
        
        if (!isPowerOfTwo(n)) {
            throw std::invalid_argument("FFT size must be a power of 2");
        }
        
        // Bit-reverse permutation
        bitReversePermutation(data);
        
        // FFT computation
        double sign = inverse ? 1.0 : -1.0;
        
        for (size_t length = 2; length <= n; length *= 2) {
            double angle = sign * 2.0 * PI / length;
            Complex wlen(std::cos(angle), std::sin(angle));
            
            for (size_t i = 0; i < n; i += length) {
                Complex w(1.0, 0.0);
                
                for (size_t j = 0; j < length / 2; ++j) {
                    Complex u = data[i + j];
                    Complex v = data[i + j + length / 2] * w;
                    
                    data[i + j] = u + v;
                    data[i + j + length / 2] = u - v;
                    
                    w *= wlen;
                }
            }
        }
        
        // Normalize for inverse transform
        if (inverse) {
            for (auto& x : data) {
                x /= static_cast<double>(n);
            }
        }
    }
    
    // Fast Fourier Transform
    ComplexVector fft(const ComplexVector& input) {
        ComplexVector data = input;
        fftRadix2(data, false);
        return data;
    }
    
    // Inverse Fast Fourier Transform
    ComplexVector ifft(const ComplexVector& input) {
        ComplexVector data = input;
        fftRadix2(data, true);
        return data;
    }
    
    // FFT for real-valued input (returns full complex result)
    ComplexVector fftReal(const std::vector<double>& input) {
        ComplexVector complex_input(input.begin(), input.end());
        return fft(complex_input);
    }
    
    // Real-valued FFT optimized (returns only positive frequencies)
    ComplexVector rfft(const std::vector<double>& input) {
        ComplexVector full_result = fftReal(input);
        size_t n = input.size();
        
        // Return only first (n/2 + 1) elements for real input
        ComplexVector result(full_result.begin(), full_result.begin() + n/2 + 1);
        return result;
    }
    
    // Inverse of real FFT
    std::vector<double> irfft(const ComplexVector& input, size_t output_size) {
        // Reconstruct full complex spectrum (conjugate symmetry)
        ComplexVector full_input(output_size);
        
        // Copy positive frequencies
        size_t half_size = input.size();
        for (size_t i = 0; i < half_size; ++i) {
            full_input[i] = input[i];
        }
        
        // Mirror negative frequencies (conjugate symmetry)
        for (size_t i = 1; i < half_size - 1; ++i) {
            full_input[output_size - i] = std::conj(input[i]);
        }
        
        ComplexVector complex_result = ifft(full_input);
        
        // Extract real part
        std::vector<double> result;
        result.reserve(output_size);
        for (const auto& c : complex_result) {
            result.push_back(c.real());
        }
        
        return result;
    }
    
    // ===== ADVANCED FFT FUNCTIONS =====
    
    // Zero-pad input to next power of 2
    ComplexVector zeroPad(const ComplexVector& input) {
        size_t padded_size = nextPowerOfTwo(input.size());
        ComplexVector padded(padded_size, Complex(0.0, 0.0));
        std::copy(input.begin(), input.end(), padded.begin());
        return padded;
    }
    
    // FFT with automatic zero-padding
    ComplexVector fftAutoPad(const ComplexVector& input) {
        ComplexVector padded = zeroPad(input);
        return fft(padded);
    }
    
    // Convolution using FFT (faster for large sequences)
    ComplexVector convolution(const ComplexVector& a, const ComplexVector& b) {
        size_t result_size = a.size() + b.size() - 1;
        size_t fft_size = nextPowerOfTwo(result_size);
        
        ComplexVector a_padded(fft_size, Complex(0.0, 0.0));
        ComplexVector b_padded(fft_size, Complex(0.0, 0.0));
        
        std::copy(a.begin(), a.end(), a_padded.begin());
        std::copy(b.begin(), b.end(), b_padded.begin());
        
        ComplexVector a_fft = fft(a_padded);
        ComplexVector b_fft = fft(b_padded);
        
        // Element-wise multiplication in frequency domain
        ComplexVector conv_fft(fft_size);
        for (size_t i = 0; i < fft_size; ++i) {
            conv_fft[i] = a_fft[i] * b_fft[i];
        }
        
        ComplexVector result = ifft(conv_fft);
        result.resize(result_size); // Trim to actual result size
        
        return result;
    }
    
    // Cross-correlation using FFT
    ComplexVector crossCorrelation(const ComplexVector& a, const ComplexVector& b) {
        size_t result_size = a.size() + b.size() - 1;
        size_t fft_size = nextPowerOfTwo(result_size);
        
        ComplexVector a_padded(fft_size, Complex(0.0, 0.0));
        ComplexVector b_padded(fft_size, Complex(0.0, 0.0));
        
        std::copy(a.begin(), a.end(), a_padded.begin());
        std::copy(b.begin(), b.end(), b_padded.begin());
        
        ComplexVector a_fft = fft(a_padded);
        ComplexVector b_fft = fft(b_padded);
        
        // Cross-correlation: conj(A) * B
        ComplexVector corr_fft(fft_size);
        for (size_t i = 0; i < fft_size; ++i) {
            corr_fft[i] = std::conj(a_fft[i]) * b_fft[i];
        }
        
        ComplexVector result = ifft(corr_fft);
        result.resize(result_size);
        
        return result;
    }
    
    // ===== SPECTRAL ANALYSIS FUNCTIONS =====
    
    // Compute magnitude spectrum
    std::vector<double> magnitudeSpectrum(const ComplexVector& fft_result) {
        std::vector<double> magnitude;
        magnitude.reserve(fft_result.size());
        
        for (const auto& c : fft_result) {
            magnitude.push_back(std::abs(c));
        }
        
        return magnitude;
    }
    
    // Compute phase spectrum
    std::vector<double> phaseSpectrum(const ComplexVector& fft_result) {
        std::vector<double> phase;
        phase.reserve(fft_result.size());
        
        for (const auto& c : fft_result) {
            phase.push_back(std::arg(c));
        }
        
        return phase;
    }
    
    // Compute power spectrum (magnitude squared)
    std::vector<double> powerSpectrum(const ComplexVector& fft_result) {
        std::vector<double> power;
        power.reserve(fft_result.size());
        
        for (const auto& c : fft_result) {
            double mag = std::abs(c);
            power.push_back(mag * mag);
        }
        
        return power;
    }
    
    // ===== WINDOWING FUNCTIONS =====
    
    // Hann window (good general-purpose window)
    std::vector<double> hannWindow(size_t length) {
        std::vector<double> window(length);
        for (size_t i = 0; i < length; ++i) {
            window[i] = 0.5 * (1.0 - std::cos(2.0 * PI * i / (length - 1)));
        }
        return window;
    }
    
    // Hamming window
    std::vector<double> hammingWindow(size_t length) {
        std::vector<double> window(length);
        for (size_t i = 0; i < length; ++i) {
            window[i] = 0.54 - 0.46 * std::cos(2.0 * PI * i / (length - 1));
        }
        return window;
    }
    
    // Blackman window (low side lobes)
    std::vector<double> blackmanWindow(size_t length) {
        std::vector<double> window(length);
        for (size_t i = 0; i < length; ++i) {
            double n = static_cast<double>(i);
            double N = static_cast<double>(length - 1);
            window[i] = 0.42 - 0.5 * std::cos(2.0 * PI * n / N) + 0.08 * std::cos(4.0 * PI * n / N);
        }
        return window;
    }
    
    // Apply window to signal
    std::vector<double> applyWindow(const std::vector<double>& signal, const std::vector<double>& window) {
        if (signal.size() != window.size()) {
            throw std::invalid_argument("Signal and window must be same size");
        }
        
        std::vector<double> windowed(signal.size());
        for (size_t i = 0; i < signal.size(); ++i) {
            windowed[i] = signal[i] * window[i];
        }
        
        return windowed;
    }
    
    // ===== SIGNAL GENERATION (for testing and demos) =====
    
    // Generate sine wave
    std::vector<double> generateSine(size_t length, double frequency, double sample_rate = 1.0, double amplitude = 1.0, double phase = 0.0) {
        std::vector<double> signal(length);
        for (size_t i = 0; i < length; ++i) {
            double t = static_cast<double>(i) / sample_rate;
            signal[i] = amplitude * std::sin(2.0 * PI * frequency * t + phase);
        }
        return signal;
    }
    
    // Generate cosine wave
    std::vector<double> generateCosine(size_t length, double frequency, double sample_rate = 1.0, double amplitude = 1.0, double phase = 0.0) {
        std::vector<double> signal(length);
        for (size_t i = 0; i < length; ++i) {
            double t = static_cast<double>(i) / sample_rate;
            signal[i] = amplitude * std::cos(2.0 * PI * frequency * t + phase);
        }
        return signal;
    }
    
    // Generate white noise
    std::vector<double> generateWhiteNoise(size_t length, double amplitude = 1.0, unsigned seed = 0) {
        std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
        std::normal_distribution<double> dist(0.0, amplitude);
        
        std::vector<double> signal(length);
        for (size_t i = 0; i < length; ++i) {
            signal[i] = dist(rng);
        }
        return signal;
    }
    
} // namespace FourierTransforms

// ===== DEMONSTRATION AND USAGE =====

int main() {
    using namespace FourierTransforms;
    
    std::cout << "=== FOURIER TRANSFORMS MODULE DEMO ===" << std::endl << std::endl;
    
    // Parameters for demo
    const size_t N = 64;           // Signal length (power of 2 for efficient FFT)
    const double fs = 100.0;       // Sample rate
    const double f1 = 5.0;         // Frequency 1
    const double f2 = 15.0;        // Frequency 2
    
    std::cout << "Signal parameters:" << std::endl;
    std::cout << "  Length: " << N << " samples" << std::endl;
    std::cout << "  Sample rate: " << fs << " Hz" << std::endl;
    std::cout << "  Frequencies: " << f1 << " Hz and " << f2 << " Hz" << std::endl << std::endl;
    
    // === SIGNAL GENERATION ===
    std::cout << "1. SIGNAL GENERATION" << std::endl;
    
    // Generate test signals
    auto sine1 = generateSine(N, f1, fs, 1.0);
    auto sine2 = generateSine(N, f2, fs, 0.5);
    auto noise = generateWhiteNoise(N, 0.1, 42);
    
    // Combine signals
    std::vector<double> signal(N);
    for (size_t i = 0; i < N; ++i) {
        signal[i] = sine1[i] + sine2[i] + noise[i];
    }
    
    std::cout << "Generated composite signal: 1.0*sin(2π*" << f1 << "*t) + 0.5*sin(2π*" << f2 << "*t) + noise" << std::endl;
    std::cout << "First 8 samples: ";
    for (size_t i = 0; i < 8; ++i) {
        std::cout << signal[i] << " ";
    }
    std::cout << std::endl << std::endl;
    
    // === BASIC FFT ===
    std::cout << "2. BASIC FFT ANALYSIS" << std::endl;
    
    // Compute FFT
    auto fft_result = fftReal(signal);
    auto magnitude = magnitudeSpectrum(fft_result);
    auto phase = phaseSpectrum(fft_result);
    
    std::cout << "FFT computed. Result size: " << fft_result.size() << std::endl;
    std::cout << "Peak magnitudes (first 16 bins):" << std::endl;
    for (size_t i = 0; i < std::min(size_t(16), magnitude.size()); ++i) {
        double freq = i * fs / N;
        std::cout << "  Bin " << i << " (" << freq << " Hz): " << magnitude[i] << std::endl;
    }
    std::cout << std::endl;
    
    // === INVERSE FFT VERIFICATION ===
    std::cout << "3. INVERSE FFT VERIFICATION" << std::endl;
    
    auto reconstructed_complex = ifft(fft_result);
    std::vector<double> reconstructed_real;
    for (const auto& c : reconstructed_complex) {
        reconstructed_real.push_back(c.real());
    }
    
    // Calculate reconstruction error
    double max_error = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double error = std::abs(signal[i] - reconstructed_real[i]);
        max_error = std::max(max_error, error);
    }
    
    std::cout << "Maximum reconstruction error: " << max_error << std::endl;
    std::cout << "Reconstruction " << (max_error < 1e-10 ? "SUCCESSFUL" : "FAILED") << std::endl;
    std::cout << std::endl;
    
    // === WINDOWING DEMO ===
    std::cout << "4. WINDOWING FUNCTIONS" << std::endl;
    
    auto hann_win = hannWindow(N);
    auto hamming_win = hammingWindow(N);
    auto blackman_win = blackmanWindow(N);
    
    auto windowed_signal = applyWindow(signal, hann_win);
    auto windowed_fft = fftReal(windowed_signal);
    auto windowed_magnitude = magnitudeSpectrum(windowed_fft);
    
    std::cout << "Applied Hann window to signal" << std::endl;
    std::cout << "Windowed signal peak magnitudes (first 8 bins):" << std::endl;
    for (size_t i = 0; i < 8; ++i) {
        double freq = i * fs / N;
        std::cout << "  " << freq << " Hz: " << windowed_magnitude[i] << std::endl;
    }
    std::cout << std::endl;
    
    // === CONVOLUTION DEMO ===
    std::cout << "5. CONVOLUTION USING FFT" << std::endl;
    
    // Create simple impulse response (low-pass filter)
    std::vector<double> impulse_response = {0.25, 0.5, 0.25}; // Simple 3-point average
    ComplexVector impulse_complex(impulse_response.begin(), impulse_response.end());
    ComplexVector signal_complex(signal.begin(), signal.end());
    
    auto convolved = convolution(signal_complex, impulse_complex);
    
    std::cout << "Convolved signal with 3-point averaging filter" << std::endl;
    std::cout << "Original vs Filtered (first 8 samples):" << std::endl;
    for (size_t i = 0; i < 8 && i < convolved.size(); ++i) {
        std::cout << "  " << signal[i] << " -> " << convolved[i].real() << std::endl;
    }
    std::cout << std::endl;
    
    // === REAL FFT DEMO ===
    std::cout << "6. OPTIMIZED REAL FFT" << std::endl;
    
    auto rfft_result = rfft(signal);
    std::cout << "Real FFT size: " << rfft_result.size() << " (vs full FFT: " << fft_result.size() << ")" << std::endl;
    
    auto reconstructed_from_rfft = irfft(rfft_result, N);
    double rfft_error = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double error = std::abs(signal[i] - reconstructed_from_rfft[i]);
        rfft_error = std::max(rfft_error, error);
    }
    
    std::cout << "Real FFT reconstruction error: " << rfft_error << std::endl;
    std::cout << "Real FFT " << (rfft_error < 1e-10 ? "SUCCESSFUL" : "FAILED") << std::endl;
    std::cout << std::endl;
    
    // === PERFORMANCE NOTE ===
    std::cout << "7. PERFORMANCE NOTES" << std::endl;
    std::cout << "• FFT is most efficient for power-of-2 lengths" << std::endl;
    std::cout << "• Real FFT saves ~50% computation for real-valued signals" << std::endl;
    std::cout << "• Windowing reduces spectral leakage" << std::endl;
    std::cout << "• FFT-based convolution is faster than direct convolution for long sequences" << std::endl;
    std::cout << std::endl;
    
    std::cout << "=== FOURIER TRANSFORMS MODULE COMPLETE ===" << std::endl;
    
    return 0;
}