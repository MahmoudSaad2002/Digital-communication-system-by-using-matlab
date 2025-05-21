# Performance Analysis of Digital Modulation Techniques (BPSK, BFSK, and QASK) under AWGN Channels

## 📘 Overview
This project presents the MATLAB-based simulation and performance analysis of three key digital modulation techniques — **BPSK (Binary Phase Shift Keying)**, **BFSK (Binary Frequency Shift Keying)**, and **QASK (Quadrature Amplitude Shift Keying / 4-ASK)** — under Additive White Gaussian Noise (AWGN) channels. It compares the theoretical and simulated **Bit Error Rate (BER)** for each scheme at different SNR levels.

---

## 🎯 Objective
The primary goal is to evaluate how these digital modulation schemes perform under noisy conditions and verify how well the simulated Bit Error Rate (BER) matches theoretical expectations. This serves as a foundational study in digital communication systems.

---

## 🛠️ Tools and Technologies
- **MATLAB R2023a**
- Simulation scripts for each modulation type
- Plotting tools for BER, spectrum, and constellation diagrams

---

## 📂 Project Structure

---

## ⚙️ Modulation Techniques

### 1. BPSK (Binary Phase Shift Keying)
- Maps `0 → -1` and `1 → +1`
- Simple phase modulation with carrier
- High noise immunity
- Simulated and theoretical BER are closely matched

### 2. BFSK (Binary Frequency Shift Keying)
- Maps `0 → f₀`, `1 → f₁` using two separate frequencies
- Robust frequency modulation
- Receiver detects frequency shifts to determine bits
- Shows reliable performance across SNR values

### 3. QASK (Quadrature ASK or 4-ASK)
- 2-bit symbols represented using 4 different amplitudes
- Efficient bandwidth utilization
- Requires higher SNR for accurate detection
- Simulated BER confirms expected trade-offs between efficiency and noise sensitivity

---

## 📊 Key Results

- **BER vs SNR Plots**: Each modulation scheme includes BER curves showing how simulated results compare with theoretical values.
- **Constellation Diagrams**: Visual validation of symbol spacing and clustering under noise.
- **Spectrum Analysis**: Frequency domain visualization for each modulation type.
- **Bandwidth Estimation**: Calculated using `obw()` MATLAB function.

---

## 📌 How to Run
1. Open the `.m` files (e.g., `BPSK_simulation.m`) in MATLAB.
2. Run the script to:
   - Generate random data
   - Modulate the signal
   - Add AWGN noise
   - Demodulate and calculate BER
   - Plot relevant figures

Each script is self-contained and produces:
- Signal plots (time domain)
- Constellation diagrams
- Spectral analysis
- BER comparison graph

---

## 📈 Sample Output
- **BPSK**: BER plot confirms excellent error performance at low SNR
- **BFSK**: Demonstrates robustness using distinct frequency mapping
- **QASK**: Provides greater data efficiency but requires better SNR for reliability

---

## 📚 Conclusion
- **BPSK**: Best performance under noisy conditions; simple and reliable.
- **BFSK**: Suitable where frequency diversity is advantageous.
- **QASK (4-ASK)**: Efficient for high data rate systems but sensitive to noise.

Each modulation method is suited for different real-world applications depending on system requirements like **bandwidth**, **robustness**, and **complexity**.

---

## 👤 Author
**Mahmoud Saad Mostafa**  
Benha University  
Electrical Engineering, Communications and Computers Engineering

Supervisor: **Dr. Ashraf Yahia Hassan**

---

## 🏁 License
This project is developed for academic use. You may use or adapt it for educational or research purposes.
