# Real-time SDR for mono/stereo FM and RDS
Software-Defined Radio developed using C++ and running based on Linux in Respberry Pi

## Software Specification
- Linux
  - Aplay
  - C++ compiler

## Hardware specification (Optional)
- Respberry Pi 4 with Linux
- [NESDR Smart SDRs](https://www.nooelec.com/store/sdr/sdr-receivers/smart.html)


## Quick Start
Only with Linux
1. Clone the code into ResPi
2. Compile using `src/Makefile`
3. Run with command `cat [.raw file] | [compiled executable]  [mode] | aplay -c 2 -f S16_LE -r [audio sampling frequency in Hz]`
4. Sample .raw files could be found in Release

ResPi with NESDR Smart SDRs
1. Clone the code into ResPi
2. Compile using `src/Makefile`
3. Run with command `rtl_sdr -f 107.1M -s [radio sampling freq in Hz] | [compiled executable] [mode] | aplay -c 2 -f S16_LE -r [audio sampling freq in Hz]`

Mode Option

Different mode support different sampling rate.
|                | Mode 0 | Mode 1 | Mode 2 | Mode 3 |
|----------------|--------|--------|--------|--------|
| Radio Fs (kHz) | 2400   | 1440   | 2400   | 1920   |
| AUdio Fs (KHz) | 48     | 48     | 44.1   | 44.1   |
|                |        |        |        |        |


## Feature
- Signal Processing
    - Radio SIgnal was processed with software-defined filters for audio recovery.
- Support Different Channel with Different Sampling Frequency 
    - Supported Channel
      - Mono: main audio
      - Stereo: for seperate left and right channel. Test with `stereo_l0_r9.raw` and `stereo_l1_r8.raw` in sample raw files
      - RDS: radio info
    - Option of sampling freq
      - Use can choose to use more stable option or higher quality option
      - see mode option above
- Runtime Optimized
    - Using multithread to take advantage of multicore
    - Algorithm optimized mathmatically include but not only with
      - Filters and resamples combines and intruction deduction
      - Reduce execution on unwanted data due to sampling frequency change

## Useage
### C++ Library
- cmath
- chrono
- algorithm
- atomic
- thread

### Folder Structure
```
SDR/
├─ build/
├─ data/
├─ doc/
├─ include/
├─ model/
├─ src/
├─ test/
├─ .gitignore
├─ CMakeLists.txt
├─ README.md
```
- data: contain sample .raw files
- include: .h files
- model: class or object model
- src: source code files

## About
This is final course prject of McMaster University COMPENG 3DY4.
### **Do Not** copy any code from this repository.
### Please follow the Mcmaster Academic Integrity Policy and the Code of Conduct of the Professional Engineers of Ontario.
Copying my code will have you results in academic dishonest. 

Partial code credits to McMaster Faculty of Engineering. 
Partial Code Credit to Yilong Wang, Chuchu Zeng, Hairong Bao.