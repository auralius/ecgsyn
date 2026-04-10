# ecgsyn

<p align="center">
    <img width="600"alt="demo" src="https://github.com/user-attachments/assets/b47e1336-6a0a-4224-be92-bfcab892bdc7" />
    <img width="600" alt="fig2" src="https://github.com/user-attachments/assets/08868e80-886f-4416-94fc-b389f0f5759c" />
</p>

---

This repository provides implementations of ECG signal generation algorithms based on nonlinear dynamical models. The focus is on embedded deployment, enabling real-time ECG synthesis on microcontrollers.

Two numerical integration methods are implemented:

- Explicit RK4 (Runge–Kutta 4th Order)
- Implicit Tustin (Bilinear Transform).

Implementation hardware:

- MCU:  [WeActStudio STM32H523](https://github.com/WeActStudio/WeActStudio.STM32H523CoreBoard/tree/master)
    - Peripherals: 12-bit true DAC
        - PA4 as reference DAC
        - PA5 as signal DAC
- OLED SSD1306 0.96"
- Mini oscilloscope


---

Current implementations:

- A Dynamical Model for Generating Synthetic Electrocardiogram Signals by McSharry et al. 

- Generation of ECG signals from a reaction-diffusion model spatially discretized by Quiroz-Juárez et al.

---

Implementations start with MATLAB / Octave, here is an example form the reaction-diffusion ECG model:
  
__Explicit RK4__

<img width="650" alt="image" src="https://github.com/user-attachments/assets/4ff2cd60-d0e0-4356-a5e2-f61690b9847b" />
  
  
__Implicit Tustin__

<img width="650" alt="image" src="https://github.com/user-attachments/assets/dc011890-06a7-4f2d-9079-65848f88b782" />




