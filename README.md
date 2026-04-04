# ecgsyn

This repository explores the reconstruction of ECG signal generation algorithms based on dynamical 
systems, targeting deployment on embedded platforms.

Hardware:

- MCU:  [WeActStudio STM32U585CI](https://github.com/WeActStudio/WeActStudio.STM32U585Cx_CoreBoard)
- Peripherals: 12-bit ADC and 12-bit true DAC

Signal loopback setup:

```
PA5 ──[1kΩ]──> A0
```

<img width="578" height="260" alt="system" src="https://github.com/user-attachments/assets/9b1cf345-e2cf-44f1-9c8a-e61b376b41ad" />

---

Current implementations:

- A Dynamical Model for Generating Synthetic Electrocardiogram Signals by McSharry et al. 

The implementation was set to be as faithful as possible with the original paper. The original 
reference implementation (`do_run(...)`) is designed for desktop environments and assumes abundant 
memory.

To make the algorithm suitable for MCU deployment, the core routine has been adapted into:
 `build_block_mv(...)`  

These modifications were developed through an iterative process of analysis and refinement, with 
assistance from ChatGPT for identifying structural issues in the original algorithm (particularly 
memory scaling and RR expansion) and proposing MCU-safe adaptations.

![demo](https://github.com/user-attachments/assets/75a439af-6e92-4347-a88d-b26342810378)

