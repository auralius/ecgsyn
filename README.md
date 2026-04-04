# ecgsyn

This repository explores the reconstruction of ECG signal generation algorithms based on dynamical 
systems, targeting deployment on embedded platforms.

Hardware:

- MCU: STM32U585
- Peripherals: 12-bit ADC and 12-bit true DAC

Signal loopback setup:

```
PA5 ──[1kΩ]──> A0
```

---

![demo1](https://github.com/user-attachments/assets/409bd2aa-b182-416f-8240-856a8f459e76)

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

