# ecgsyn

This repository contains some attempts to reconstruct ECG signal generation algorithms that are based 
on dynamical equations on an embedded system.

Hardware: STM32U585 with 12-bit ADC and 12-bit True CAD

```
PA5 ──[1kΩ]──> A0
```

<img width="578" height="260" alt="system" src="https://github.com/user-attachments/assets/b1590ddd-1762-41bd-aea3-be34c6403f60" />


![demo1](https://github.com/user-attachments/assets/409bd2aa-b182-416f-8240-856a8f459e76)

Current implementation:

- A Dynamical Model for Generating Synthetic Electrocardiogram Signals by McSharry et.al.  