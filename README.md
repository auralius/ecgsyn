# ecgsyn

This repository explores the reconstruction of ECG signal generation algorithms based on dynamical 
systems, targeting deployment on embedded platforms.

Hardware:

- MCU:  [WeActStudio STM32H523](https://github.com/WeActStudio/WeActStudio.STM32U585Cx_CoreBoard)
    - Peripherals: 12-bit true DAC
        - PA4 as reference DAC
        - PA5 as signal DAC
- OLED SSD1306 0.96"
- Mini oscilloscope


---

Current implementations:

- A Dynamical Model for Generating Synthetic Electrocardiogram Signals by McSharry et al. 

The implementation is set to be as faithful as possible with the original paper. The original 
reference implementation (`do_run(...)`) is designed for desktop environments and assumes abundant 
memory.

To make the algorithm suitable for MCU deployment, the core routine has been adapted into:
 `build_block_mv(...)`  

These modifications were developed through an iterative process of analysis and refinement, with 
assistance from ChatGPT for identifying structural issues in the original algorithm (particularly 
memory scaling and RR expansion) and proposing MCU-safe adaptations.

__Besides using the explicit Runge-Kutta solver as in the original paper, we also propose and implement an implicit bilinear (Tustin) solver.__

![demo](https://github.com/user-attachments/assets/75a439af-6e92-4347-a88d-b26342810378)

---

- Generation of ECG signals from a reaction-diffusion model spatially discretized by Quiroz-Juárez et al.

The paper does not provide accompanying software; however, the implementation is relatively straightforward. The nonlinear systems introduced—particularly in the quasiperiodic and chaotic regimes—are inherently sensitive to numerical methods. As a result, different solvers can produce noticeably different waveforms. In addition, the paper does not specify the value of `HRbpm` used in the simulations, which further contributes to variability in reproduction.

__Similar to McSharry ECG model, this ECG model also uses explicit Runge-Kutta solver. Thus, we also propose and implement an implicit bilinear (Tustin) solver.__


```octave
 % Define the dynamic system in (x1, x2, x3, x4)
 % dx/dt as a column vector

 % Factor Gamma_t
 Gamma_t = 0.08804 * HR_bpm  - 0.06754;

 ode_fun = @(t, x) Gamma_t .* [ ...
     x(1) - x(2) - C*x(1)*x(2) - x(1)*x(2)^2; ...                                  % dx1/dt
     H*x(1) - 3*x(2) + C*x(1)*x(2) + x(1)*x(2)^2 + beta*(x(4) - x(2)); ...         % dx2/dt
     x(3) - x(4) - C*x(3)*x(4) - x(3)*x(4)^2; ...                                  % dx3/dt
     H*x(3) - 3*x(4) + C*x(3)*x(4) + x(3)*x(4)^2 + 2*beta*(x(2) - x(4)) ...        % dx4/dt
 ];

 % Solve the ODE system
 [t, x_out] = ode15s(ode_fun, tspan, x0, options);

 % Apply Equation (4): Linear combination of the four states
 ECG = alpha(1)*x_out(:,1) + alpha(2)*x_out(:,2) + alpha(3)*x_out(:,3) + alpha(4)*x_out(:,4);
```

__Continous System with `ode15s`__

<img width="650" alt="image" src="https://github.com/user-attachments/assets/4ff2cd60-d0e0-4356-a5e2-f61690b9847b" />

__Discrete System with implicit Tustin__

<img width="650" alt="image" src="https://github.com/user-attachments/assets/dc011890-06a7-4f2d-9079-65848f88b782" />

__Implementation example in STM32U585CI__

![nyutnyut-ezgif com-crop](https://github.com/user-attachments/assets/1fc92dd8-1633-4237-8016-77b70524828f)


