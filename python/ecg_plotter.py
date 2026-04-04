import serial
import matplotlib.pyplot as plt
from collections import deque

PORT = "/dev/ttyACM0"
BAUD = 921600
MAX_POINTS = 1000
PLOT_EVERY = 25   # update plot every 25 received lines

ser = serial.Serial(PORT, BAUD, timeout=1)

idx_buf = deque(maxlen=MAX_POINTS)
ecg_mv_buf = deque(maxlen=MAX_POINTS)
dac_buf = deque(maxlen=MAX_POINTS)
adc_buf = deque(maxlen=MAX_POINTS)

plt.ion()

fig1, ax1 = plt.subplots()
line_mv, = ax1.plot([], [])
ax1.set_xlabel("Sample")
ax1.set_ylabel("ECG (mV)")
ax1.set_title("Generated ECG")
ax1.grid(True)
ax1.set_xlim(0, MAX_POINTS)
ax1.set_ylim(-0.6, 1.4)   # fixed range for ECG

fig2, ax2 = plt.subplots()
line_dac, = ax2.plot([], [], label="DAC code")
line_adc, = ax2.plot([], [], label="ADC code")
ax2.set_xlabel("Sample")
ax2.set_ylabel("Code")
ax2.set_title("DAC and ADC loopback")
ax2.legend()
ax2.grid(True)
ax2.set_xlim(0, MAX_POINTS)
ax2.set_ylim(0, 4095)     # fixed ADC/DAC range

# Wait for CSV header
while True:
    line = ser.readline().decode(errors="ignore").strip()
    if line:
        print("Received line: ", line)
    if line == "idx,ecg_mv,dac_code,adc_code":
        break

count = 0

try:
    while True:
        line = ser.readline().decode(errors="ignore").strip()
        if not line:
            print("Empty line received")
            continue

        parts = line.split(",")
        if len(parts) != 4:
            print("Invalid line: ", line)
            continue

        try:
            idx = int(parts[0])
            ecg_mv = float(parts[1])
            dac_code = int(parts[2])
            adc_code = int(parts[3])
        except ValueError:
            continue

        idx_buf.append(idx)
        ecg_mv_buf.append(ecg_mv)
        dac_buf.append(dac_code)
        adc_buf.append(adc_code)

        count += 1
        if count < PLOT_EVERY:
            continue
        count = 0

        x = range(len(idx_buf))

        line_mv.set_data(x, ecg_mv_buf)
        line_dac.set_data(x, dac_buf)
        line_adc.set_data(x, adc_buf)

        fig1.canvas.draw_idle()
        fig2.canvas.draw_idle()
        plt.pause(0.001)

except KeyboardInterrupt:
    pass
finally:
    ser.close()