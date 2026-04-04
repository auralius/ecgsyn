#include <Arduino.h>
#include <math.h>
#include "ecgsyn_model.h"

using namespace ecgsyn;

static const double ECG_CENTER_MV = 0.4;
static const double ECG_SCALE_X   = 1000.0;
static const double VDDA          = 3.3;
static const int    DAC_BITS      = 12;
static const int    ADC_BITS      = 12;
static const int    DAC_MAX       = (1 << DAC_BITS) - 1;

#ifndef DAC_PIN
#define DAC_PIN PA5
#endif

#ifndef ADC_PIN
#define ADC_PIN A0
#endif

EcgSynParams params;
EcgSynContext ecgCtx;

double* ecgBlockMv = nullptr;
uint16_t* dacBuf = nullptr;
int dacLen = 0;
volatile int playIndex = 0;

static uint16_t ecg_mv_to_dac(double ecg_mv) {
  double volts = (ecg_mv + ECG_CENTER_MV) * ECG_SCALE_X * 0.001;
  if (volts < 0.0) volts = 0.0;
  if (volts > VDDA) volts = VDDA;

  int code = (int)lround((volts / VDDA) * DAC_MAX);
  if (code < 0) code = 0;
  if (code > DAC_MAX) code = DAC_MAX;
  return (uint16_t)code;
}

static bool build_dac_block() {
  if (!build_block_mv(params, ecgBlockMv, dacLen, ecgCtx)) {
    return false;
  }

  dacBuf = (uint16_t*)malloc((size_t)dacLen * sizeof(uint16_t));
  if (!dacBuf) return false;

  for (int i = 0; i < dacLen; i++) {
    dacBuf[i] = ecg_mv_to_dac(ecgBlockMv[i]);
  }
  return true;
}

void setup() {
  Serial.begin(115200);
  delay(3000);

  analogWriteResolution(DAC_BITS);
  analogReadResolution(ADC_BITS);
  pinMode(DAC_PIN, OUTPUT);
  pinMode(ADC_PIN, INPUT);

  params.ecg_fs      = 256;
  params.internal_fs = 256;
  params.n_beats     = 16;
  params.hr_mean     = 60.0;
  params.hr_std      = 1.0;
  params.noise_mv    = 0.0;
  params.seed_init   = 1;

  Serial.println("Building faithful ECG block...");
  if (!build_dac_block()) {
    Serial.println("ECG build failed.");
    while (1) delay(1000);
  }

  Serial.print("dacLen = ");
  Serial.println(dacLen);
  Serial.println("Playback start.");
  Serial.println("idx,ecg_mv,dac_code,adc_code");
}

void loop() {
  static uint32_t last_us = 0;
  const uint32_t period_us = 1000000UL / params.ecg_fs;
  uint32_t now = micros();

  if ((uint32_t)(now - last_us) >= period_us) {
    last_us += period_us;

    const int idx = playIndex;
    const double ecg_mv = ecgBlockMv[idx];
    const uint16_t dac_code = dacBuf[idx];

    analogWrite(DAC_PIN, dac_code);
    delayMicroseconds(20);
    int adc_code = analogRead(ADC_PIN);

    Serial.print(idx);
    Serial.print(",");
    Serial.print(ecg_mv, 6);
    Serial.print(",");
    Serial.print(dac_code);
    Serial.print(",");
    Serial.println(adc_code);

    playIndex++;
    if (playIndex >= dacLen) playIndex = 0;
  }
}