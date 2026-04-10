/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2026 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */

/* Includes ------------------------------------------------------------------*/
#include "main.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include "rd_model.h"
#include "ssd1306.h"
#include "ssd1306_fonts.h"
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
#define RD_GAIN_MV_PER_MV   (1200.0f)
#define RD_OFFSET_MV        (1700.0f)
#define VDDA_MV             (3300.0f)
#define DAC_BITS            (12)
#define DAC_MAX             ((1U << DAC_BITS) - 1U)
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
ADC_HandleTypeDef hadc1;
DAC_HandleTypeDef hdac1;
I2C_HandleTypeDef hi2c1;

/* USER CODE BEGIN PV */
static RDParams params;
static RDContext rdCtx;

static uint16_t *dacBuf = NULL;
static int dacLen = 0;
static int playIndex = 0;
static uint16_t refDacCode = 0;

static uint32_t last_sample_us = 0;

static int last_adc_sig = 0;
static int last_adc_ref = 0;
static int last_adc_diff = 0;

static uint32_t build_time_ms = 0;
/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_ADC1_Init(void);
static void MX_DAC1_Init(void);
static void MX_I2C1_Init(void);

/* USER CODE BEGIN PFP */
static uint16_t rd_mv_to_dac(float mv);
static bool build_dac_block(void);
static uint16_t read_adc_channel(uint32_t channel);
static void oled_show_build_info(void);
static void oled_show_error(const char *msg);

static void dwt_delay_init(void);
static uint32_t micros_dwt(void);
static const char* rd_solver_name(void);
/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

static void dwt_delay_init(void)
{
  CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
  DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
  DWT->CYCCNT = 0;
}

static uint32_t micros_dwt(void)
{
  return (uint32_t)(DWT->CYCCNT / (SystemCoreClock / 1000000UL));
}

static const char* rd_solver_name(void)
{
#if defined(RD_SOLVER_RK4)
  return "RD RK4";
#elif defined(RD_SOLVER_TUSTIN)
  return "RD Tustin";
#else
  return "RD ?";
#endif
}

static uint16_t rd_mv_to_dac(float mv)
{
  float v_mv = RD_OFFSET_MV + RD_GAIN_MV_PER_MV * mv;

  if (v_mv < 0.0f) v_mv = 0.0f;
  if (v_mv > VDDA_MV) v_mv = VDDA_MV;

  int code = (int)lroundf((v_mv / VDDA_MV) * (float)DAC_MAX);
  if (code < 0) code = 0;
  if (code > (int)DAC_MAX) code = (int)DAC_MAX;
  return (uint16_t)code;
}

static bool build_dac_block(void)
{
  float *rdBlockMv = NULL;
  dacLen = 0;

  if (!build_rd_block_mv(&params, &rdBlockMv, &dacLen, &rdCtx)) {
    return false;
  }

  dacBuf = (uint16_t*)malloc((size_t)dacLen * sizeof(uint16_t));
  if (!dacBuf) {
    if (rdBlockMv) free(rdBlockMv);
    return false;
  }

  for (int i = 0; i < dacLen; i++) {
    dacBuf[i] = rd_mv_to_dac(rdBlockMv[i]);
  }

  free(rdBlockMv);
  return true;
}

static uint16_t read_adc_channel(uint32_t channel)
{
  ADC_ChannelConfTypeDef sConfig = {0};

  sConfig.Channel = channel;
  sConfig.Rank = ADC_REGULAR_RANK_1;
  sConfig.SamplingTime = ADC_SAMPLETIME_2CYCLES_5;
  sConfig.SingleDiff = ADC_SINGLE_ENDED;
  sConfig.OffsetNumber = ADC_OFFSET_NONE;
  sConfig.Offset = 0;

  if (HAL_ADC_ConfigChannel(&hadc1, &sConfig) != HAL_OK) {
    Error_Handler();
  }

  if (HAL_ADC_Start(&hadc1) != HAL_OK) {
    Error_Handler();
  }

  if (HAL_ADC_PollForConversion(&hadc1, 10) != HAL_OK) {
    Error_Handler();
  }

  uint16_t val = (uint16_t)HAL_ADC_GetValue(&hadc1);

  if (HAL_ADC_Stop(&hadc1) != HAL_OK) {
    Error_Handler();
  }

  return val;
}

static void oled_show_build_info(void)
{
  char line[32];
  uint32_t sysclk_mhz;
  uint32_t hclk_mhz;
  uint32_t pclk1_mhz;

  SystemCoreClockUpdate();
  sysclk_mhz = SystemCoreClock / 1000000UL;
  hclk_mhz   = HAL_RCC_GetHCLKFreq() / 1000000UL;
  pclk1_mhz  = HAL_RCC_GetPCLK1Freq() / 1000000UL;

  ssd1306_Fill(Black);

  ssd1306_SetCursor(0, 0);
  snprintf(line, sizeof(line), "SYS:%3lu H:%3lu",
           (unsigned long)sysclk_mhz,
           (unsigned long)hclk_mhz);
  ssd1306_WriteString(line, Font_7x10, White);

  ssd1306_SetCursor(0, 16);
  snprintf(line, sizeof(line), "P1:%3lu B:%4lu",
           (unsigned long)pclk1_mhz,
           (unsigned long)build_time_ms);
  ssd1306_WriteString(line, Font_7x10, White);

  ssd1306_SetCursor(0, 32);
  snprintf(line, sizeof(line), "Len:%d", dacLen);
  ssd1306_WriteString(line, Font_7x10, White);

  ssd1306_SetCursor(0, 48);
  snprintf(line, sizeof(line), "%s", rd_solver_name());
  ssd1306_WriteString(line, Font_7x10, White);

  ssd1306_UpdateScreen();
}

static void oled_show_error(const char *msg)
{
  ssd1306_Fill(Black);
  ssd1306_SetCursor(0, 0);
  ssd1306_WriteString("ERROR", Font_7x10, White);
  ssd1306_SetCursor(0, 16);
  ssd1306_WriteString((char*)msg, Font_7x10, White);
  ssd1306_UpdateScreen();
}

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
  /* USER CODE BEGIN 1 */
  char line[32];
  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_ADC1_Init();
  MX_DAC1_Init();
  MX_I2C1_Init();

  /* USER CODE BEGIN 2 */
  HAL_Delay(100);
  SystemCoreClockUpdate();
  dwt_delay_init();
  ssd1306_Init();

  ssd1306_Fill(Black);
  ssd1306_SetCursor(0, 0);
  snprintf(line, sizeof(line), "Building %s...", rd_solver_name());
  ssd1306_WriteString(line, Font_7x10, White);
  ssd1306_UpdateScreen();

  rd_init_default_params(&params);

  if (!rd_init(&params, &rdCtx)) {
    oled_show_error("RD init fail");
    Error_Handler();
  }

  uint32_t t0 = HAL_GetTick();
  if (!build_dac_block()) {
    oled_show_error("Build failed");
    Error_Handler();
  }
  build_time_ms = HAL_GetTick() - t0;

  oled_show_build_info();

  if (HAL_DAC_Start(&hdac1, DAC_CHANNEL_1) != HAL_OK) {
    oled_show_error("DAC1 CH1 fail");
    Error_Handler();
  }

  if (HAL_DAC_Start(&hdac1, DAC_CHANNEL_2) != HAL_OK) {
    oled_show_error("DAC1 CH2 fail");
    Error_Handler();
  }

  refDacCode = rd_mv_to_dac(0.0f);

  /* DAC1 CH1 = PA4 = reference */
  if (HAL_DAC_SetValue(&hdac1, DAC_CHANNEL_1, DAC_ALIGN_12B_R, refDacCode) != HAL_OK) {
    oled_show_error("Set DAC ref fail");
    Error_Handler();
  }

  /* DAC1 CH2 = PA5 = signal */
  if (HAL_DAC_SetValue(&hdac1, DAC_CHANNEL_2, DAC_ALIGN_12B_R, dacBuf[0]) != HAL_OK) {
    oled_show_error("Set DAC sig fail");
    Error_Handler();
  }

  last_sample_us = micros_dwt();
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  const uint32_t period_us = (uint32_t)lroundf(1000000.0f / params.ecg_fs);

  while (1)
  {
    if ((uint32_t)(micros_dwt() - last_sample_us) >= period_us)
    {
      last_sample_us += period_us;

      if (HAL_DAC_SetValue(&hdac1, DAC_CHANNEL_2, DAC_ALIGN_12B_R, dacBuf[playIndex]) != HAL_OK) {
        oled_show_error("DAC update fail");
        Error_Handler();
      }

      for (volatile int k = 0; k < 200; k++) { __NOP(); }

      last_adc_sig  = (int)read_adc_channel(ADC_CHANNEL_1); /* PA0 / ADC1_INP1 */
      last_adc_ref  = (int)read_adc_channel(ADC_CHANNEL_2); /* PA1 / ADC1_INP2 */
      last_adc_diff = last_adc_sig - last_adc_ref;

      playIndex++;
      if (playIndex >= dacLen) {
        playIndex = 0;
      }
    }
  }
  /* USER CODE END WHILE */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE0);
  while(!__HAL_PWR_GET_FLAG(PWR_FLAG_VOSRDY)) {}

  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_LSI|RCC_OSCILLATORTYPE_CSI;
  RCC_OscInitStruct.LSIState = RCC_LSI_ON;
  RCC_OscInitStruct.CSIState = RCC_CSI_ON;
  RCC_OscInitStruct.CSICalibrationValue = RCC_CSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLL1_SOURCE_CSI;
  RCC_OscInitStruct.PLL.PLLM = 1;
  RCC_OscInitStruct.PLL.PLLN = 125;
  RCC_OscInitStruct.PLL.PLLP = 2;
  RCC_OscInitStruct.PLL.PLLQ = 2;
  RCC_OscInitStruct.PLL.PLLR = 2;
  RCC_OscInitStruct.PLL.PLLRGE = RCC_PLL1_VCIRANGE_2;
  RCC_OscInitStruct.PLL.PLLVCOSEL = RCC_PLL1_VCORANGE_WIDE;
  RCC_OscInitStruct.PLL.PLLFRACN = 0;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2
                              |RCC_CLOCKTYPE_PCLK3;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV1;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;
  RCC_ClkInitStruct.APB3CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_5) != HAL_OK)
  {
    Error_Handler();
  }

  __HAL_FLASH_SET_PROGRAM_DELAY(FLASH_PROGRAMMING_DELAY_2);
}

/* USER CODE BEGIN 4 */
static void MX_ADC1_Init(void)
{
  ADC_ChannelConfTypeDef sConfig = {0};

  hadc1.Instance = ADC1;
  hadc1.Init.ClockPrescaler = ADC_CLOCK_ASYNC_DIV1;
  hadc1.Init.Resolution = ADC_RESOLUTION_12B;
  hadc1.Init.DataAlign = ADC_DATAALIGN_RIGHT;
  hadc1.Init.ScanConvMode = ADC_SCAN_DISABLE;
  hadc1.Init.EOCSelection = ADC_EOC_SINGLE_CONV;
  hadc1.Init.LowPowerAutoWait = DISABLE;
  hadc1.Init.ContinuousConvMode = DISABLE;
  hadc1.Init.NbrOfConversion = 1;
  hadc1.Init.DiscontinuousConvMode = DISABLE;
  hadc1.Init.ExternalTrigConv = ADC_SOFTWARE_START;
  hadc1.Init.ExternalTrigConvEdge = ADC_EXTERNALTRIGCONVEDGE_NONE;
  hadc1.Init.DMAContinuousRequests = DISABLE;
  hadc1.Init.SamplingMode = ADC_SAMPLING_MODE_NORMAL;
  hadc1.Init.Overrun = ADC_OVR_DATA_PRESERVED;
  hadc1.Init.OversamplingMode = DISABLE;
  if (HAL_ADC_Init(&hadc1) != HAL_OK)
  {
    Error_Handler();
  }

  sConfig.Channel = ADC_CHANNEL_1;
  sConfig.Rank = ADC_REGULAR_RANK_1;
  sConfig.SamplingTime = ADC_SAMPLETIME_2CYCLES_5;
  sConfig.SingleDiff = ADC_SINGLE_ENDED;
  sConfig.OffsetNumber = ADC_OFFSET_NONE;
  sConfig.Offset = 0;
  if (HAL_ADC_ConfigChannel(&hadc1, &sConfig) != HAL_OK)
  {
    Error_Handler();
  }
}

static void MX_DAC1_Init(void)
{
  DAC_ChannelConfTypeDef sConfig = {0};

  hdac1.Instance = DAC1;
  if (HAL_DAC_Init(&hdac1) != HAL_OK)
  {
    Error_Handler();
  }

  sConfig.DAC_HighFrequency = DAC_HIGH_FREQUENCY_INTERFACE_MODE_DISABLE;
  sConfig.DAC_DMADoubleDataMode = DISABLE;
  sConfig.DAC_SignedFormat = DISABLE;
  sConfig.DAC_SampleAndHold = DAC_SAMPLEANDHOLD_DISABLE;
  sConfig.DAC_Trigger = DAC_TRIGGER_NONE;
  sConfig.DAC_OutputBuffer = DAC_OUTPUTBUFFER_ENABLE;
  sConfig.DAC_ConnectOnChipPeripheral = DAC_CHIPCONNECT_EXTERNAL;
  sConfig.DAC_UserTrimming = DAC_TRIMMING_FACTORY;
  if (HAL_DAC_ConfigChannel(&hdac1, &sConfig, DAC_CHANNEL_1) != HAL_OK)
  {
    Error_Handler();
  }

  if (HAL_DAC_ConfigChannel(&hdac1, &sConfig, DAC_CHANNEL_2) != HAL_OK)
  {
    Error_Handler();
  }
}

static void MX_I2C1_Init(void)
{
  hi2c1.Instance = I2C1;
  hi2c1.Init.Timing = 0x00707CBB;
  hi2c1.Init.OwnAddress1 = 0;
  hi2c1.Init.AddressingMode = I2C_ADDRESSINGMODE_7BIT;
  hi2c1.Init.DualAddressMode = I2C_DUALADDRESS_DISABLE;
  hi2c1.Init.OwnAddress2 = 0;
  hi2c1.Init.OwnAddress2Masks = I2C_OA2_NOMASK;
  hi2c1.Init.GeneralCallMode = I2C_GENERALCALL_DISABLE;
  hi2c1.Init.NoStretchMode = I2C_NOSTRETCH_DISABLE;
  if (HAL_I2C_Init(&hi2c1) != HAL_OK)
  {
    Error_Handler();
  }

  if (HAL_I2CEx_ConfigAnalogFilter(&hi2c1, I2C_ANALOGFILTER_ENABLE) != HAL_OK)
  {
    Error_Handler();
  }

  if (HAL_I2CEx_ConfigDigitalFilter(&hi2c1, 0) != HAL_OK)
  {
    Error_Handler();
  }
}

static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};

  __HAL_RCC_GPIOA_CLK_ENABLE();
  __HAL_RCC_GPIOB_CLK_ENABLE();

  GPIO_InitStruct.Pin = GPIO_PIN_0 | GPIO_PIN_1;
  GPIO_InitStruct.Mode = GPIO_MODE_ANALOG;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(GPIOA, &GPIO_InitStruct);
}
/* USER CODE END 4 */

void Error_Handler(void)
{
  __disable_irq();
  while (1)
  {
  }
}

#ifdef USE_FULL_ASSERT
void assert_failed(uint8_t *file, uint32_t line)
{
  (void)file;
  (void)line;
}
#endif
