#ifndef MOTOR_STATES_H
#define MOTOR_STATES_H

struct EncVelState {
  bool      initialized = false;
  long      prev_raw19  = 0;
  long long acc         = 0;
  long long acc_prev1   = 0;
  long long acc_prev2   = 0;
  uint8_t   warm        = 0;
  float     omega_lpf   = 0.0f;
  float     omega       = 0.0f;
};

struct CtrlState {
  float tau_bias      = 0.0f;
  float cur_cmd_prev  = 0.0f;
  int   last_sign     = 0;
  float i_meas_lpf    = 0.0f;
};

#endif
