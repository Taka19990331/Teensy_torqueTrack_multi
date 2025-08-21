// ========================== Reference torque profile ================
/*
float refT[] = {
  -0.024044409379305, -0.057379170811646, -0.082347726221913, -0.101766214933962,
  -0.113655220028023, -0.117663416354635, -0.115056381548052, -0.108405281732599,
  -0.099988540207362, -0.090879528747806, -0.081642672434475, -0.072527546212679,
  -0.063458960718276, -0.054249479016750, -0.044811952196179, -0.035216063994700,
  -0.025540612547182, -0.015911614705552, -0.006501919148448,  0.002432548302564,
   0.010745191351587,  0.018194139928874,  0.024496791284796,  0.029439389639198,
   0.032744632722659,  0.034241708096117,  0.034058915272989,  0.032568936986831,
   0.030153271810964,  0.027117638786424,  0.023741043145104,  0.020194030904439,
   0.016649395025736,  0.013355774514323,  0.010522464005195,  0.008247017600510,
   0.006424851063036,  0.004835031486204,  0.003273337451531,  0.001722633564767,
   0.000388970497158, -0.000470791700955, -0.000641985659551, -0.000067311244177,
   0.001179679523487,  0.002736159666071,  0.004406850039273,  0.007230227976409,
   0.009152320921940, -0.006956574288570
}; // example: hip profile

float refT[] = {
  -0.001627031189353, -0.010078134806660, -0.019665055287961, -0.028449318847842,
  -0.035490434652401, -0.040735705640484, -0.044278299173914, -0.046483940264661,
  -0.047862586378117, -0.048865240917389, -0.049824949487050, -0.050939523812327,
  -0.052302587460691, -0.053956909666418, -0.055916233613297, -0.058138984270470,
  -0.060535773640330, -0.062966925474941, -0.065243978655166, -0.067132944173804,
  -0.068349359525542, -0.068588261887857, -0.067571582397193, -0.064978413580525,
  -0.060575065155178, -0.054480797499138, -0.047138104002056, -0.039079301695244,
  -0.030977734228921, -0.023438668348530, -0.016788037942653, -0.011143652372293,
  -0.006623949318000, -0.003236190508745, -0.000910088917431,  0.000458748984804,
   0.001073982824656,  0.001214655788092,  0.001175950182825,  0.001181999947384,
   0.001345780427386,  0.001695567089574,  0.002206015874143,  0.002808956358927,
   0.003390201082491,  0.003804974224416,  0.003985610244484,  0.004088539502501,
   0.004036081121161,  0.001432784028159
}; // example: knee profile

int   indexRefT = 0;
int   refTsize  = sizeof(refT)/sizeof(refT[0]);

  // Update torque reference from the profile (example: scaling by 0.5)
  // Tp_ref = refT[indexRefT] * 0.75f;
  // indexRefT++;
  // if (indexRefT >= refTsize) indexRefT = 0;
*/

/*
Feedforward friction control (multi-motor) with per-motor velocity unwrap+LPF
and per-motor current command: friction FF, bias learning, damping, rate limit.

- Motors are addressed by names: 1h, 1k, 2h, 2k (M1H..M2K)
- Each motor has its own encoder state, control state, command and measurements
- No need to declare "how many motors are connected": we try all; if conn=false we skip
- Velocity LPF is inside computeOmegaFromEncLPF() (avoid double LPF outside)
*/
struct EncVelState;
struct CtrlState;
#include <SPI.h>
//#include <TimerOne.h>
#include <IntervalTimer.h>
IntervalTimer ctrlTimer;
#include <math.h>
// #pragma GCC diagnostic error "-Wshadow"   // (optional) turn shadowing into a compile error

// ========================== Pin definitions ==========================
#define LED1h 23
#define LED1k 22
#define LED2h 19
#define LED2k 18

#define SS1h 38
#define SS1k 7
#define SS2h 6
#define SS2k 5

// ========================== Motor identity ==========================
enum MotorId : uint8_t { M1H=0, M1K=1, M2H=2, M2K=3, MOTOR_COUNT=4 };

inline const char* motorName(MotorId m) {
  static const char* N[] = {"1h","1k","2h","2k"};
  return N[m];
}

// Chip select pins (match your wiring)
const int SS_PIN[MOTOR_COUNT] = { SS1h, SS1k, SS2h, SS2k };

// ========================== Motor / model params =====================
// Motor torque constant
const float Kt = 0.043f;        // [Nm/A]

// Encoder (19-bit SSI)
const int   ENC_BITS = 19;
const long  ENC_MOD  = 1L << ENC_BITS;            // 524288
const long  ENC_HALF = ENC_MOD >> 1;              // 262144
const float ENC2RAD  = 2.0f * PI / 524288.0f;     // [rad/count]

// Identified friction parameters (tune to your plant)
const float Tc_fric  = 0.011f;     // Coulomb friction [Nm]
const float Ts_fric  = 0.000f;     // Static/Stribeck peak [Nm] (>= Tc_fric for breakaway) — keep 0 if you don't want breakaway
const float omega_s  = 0.3469f;    // Stribeck velocity [rad/s]
const float alpha_s  = 0.444f;     // Stribeck shape factor [-]
const float B        = 0.000111f;  // Viscous friction [Nm/(rad/s)]

// ========================== Control frequency ========================
const int   control_freq_hz = 300;                 // [Hz]
const float dt              = 1.0f / control_freq_hz; // [s]

// ========================== Control tuning ===========================
// Feedforward-only control; rely on driver's inner current loop.
// We add: friction FF, bias learning, damping, current rate limit, and filtering.
const float K_fric         = 0.9f;     // Friction FF gain (0.7–1.2 typical)
const float I_MAX          = 1.0f;     // Current limit [A]
const float DI_MAX         = 0.02f;    // Max current step per ISR [A/tick] (~6 A/s @300Hz)
const float OMEGA_DB       = 0.01f;    // Deadband for sign stabilization [rad/s]
const float TAU_BIAS_MAX   = 0.06f;    // Clamp for bias integrator [Nm]

// Bias learning (slow offset rejection near neutral)
float bias_alpha           = 0.0f;     // IIR coefficient for slow bias (set in setup)
const float bias_cutoff_hz = 0.5f;     // 0.3–0.8 Hz reasonable
const float tau_ref_th     = 0.005f;   // |Tp_ref| threshold [Nm]
const float omega_th       = 0.1f;     // |omega| threshold [rad/s] (~5.7 deg/s)

// Virtual damping
const float Kd_visc = 0.00015f;        // [Nm/(rad/s)] start 0.00005–0.0001, here a bit higher
const float Tgate    = 0.005f;         // [Nm] gate center for |Tp_ref|; g=0.5 at Tgate
const float gate_p   = 2.0f;           // gate shape (>=1), larger => sharper

// ========================== Velocity filtering =======================
const float omega_fc_hz = 25.0f;       // 20–30 Hz typical @300Hz
float omega_alpha = 0.0f;              // IIR coefficient for velocity LPF (set in setup)

// ========================== Outlier rejection =======================
const float OMEGA_CAP   = 50.0f; // Physical limit [rad/s] (tame absurd spikes)
const long  DCOUNTS_CAP = (long)((OMEGA_CAP * dt) / ENC2RAD + 0.5f); // max counts per ISR

// ========================== SPI comm ================================
SPISettings settings(100000, MSBFIRST, SPI_MODE3); // 100 kHz
byte last_rx_buf[12];
byte last_tx_buf[12];

// ========================== Per-motor state =========================
// Encoder unwrap + velocity + LPF state
struct EncVelState {
  bool      initialized = false;
  long      prev_raw19  = 0;     // last 19-bit raw encoder sample
  long long acc         = 0;     // unwrapped cumulative counts
  long long acc_prev1   = 0;     // cum[n-1]
  long long acc_prev2   = 0;     // cum[n-2]
  uint8_t   warm        = 0;     // warm-up counter for 3-pt diff
  float     omega_lpf   = 0.0f;  // LPF state
  float     omega       = 0.0f;  // last filtered velocity [rad/s]
};

// Control-side state per motor
struct CtrlState {
  float tau_bias      = 0.0f;  // learned torque bias [Nm]
  float cur_cmd_prev  = 0.0f;  // previous command for rate limiting [A]
  int   last_sign     = 0;     // sign memory for stabilized friction sign
  float i_meas_lpf    = 0.0f;  // measured current LPF [A]
};

// Arrays indexed by MotorId
EncVelState enc_state[MOTOR_COUNT];
CtrlState   ctrl_state[MOTOR_COUNT];

bool  conn_m   [MOTOR_COUNT] = {};  // SPI connection (checksum match)
long  enc_raw  [MOTOR_COUNT] = {};  // raw encoder (int32 from driver)
float cur_meas [MOTOR_COUNT] = {};  // measured current [A] from driver
float omega_m  [MOTOR_COUNT] = {};  // filtered angular velocity [rad/s]
float cur_cmd_m[MOTOR_COUNT] = {};  // command current [A] (to send)

// Reference torque per motor [Nm] (set these from your application/loop)
float Tp_ref_m[MOTOR_COUNT] = {0,0,0,0};
#define TPREF_1H Tp_ref_m[M1H]
#define TPREF_1K Tp_ref_m[M1K]
#define TPREF_2H Tp_ref_m[M2H]
#define TPREF_2K Tp_ref_m[M2K]

// Named aliases for readability (optional)
#define OMEGA_1H   (omega_m[M1H])
#define OMEGA_1K   (omega_m[M1K])
#define OMEGA_2H   (omega_m[M2H])
#define OMEGA_2K   (omega_m[M2K])
#define CURCMD_1H  (cur_cmd_m[M1H])
#define CURCMD_1K  (cur_cmd_m[M1K])
#define CURCMD_2H  (cur_cmd_m[M2H])
#define CURCMD_2K  (cur_cmd_m[M2K])

// ========================== Helpers ================================
inline float iirAlpha(float cutoff_hz, float Ts) {
  // First-order IIR in exponential form: y = a*y + (1-a)*x, with a = exp(-2*pi*fc*Ts)
  return expf(-2.0f * PI * cutoff_hz * Ts);
}

inline long unwrapEncDelta(long now, long prev) {
  // Unwrap delta between two modular encoder samples
  long d = now - prev;
  if (d >  ENC_HALF) d -= ENC_MOD;
  if (d < -ENC_HALF) d += ENC_MOD;
  return d;
}

// Return only friction magnitude of Stribeck + Coulomb (no viscous part)
inline float frictionMag(float w_abs) {
  // (Ts - Tc)*exp(-( |w|/ws )^alpha) + Tc
  float strb = (Ts_fric - Tc_fric) * expf(-powf(w_abs / omega_s, alpha_s));
  return (Tc_fric + strb);
}

// Smooth gate: g(|Tp_ref|) in [0..1], 1 near zero torque, ->0 for large |Tp_ref|
inline float gate_by_ref(float Tp_ref_abs) {
  float r = Tp_ref_abs / Tgate;
  return 1.0f / (1.0f + powf(r, gate_p));   // 1/(1+r^p)
}

// ========================== Velocity function ======================
// Returns filtered angular velocity [rad/s] for a motor.
inline float computeOmegaFromEncLPF(long enc_raw_32, EncVelState& st, float dt_, float alpha_lpf) {
  long enc_now = (long)(enc_raw_32 & (ENC_MOD - 1));
  if (!st.initialized) {
    st.prev_raw19 = enc_now;
    st.acc = st.acc_prev1 = st.acc_prev2 = 0;
    st.warm = 0;
    st.omega_lpf = st.omega = 0.0f;
    st.initialized = true;
  }

  long d_counts = unwrapEncDelta(enc_now, st.prev_raw19);
  st.prev_raw19 = enc_now;

  // Reject non-physical jumps (counts per ISR)
  if (labs(d_counts) > DCOUNTS_CAP) d_counts = 0;

  // Update cumulative counts
  st.acc += d_counts;

  // 3-point backward difference (2nd-order); simple diff during warm-up
  float omega_raw;
  if (st.warm < 2) {
    omega_raw = (d_counts * ENC2RAD) / dt_;
    st.warm++;
  } else {
    long long num_counts = 3LL * st.acc - 4LL * st.acc_prev1 + st.acc_prev2;
    omega_raw = ((float)num_counts * ENC2RAD) / (2.0f * dt_);
  }

  // Update history (exactly once)
  st.acc_prev2 = st.acc_prev1;
  st.acc_prev1 = st.acc;

  // Optional clamp to tame absurd spikes
  #ifdef OMEGA_CAP
  if (fabsf(omega_raw) > OMEGA_CAP) {
    omega_raw = copysignf(OMEGA_CAP, omega_raw);
  }
  #endif

  // Single-pole IIR LPF
  st.omega     = alpha_lpf * st.omega_lpf + (1.0f - alpha_lpf) * omega_raw;
  st.omega_lpf = st.omega;
  return st.omega; // filtered ω
}

// ========================== Current command ========================
// Computes current command [A] for one motor.
// Inputs: Tp_ref [Nm], filtered omega [rad/s], measured current (raw) [A].
// Uses per-motor CtrlState for bias-learning, sign memory, rate limiting.
inline float computeCurrentCommand(float Tp_ref, float omega, float i_meas_raw,
                                   CtrlState& cs, float alpha_i_meas) {
  // 1) Measured torque (LPF current for robustness)
  cs.i_meas_lpf = alpha_i_meas * cs.i_meas_lpf + (1.0f - alpha_i_meas) * i_meas_raw;
  float Tm_meas = Kt * cs.i_meas_lpf;

  // 2) Stabilized sign for Coulomb/Stribeck
  int sgn_w;
  if (fabsf(omega) > OMEGA_DB) {
    sgn_w = (omega > 0) - (omega < 0);
    cs.last_sign = sgn_w;
  } else {
    int sgn_ref = (Tp_ref > 0) - (Tp_ref < 0);
    sgn_w = (sgn_ref != 0) ? sgn_ref : cs.last_sign;
  }

  // 3) Friction estimate (Coulomb/Stribeck + viscous)
  float Tfric = sgn_w * frictionMag(fabsf(omega)) + B * omega;

  // 4) Output-side torque estimate (for bias learning)
  float Tp_meas = Tm_meas - Tfric;

  // 5) Bias learning near neutral (very slow IIR)
  if (fabsf(Tp_ref) < tau_ref_th && fabsf(omega) < omega_th) {
    float err_slow = (Tp_meas - Tp_ref);
    cs.tau_bias = bias_alpha * cs.tau_bias + (1.0f - bias_alpha) * err_slow;
    // Clamp to avoid runaway
    if (cs.tau_bias >  TAU_BIAS_MAX) cs.tau_bias =  TAU_BIAS_MAX;
    if (cs.tau_bias < -TAU_BIAS_MAX) cs.tau_bias = -TAU_BIAS_MAX;
  }

  // 6) Virtual damping (gated by |Tp_ref| → 0 near zero)
  float g = gate_by_ref(fabsf(Tp_ref));
  float T_damp = -Kd_visc * omega * g;

  // 7) Torque → current target
  float Tm_cmd_target  = Tp_ref + K_fric * Tfric - cs.tau_bias + T_damp;
  float cur_cmd_target = Tm_cmd_target / Kt;

  // 8) Rate limit & saturation
  float di = cur_cmd_target - cs.cur_cmd_prev;
  if (di >  DI_MAX) di =  DI_MAX;
  if (di < -DI_MAX) di = -DI_MAX;
  float cur_cmd = cs.cur_cmd_prev + di;
  if (cur_cmd >  I_MAX) cur_cmd =  I_MAX;
  if (cur_cmd < -I_MAX) cur_cmd = -I_MAX;
  cs.cur_cmd_prev = cur_cmd;

  return cur_cmd;
}

// ========================== SPI packet I/O ==========================
void send_current_command(int SS, bool servo_on, float current_cmd, bool* connection, long* val_enc, float* val_cur) {
  byte tx_buf[12] = {0};
  byte rx_buf[12] = {0};

  // Build TX
  tx_buf[0] = servo_on ? 0x08 : 0x00;                      // header
  memcpy(&tx_buf[3], &current_cmd, sizeof(float));         // current cmd [3..6]

  // Checksum (sum of 16-bit words over first 10 bytes)
  uint16_t checksum = 0;
  for (int i = 0; i < 10; i += 2) {
    uint16_t w = tx_buf[i] | (tx_buf[i + 1] << 8);
    checksum += w;
  }
  checksum &= 0x7FFF;
  checksum |= 0x8000;
  tx_buf[10] = checksum & 0xFF;
  tx_buf[11] = (checksum >> 8) & 0xFF;

  // SPI transfer
  SPI.beginTransaction(settings);
  delayMicroseconds(5);
  digitalWrite(SS, LOW);
  delayMicroseconds(10);
  for (int i = 0; i < 12; i++) {
    rx_buf[i] = SPI.transfer(tx_buf[i]);
  }
  delayMicroseconds(10);
  digitalWrite(SS, HIGH);
  SPI.endTransaction();

  // Save last TX/RX (optional)
  memcpy(last_tx_buf, tx_buf, 12);
  memcpy(last_rx_buf, rx_buf, 12);

  // Parse checksum
  uint16_t rx_checksum = rx_buf[10] | (rx_buf[11] << 8);
  uint16_t computed = 0;
  for (int i = 0; i < 10; i += 2) {
    uint16_t w = rx_buf[i] | (rx_buf[i + 1] << 8);
    computed += w;
  }
  computed &= 0x7FFF;
  computed |= 0x8000;

  *connection = (rx_checksum == computed);

  // Parse encoder (bytes [1..4])
  int32_t enc_raw = *((int32_t*)&rx_buf[1]);
  *val_enc = enc_raw;

  // Parse measured current as float (bytes [5..8])
  memcpy(val_cur, &rx_buf[5], sizeof(float));
}

// ========================== Setup ==================================
void setup() {
  pinMode(10, OUTPUT);           // SS must be OUTPUT to keep SPI master mode
  pinMode(11, OUTPUT);           // MOSI
  pinMode(12, INPUT);            // MISO
  pinMode(13, OUTPUT);           // SCK

  // Set all CS pins OUTPUT/HIGH (safe idle)
  for (uint8_t m=0; m<MOTOR_COUNT; ++m) {
    pinMode(SS_PIN[m], OUTPUT);
    digitalWrite(SS_PIN[m], HIGH);
  }

  pinMode(LED1h, OUTPUT);
  pinMode(LED1k, OUTPUT);
  pinMode(LED2h, OUTPUT);
  pinMode(LED2k, OUTPUT);

  // Precompute IIR alphas
  omega_alpha = iirAlpha(omega_fc_hz, dt);
  bias_alpha  = iirAlpha(bias_cutoff_hz, dt);

  Serial.begin(57600);
  SPI.begin();
  delay(100);

  // (Optional) quick connection warm-up: try once per motor
  for (uint8_t m=0; m<MOTOR_COUNT; ++m) {
    bool ok=false; long e=0; float i=0;
    send_current_command(SS_PIN[m], 0, 0.0f, &ok, &e, &i);
    conn_m[m] = ok;
  }

  // Start control ISR
  ctrlTimer.begin(timerISR, 1000000 / control_freq_hz); // period in us
  //Timer1.initialize(1000000 / control_freq_hz);
  //Timer1.attachInterrupt(timerISR);
}

// ========================== Main loop ===============================
void loop() {
  // Example: update per-motor torque references here if you want
  // TPREF_1H = 0.00f;
  TPREF_1K = 0.0f;
  // TPREF_2H = -0.01f;
  // TPREF_2K = 0.00f;

  delay(100);
}

// ========================== Control ISR =============================
void timerISR() {
  // 1) IO round: send PREVIOUS command and read back for ALL motors
  for (uint8_t m = 0; m < MOTOR_COUNT; ++m) {
    send_current_command(SS_PIN[m], 1, cur_cmd_m[m],
                         &conn_m[m], &enc_raw[m], &cur_meas[m]);
  }

  // 2) Per-motor compute: velocity and control (skip if not connected)
  const float i_alpha = iirAlpha(30.0f, dt); // current measurement LPF ~30 Hz
  for (uint8_t m = 0; m < MOTOR_COUNT; ++m) {
    if (!conn_m[m]) {            // hardware temporarily disabled/disconnected
      cur_cmd_m[m] = 0.0f;       // keep safe neutral
      continue;                  // skip this motor entirely
    }

    // Filtered velocity
    omega_m[m] = computeOmegaFromEncLPF(enc_raw[m], enc_state[m], dt, omega_alpha);

    // Compute command current from per-motor reference torque
    float Tp_ref = Tp_ref_m[m]; // [Nm]
    cur_cmd_m[m] = computeCurrentCommand(Tp_ref, omega_m[m], cur_meas[m],
                                         ctrl_state[m], i_alpha);
  }

  // 3) (Optional) lightweight debug at ~50 Hz
  // static uint16_t dbg=0;
  // if(++dbg >= 6) {
  //   dbg = 0;
  //   Serial.print("1h: conn="); Serial.print(conn_m[M1H]);
  //   Serial.print(" w="); Serial.print(OMEGA_1H,4);
  //   Serial.print(" i_cmd="); Serial.println(CURCMD_1H,4);
  // }
}
