#include "rthybrid_komendantov_kononenko_1996_neuron_c.h"

#include <math.h>
#include <string.h>

enum {
  V = 0,
  M_B = 1,
  H_B = 2,
  M = 3,
  H = 4,
  N = 5,
  M_CA = 6,
  CA = 7,
};

enum {
  DT = 0,
  I = 1,
  SYN = 2,
  CM = 3,
  G_NA_V = 4,
  V_NA = 5,
  G_K = 6,
  V_K = 7,
  G_NA = 8,
  G_B = 9,
  V_B = 10,
  G_NA_TTX = 11,
  G_K_TEA = 12,
  G_CA = 13,
  V_CA = 14,
  G_CA_CA = 15,
  K_BETA = 16,
  BETA = 17,
  RHO = 18,
  K_S = 19,
};

static const double DTS[19] = {0.000010, 0.000020, 0.000030, 0.000040, 0.000050,
                               0.000060, 0.000070, 0.000080, 0.000090, 0.000100,
                               0.000200, 0.000300, 0.000400, 0.000500, 0.000600,
                               0.000700, 0.000800, 0.000900, 0.001000};

static const double PTS[19] = {
    489646.000000, 250715.000000, 172155.000000, 131029.714286, 105694.000000,
    88541.272727,  76167.307692,  66740.857143,  59461.687500,  53611.000000,
    26996.527778,  18041.654545,  13546.534247,  10844.543478,  9040.509091,
    7750.382812,   6783.319728,   6030.793939,   5428.152174};

typedef void (*rtsyn_rk4_deriv_fn_t)(const double *state, double *deriv, size_t n, void *user_data);
extern void rtsyn_plugin_rk4_step_n(double *state, size_t n, double dt, rtsyn_rk4_deriv_fn_t deriv_fn, void *user_data);

static double sanitize_syn(double v) {
  if (!isfinite(v)) {
    return 0.0;
  }
  return v;
}

static int vars_finite(const double vars[8]) {
  for (int i = 0; i < 8; ++i) {
    if (!isfinite(vars[i])) {
      return 0;
    }
  }
  return 1;
}

static double i_na_v(const double vars[8], const double p[20]) {
  return p[G_NA_V] * (1.0 / (1.0 + exp(-0.2 * (vars[V] + 45.0)))) *
         (vars[V] - p[V_NA]);
}

static double i_k(const double vars[8], const double p[20]) {
  return p[G_K] * (vars[V] - p[V_K]);
}

static double i_na(const double vars[8], const double p[20]) {
  return p[G_NA] * (vars[V] - p[V_NA]);
}

static double i_b(const double vars[8], const double p[20]) {
  return p[G_B] * vars[M_B] * vars[H_B] * (vars[V] - p[V_B]);
}

static double f_m_b(const double vars[8]) {
  return (1.0 / (1.0 + exp(0.4 * (vars[V] + 34.0))) - vars[M_B]) / 0.05;
}

static double f_h_b(const double vars[8]) {
  return (1.0 / (1.0 + exp(-0.55 * (vars[V] + 43.0))) - vars[H_B]) / 1.5;
}

static double i_na_ttx(const double vars[8], const double p[20]) {
  return p[G_NA_TTX] * vars[M] * vars[M] * vars[M] * vars[H] *
         (vars[V] - p[V_NA]);
}

static double i_k_tea(const double vars[8], const double p[20]) {
  return p[G_K_TEA] * vars[N] * vars[N] * vars[N] * vars[N] *
         (vars[V] - p[V_K]);
}

static double f_m(const double vars[8]) {
  return (1.0 / (1.0 + exp(-0.4 * (vars[V] + 31.0))) - vars[M]) / 0.0005;
}

static double f_h(const double vars[8]) {
  return (1.0 / (1.0 + exp(0.25 * (vars[V] + 45.0))) - vars[H]) / 0.01;
}

static double f_n(const double vars[8]) {
  return (1.0 / (1.0 + exp(-0.18 * (vars[V] + 25.0))) - vars[N]) / 0.015;
}

static double i_ca(const double vars[8], const double p[20]) {
  return p[G_CA] * vars[M_CA] * vars[M_CA] * (vars[V] - p[V_CA]);
}

static double f_m_ca(const double vars[8]) {
  return (1.0 / (1.0 + exp(-0.2 * vars[V])) - vars[M_CA]) / 0.01;
}

static double i_ca_ca(const double vars[8], const double p[20]) {
  return p[G_CA_CA] * (1.0 / (1.0 + exp(-0.06 * (vars[V] + 45.0)))) *
         (1.0 / (1.0 + exp(p[K_BETA] * (vars[CA] - p[BETA])))) *
         (vars[V] - p[V_CA]);
}

static double f_ca(const double vars[8], const double p[20]) {
  return p[RHO] * ((-i_ca(vars, p) / 808.310846) - (p[K_S] * vars[CA]));
}

static void eval(const double vars[8], double ret[8], const double p[20], double syn) {
  const double i_na_ttx_v = i_na_ttx(vars, p);
  const double i_k_tea_v = i_k_tea(vars, p);
  const double i_k_v = i_k(vars, p);
  const double i_na_vv = i_na(vars, p);
  const double i_na_v_gate = i_na_v(vars, p);
  const double i_b_v = i_b(vars, p);
  const double i_ca_v = i_ca(vars, p);
  const double i_ca_ca_v = i_ca_ca(vars, p);

  const double syn_sanitized = sanitize_syn(syn);

  ret[V] = (-(i_na_ttx_v + i_k_tea_v + i_k_v + i_na_vv + i_na_v_gate + i_b_v +
              i_ca_v + i_ca_ca_v) +
            p[I] - syn_sanitized) /
           p[CM];
  ret[M_B] = f_m_b(vars);
  ret[H_B] = f_h_b(vars);
  ret[M] = f_m(vars);
  ret[H] = f_h(vars);
  ret[N] = f_n(vars);
  ret[M_CA] = f_m_ca(vars);
  ret[CA] = f_ca(vars, p);
}

typedef struct {
  const double *params;
  double syn;
} rthybrid_komendantov_kononenko_1996_neuron_c_ctx_t;

static void rthybrid_komendantov_kononenko_1996_neuron_c_deriv(const double *state, double *deriv,
                                                                size_t n, void *user_data) {
  (void)n;
  if (state == NULL || deriv == NULL || user_data == NULL) {
    return;
  }
  const rthybrid_komendantov_kononenko_1996_neuron_c_ctx_t *ctx =
      (const rthybrid_komendantov_kononenko_1996_neuron_c_ctx_t *)user_data;
  eval(state, deriv, ctx->params, ctx->syn);
}

static double select_dt_neuron_model(rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s, double pts_live,
                                     double *chosen_pts) {
  double aux = pts_live;
  double factor = 1.0;
  double pts_sel = -1.0;
  double dt_sel = -1.0;
  int flag = 0;

  while (aux < PTS[0]) {
    aux = pts_live * factor;
    factor += 1.0;

    for (int i = 18; i >= 0; --i) {
      if (PTS[i] > aux) {
        dt_sel = DTS[i];
        pts_sel = PTS[i];
        {
          double ratio = pts_sel / pts_live;
          if ((ratio - floor(ratio)) <= 0.1 * floor(ratio)) {
            flag = 1;
          }
        }
        break;
      }
    }

    if (flag) {
      break;
    }
  }

  if (!flag) {
    for (int i = 18; i >= 0; --i) {
      if (PTS[i] > aux) {
        dt_sel = DTS[i];
        pts_sel = PTS[i];
        break;
      }
    }
  }

  s->params[DT] = dt_sel;
  *chosen_pts = pts_sel;
  return dt_sel;
}

static void refresh_burst(rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s) {
  double burst = s->burst_duration_value <= -1.0 ? s->input_burst_duration
                                                  : s->burst_duration_value;
  if (burst <= 0.0) {
    burst = 1e-9;
  }
  s->burst_duration = burst;

  {
    double pts = -1.0;
    double denom;
    select_dt_neuron_model(s, burst * s->freq, &pts);
    denom = burst * s->freq;
    if (denom <= 0.0) {
      s->s_points = 1;
    } else {
      long steps = (long)(pts / denom);
      if (steps < 1) {
        steps = 1;
      }
      s->s_points = (size_t)steps;
    }
  }
}

void rthybrid_komendantov_kononenko_1996_neuron_c_init(rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s) {
  memset(s, 0, sizeof(*s));

  s->burst_duration = 1.0;
  s->burst_duration_value = 1.0;
  s->freq = 1000.0;
  s->s_points = 1;
  s->input_burst_duration = 1.0;

  s->params[I] = 0.0;
  s->params[SYN] = 0.0;
  s->params[CM] = 0.02;
  s->params[G_NA_V] = 0.11;
  s->params[V_NA] = 40.0;
  s->params[G_NA] = 0.0231;
  s->params[G_NA_TTX] = 400.0;
  s->params[G_K] = 0.25;
  s->params[G_K_TEA] = 10.0;
  s->params[V_K] = -70.0;
  s->params[G_B] = 0.165;
  s->params[V_B] = -58.0;
  s->params[G_CA] = 1.5;
  s->params[V_CA] = 150.0;
  s->params[G_CA_CA] = 0.02;
  s->params[K_BETA] = 15000.0;
  s->params[BETA] = 0.00004;
  s->params[RHO] = 0.002;
  s->params[K_S] = 50.0;

  s->vars[V] = -55.0;
  s->vars[CA] = 0.0;
  s->vars[M_B] = 1.0 / (1.0 + exp(0.4 * (s->vars[V] + 34.0)));
  s->vars[H_B] = 1.0 / (1.0 + exp(-0.55 * (s->vars[V] + 43.0)));
  s->vars[M] = 1.0 / (1.0 + exp(-0.4 * (s->vars[V] + 31.0)));
  s->vars[H] = 1.0 / (1.0 + exp(0.25 * (s->vars[V] + 45.0)));
  s->vars[N] = 1.0 / (1.0 + exp(-0.18 * (s->vars[V] + 25.0)));
  s->vars[M_CA] = 1.0 / (1.0 + exp(-0.2 * s->vars[V]));

  refresh_burst(s);
}

void rthybrid_komendantov_kononenko_1996_neuron_c_set_config(rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s, const char *key, size_t len, double value) {
  if (len == 18 && strncmp(key, "Burst duration (s)", len) == 0) {
    s->burst_duration_value = value;
  } else if ((len == 1 && strncmp(key, "i", len) == 0) ||
             (len == 1 && strncmp(key, "I", len) == 0)) {
    s->params[I] = value;
  } else if (len == 2 && strncmp(key, "cm", len) == 0) {
    s->params[CM] = value;
  } else if (len == 6 && strncmp(key, "g_na_v", len) == 0) {
    s->params[G_NA_V] = value;
  } else if (len == 4 && strncmp(key, "v_na", len) == 0) {
    s->params[V_NA] = value;
  } else if (len == 4 && strncmp(key, "g_na", len) == 0) {
    s->params[G_NA] = value;
  } else if (len == 8 && strncmp(key, "g_na_ttx", len) == 0) {
    s->params[G_NA_TTX] = value;
  } else if (len == 3 && strncmp(key, "g_k", len) == 0) {
    s->params[G_K] = value;
  } else if (len == 7 && strncmp(key, "g_k_tea", len) == 0) {
    s->params[G_K_TEA] = value;
  } else if (len == 3 && strncmp(key, "v_k", len) == 0) {
    s->params[V_K] = value;
  } else if (len == 3 && strncmp(key, "g_b", len) == 0) {
    s->params[G_B] = value;
  } else if (len == 3 && strncmp(key, "v_b", len) == 0) {
    s->params[V_B] = value;
  } else if (len == 4 && strncmp(key, "g_ca", len) == 0) {
    s->params[G_CA] = value;
  } else if (len == 4 && strncmp(key, "v_ca", len) == 0) {
    s->params[V_CA] = value;
  } else if (len == 7 && strncmp(key, "g_ca_ca", len) == 0) {
    s->params[G_CA_CA] = value;
  } else if (len == 6 && strncmp(key, "k_beta", len) == 0) {
    s->params[K_BETA] = value;
  } else if (len == 4 && strncmp(key, "beta", len) == 0) {
    s->params[BETA] = value;
  } else if (len == 3 && strncmp(key, "rho", len) == 0) {
    s->params[RHO] = value;
  } else if (len == 3 && strncmp(key, "k_s", len) == 0) {
    s->params[K_S] = value;
  } else if (len == 8 && strncmp(key, "v0 (mV)", len) == 0) {
    s->vars[V] = value;
  }
}

void rthybrid_komendantov_kononenko_1996_neuron_c_set_input(rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s, const char *key, size_t len, double value) {
  if (len == 9 && strncmp(key, "Isyn (nA)", len) == 0) {
    s->input_syn = sanitize_syn(value);
  } else if (len == 18 && strncmp(key, "Burst duration (s)", len) == 0) {
    s->input_burst_duration = value;
  }
}

void rthybrid_komendantov_kononenko_1996_neuron_c_process(rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s, double period_seconds) {
  if (period_seconds <= 0.0) {
    return;
  }

  s->input_syn = sanitize_syn(s->input_syn);
  s->freq = 1.0 / period_seconds;
  refresh_burst(s);
  if (!isfinite(s->params[DT]) || s->params[DT] <= 0.0) {
    s->params[DT] = 0.0001;
  }

  rthybrid_komendantov_kononenko_1996_neuron_c_ctx_t ctx = {
      .params = s->params,
      .syn = s->input_syn,
  };

  for (size_t i = 0; i < s->s_points; ++i) {
    double prev[8];
    memcpy(prev, s->vars, sizeof(prev));
    rtsyn_plugin_rk4_step_n(s->vars, 8, s->params[DT],
                            rthybrid_komendantov_kononenko_1996_neuron_c_deriv, &ctx);
    if (!vars_finite(s->vars)) {
      memcpy(s->vars, prev, sizeof(prev));
      break;
    }
  }
}

double rthybrid_komendantov_kononenko_1996_neuron_c_get_output(const rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s, const char *key, size_t len) {
  if (len == 6 && strncmp(key, "Vm (v)", len) == 0) {
    return s->vars[V] / 1000.0;
  }
  if (len == 7 && strncmp(key, "Vm (mV)", len) == 0) {
    return s->vars[V];
  }
  return 0.0;
}
