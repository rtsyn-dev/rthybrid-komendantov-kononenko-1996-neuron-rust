#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double vars[8];
  double params[20];
  double burst_duration;
  double burst_duration_value;
  double freq;
  size_t s_points;
  double input_syn;
  double input_burst_duration;
} rthybrid_komendantov_kononenko_1996_neuron_c_state_t;

void rthybrid_komendantov_kononenko_1996_neuron_c_init(rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s);
void rthybrid_komendantov_kononenko_1996_neuron_c_set_config(rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s, const char *key, size_t len, double value);
void rthybrid_komendantov_kononenko_1996_neuron_c_set_input(rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s, const char *key, size_t len, double value);
void rthybrid_komendantov_kononenko_1996_neuron_c_process(rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s, double period_seconds);
double rthybrid_komendantov_kononenko_1996_neuron_c_get_output(const rthybrid_komendantov_kononenko_1996_neuron_c_state_t *s, const char *key, size_t len);

#ifdef __cplusplus
}
#endif
