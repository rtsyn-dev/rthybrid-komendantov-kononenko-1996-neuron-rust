use rtsyn_plugin::prelude::*;
use serde_json::Value;

const V: usize = 0;
const M_B: usize = 1;
const H_B: usize = 2;
const M: usize = 3;
const H: usize = 4;
const N: usize = 5;
const M_CA: usize = 6;
const CA: usize = 7;

const DT: usize = 0;
const I: usize = 1;
const CM: usize = 3;
const G_NA_V: usize = 4;
const V_NA: usize = 5;
const G_K: usize = 6;
const V_K: usize = 7;
const G_NA: usize = 8;
const G_B: usize = 9;
const V_B: usize = 10;
const G_NA_TTX: usize = 11;
const G_K_TEA: usize = 12;
const G_CA: usize = 13;
const V_CA: usize = 14;
const G_CA_CA: usize = 15;
const K_BETA: usize = 16;
const BETA: usize = 17;
const RHO: usize = 18;
const K_S: usize = 19;

#[derive(Debug)]
struct RthybridKomendantovKononenko1996NeuronRust {
    input_syn: f64,
    input_burst_duration: f64,
    out_0: f64,
    out_1: f64,
    vars: [f64; 8],
    params: [f64; 20],
    burst_duration_value: f64,
    freq: f64,
    s_points: usize,
}

impl Default for RthybridKomendantovKononenko1996NeuronRust {
    fn default() -> Self {
        let v0: f64 = -55.0;
        let mut vars = [0.0; 8];
        vars[V] = v0;
        vars[CA] = 0.0;
        vars[M_B] = 1.0 / (1.0 + (0.4 * (v0 + 34.0)).exp());
        vars[H_B] = 1.0 / (1.0 + (-0.55 * (v0 + 43.0)).exp());
        vars[M] = 1.0 / (1.0 + (-0.4 * (v0 + 31.0)).exp());
        vars[H] = 1.0 / (1.0 + (0.25 * (v0 + 45.0)).exp());
        vars[N] = 1.0 / (1.0 + (-0.18 * (v0 + 25.0)).exp());
        vars[M_CA] = 1.0 / (1.0 + (-0.2 * v0).exp());

        let mut params = [0.0; 20];
        params[I] = 0.0;
        params[CM] = 0.02;
        params[G_NA_V] = 0.11;
        params[V_NA] = 40.0;
        params[G_NA] = 0.0231;
        params[G_NA_TTX] = 400.0;
        params[G_K] = 0.25;
        params[G_K_TEA] = 10.0;
        params[V_K] = -70.0;
        params[G_B] = 0.165;
        params[V_B] = -58.0;
        params[G_CA] = 1.5;
        params[V_CA] = 150.0;
        params[G_CA_CA] = 0.02;
        params[K_BETA] = 15000.0;
        params[BETA] = 0.00004;
        params[RHO] = 0.002;
        params[K_S] = 50.0;
        params[DT] = 0.001;

        Self {
            input_syn: 0.0,
            input_burst_duration: 1.0,
            out_0: 0.0,
            out_1: 0.0,
            vars,
            params,
            burst_duration_value: 1.0,
            freq: 1000.0,
            s_points: 1,
        }
    }
}

impl RthybridKomendantovKononenko1996NeuronRust {
    fn sanitize_syn(v: f64) -> f64 {
        if v.is_finite() { v } else { 0.0 }
    }

    fn i_na_v(vars: &[f64; 8], p: &[f64; 20]) -> f64 {
        p[G_NA_V] * (1.0 / (1.0 + (-0.2 * (vars[V] + 45.0)).exp())) * (vars[V] - p[V_NA])
    }

    fn i_k(vars: &[f64; 8], p: &[f64; 20]) -> f64 {
        p[G_K] * (vars[V] - p[V_K])
    }

    fn i_na(vars: &[f64; 8], p: &[f64; 20]) -> f64 {
        p[G_NA] * (vars[V] - p[V_NA])
    }

    fn i_b(vars: &[f64; 8], p: &[f64; 20]) -> f64 {
        p[G_B] * vars[M_B] * vars[H_B] * (vars[V] - p[V_B])
    }

    fn i_na_ttx(vars: &[f64; 8], p: &[f64; 20]) -> f64 {
        p[G_NA_TTX] * vars[M].powi(3) * vars[H] * (vars[V] - p[V_NA])
    }

    fn i_k_tea(vars: &[f64; 8], p: &[f64; 20]) -> f64 {
        p[G_K_TEA] * vars[N].powi(4) * (vars[V] - p[V_K])
    }

    fn i_ca(vars: &[f64; 8], p: &[f64; 20]) -> f64 {
        p[G_CA] * vars[M_CA].powi(2) * (vars[V] - p[V_CA])
    }

    fn i_ca_ca(vars: &[f64; 8], p: &[f64; 20]) -> f64 {
        p[G_CA_CA]
            * (1.0 / (1.0 + (-0.06 * (vars[V] + 45.0)).exp()))
            * (1.0 / (1.0 + (p[K_BETA] * (vars[CA] - p[BETA])).exp()))
            * (vars[V] - p[V_CA])
    }

    fn eval(vars: &[f64; 8], der: &mut [f64; 8], p: &[f64; 20], syn: f64) {
        let syn = Self::sanitize_syn(syn);
        let i_na_ttx = Self::i_na_ttx(vars, p);
        let i_k_tea = Self::i_k_tea(vars, p);
        let i_k = Self::i_k(vars, p);
        let i_na = Self::i_na(vars, p);
        let i_na_v = Self::i_na_v(vars, p);
        let i_b = Self::i_b(vars, p);
        let i_ca = Self::i_ca(vars, p);
        let i_ca_ca = Self::i_ca_ca(vars, p);

        der[V] = (-(i_na_ttx + i_k_tea + i_k + i_na + i_na_v + i_b + i_ca + i_ca_ca)
            + p[I]
            - syn)
            / p[CM];
        der[M_B] = (1.0 / (1.0 + (0.4 * (vars[V] + 34.0)).exp()) - vars[M_B]) / 0.05;
        der[H_B] = (1.0 / (1.0 + (-0.55 * (vars[V] + 43.0)).exp()) - vars[H_B]) / 1.5;
        der[M] = (1.0 / (1.0 + (-0.4 * (vars[V] + 31.0)).exp()) - vars[M]) / 0.0005;
        der[H] = (1.0 / (1.0 + (0.25 * (vars[V] + 45.0)).exp()) - vars[H]) / 0.01;
        der[N] = (1.0 / (1.0 + (-0.18 * (vars[V] + 25.0)).exp()) - vars[N]) / 0.015;
        der[M_CA] = (1.0 / (1.0 + (-0.2 * vars[V]).exp()) - vars[M_CA]) / 0.01;
        der[CA] = p[RHO] * ((-i_ca / 808.310846) - (p[K_S] * vars[CA]));
    }

    fn select_dt_neuron_model(pts_live: f64) -> (f64, f64) {
        const DTS: [f64; 19] = [
            0.000010, 0.000020, 0.000030, 0.000040, 0.000050, 0.000060, 0.000070, 0.000080,
            0.000090, 0.000100, 0.000200, 0.000300, 0.000400, 0.000500, 0.000600, 0.000700,
            0.000800, 0.000900, 0.001000,
        ];
        const PTS: [f64; 19] = [
            489646.000000,
            250715.000000,
            172155.000000,
            131029.714286,
            105694.000000,
            88541.272727,
            76167.307692,
            66740.857143,
            59461.687500,
            53611.000000,
            26996.527778,
            18041.654545,
            13546.534247,
            10844.543478,
            9040.509091,
            7750.382812,
            6783.319728,
            6030.793939,
            5428.152174,
        ];

        let mut aux = pts_live;
        let mut factor = 1.0;
        let mut chosen = (DTS[DTS.len() - 1], PTS[PTS.len() - 1]);
        let mut found = false;

        while aux < PTS[0] {
            aux = pts_live * factor;
            factor += 1.0;
            for i in (0..PTS.len()).rev() {
                if PTS[i] > aux {
                    chosen = (DTS[i], PTS[i]);
                    let ratio = chosen.1 / pts_live;
                    if (ratio - ratio.floor()) <= 0.1 * ratio.floor() {
                        found = true;
                    }
                    break;
                }
            }
            if found {
                break;
            }
        }

        if !found {
            for i in (0..PTS.len()).rev() {
                if PTS[i] > aux {
                    chosen = (DTS[i], PTS[i]);
                    break;
                }
            }
        }

        chosen
    }

    fn refresh_burst(&mut self) {
        let mut burst = if self.burst_duration_value <= -1.0 {
            self.input_burst_duration
        } else {
            self.burst_duration_value
        };
        if !burst.is_finite() || burst <= 0.0 {
            burst = 1e-9;
        }

        let pts_live = burst * self.freq;
        let (dt, pts) = Self::select_dt_neuron_model(pts_live);
        self.params[DT] = dt;

        if pts_live <= 0.0 {
            self.s_points = 1;
        } else {
            let steps = (pts / pts_live) as usize;
            self.s_points = steps.max(1);
        }
    }
}

impl PluginDescriptor for RthybridKomendantovKononenko1996NeuronRust {
    fn name() -> &'static str {
        "RTHybrid Komendantov-Kononenko 1996 Neuron"
    }

    fn kind() -> &'static str {
        "rthybrid_komendantov_kononenko_1996_neuron"
    }

    fn plugin_type() -> PluginType {
        PluginType::Computational
    }

    fn inputs() -> &'static [&'static str] {
        &["Isyn (nA)", "Burst duration (s)"]
    }

    fn outputs() -> &'static [&'static str] {
        &["Vm (v)", "Vm (mV)"]
    }

    fn internal_variables() -> &'static [&'static str] {
        &["m", "h", "n", "ca"]
    }

    fn default_vars() -> Vec<(&'static str, Value)> {
        vec![
            ("Burst duration (s)", 1.0.into()),
            ("I", 0.0.into()),
            ("cm", 0.02.into()),
            ("g_na_v", 0.11.into()),
            ("v_na", 40.0.into()),
            ("g_na", 0.0231.into()),
            ("g_na_ttx", 400.0.into()),
            ("g_k", 0.25.into()),
            ("g_k_tea", 10.0.into()),
            ("v_k", (-70.0).into()),
            ("g_b", 0.165.into()),
            ("v_b", (-58.0).into()),
            ("g_ca", 1.5.into()),
            ("v_ca", 150.0.into()),
            ("g_ca_ca", 0.02.into()),
            ("k_beta", 15000.0.into()),
            ("beta", 0.00004.into()),
            ("rho", 0.002.into()),
            ("k_s", 50.0.into()),
            ("v0 (mV)", (-55.0).into()),
        ]
    }

    fn behavior() -> PluginBehavior {
        PluginBehavior {
            supports_start_stop: true,
            supports_restart: true,
            supports_apply: false,
            extendable_inputs: ExtendableInputs::None,
            loads_started: false,
            external_window: false,
            starts_expanded: true,
            start_requires_connected_inputs: Vec::new(),
            start_requires_connected_outputs: Vec::new(),
        }
    }
}

impl PluginRuntime for RthybridKomendantovKononenko1996NeuronRust {
    fn set_config_value(&mut self, key: &str, value: &Value) {
        let Some(v) = value.as_f64() else {
            return;
        };
        match key {
            "Burst duration (s)" => self.burst_duration_value = v,
            "i" | "I" => self.params[I] = v,
            "cm" => self.params[CM] = v,
            "g_na_v" => self.params[G_NA_V] = v,
            "v_na" => self.params[V_NA] = v,
            "g_na" => self.params[G_NA] = v,
            "g_na_ttx" => self.params[G_NA_TTX] = v,
            "g_k" => self.params[G_K] = v,
            "g_k_tea" => self.params[G_K_TEA] = v,
            "v_k" => self.params[V_K] = v,
            "g_b" => self.params[G_B] = v,
            "v_b" => self.params[V_B] = v,
            "g_ca" => self.params[G_CA] = v,
            "v_ca" => self.params[V_CA] = v,
            "g_ca_ca" => self.params[G_CA_CA] = v,
            "k_beta" => self.params[K_BETA] = v,
            "beta" => self.params[BETA] = v,
            "rho" => self.params[RHO] = v,
            "k_s" => self.params[K_S] = v,
            "v0 (mV)" => self.vars[V] = v,
            _ => {}
        }
    }

    fn set_input_value(&mut self, key: &str, v: f64) {
        match key {
            "Isyn (nA)" => self.input_syn = Self::sanitize_syn(v),
            "Burst duration (s)" => self.input_burst_duration = if v.is_finite() { v } else { 1.0 },
            _ => {}
        }
    }

    fn process_tick(&mut self, _tick: u64, period_seconds: f64) {
        if !period_seconds.is_finite() || period_seconds <= 0.0 {
            return;
        }

        self.input_syn = Self::sanitize_syn(self.input_syn);
        self.freq = 1.0 / period_seconds;
        self.refresh_burst();

        for _ in 0..self.s_points {
            let prev = self.vars;
            let params = self.params;
            let syn = self.input_syn;
            rk4_step(&mut self.vars, self.params[DT], |st, der| {
                Self::eval(st, der, &params, syn);
            });
            if self.vars.iter().any(|x| !x.is_finite()) {
                self.vars = prev;
                break;
            }
        }

        self.out_0 = self.vars[V] / 1000.0;
        self.out_1 = self.vars[V];
    }

    fn get_output_value(&self, key: &str) -> f64 {
        match key {
            "Vm (v)" => self.out_0,
            "Vm (mV)" => self.out_1,
            _ => 0.0,
        }
    }

    fn get_internal_value(&self, key: &str) -> Option<f64> {
        match key {
            "m" => Some(self.vars[M]),
            "h" => Some(self.vars[H]),
            "n" => Some(self.vars[N]),
            "ca" => Some(self.vars[CA]),
            _ => None,
        }
    }
}

rtsyn_plugin::export_plugin!(RthybridKomendantovKononenko1996NeuronRust);
