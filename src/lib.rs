use rtsyn_plugin::prelude::*;
use serde_json::Value;
use std::mem::MaybeUninit;

#[repr(C)]
struct RthybridKomendantovKononenko1996NeuronCState {
    vars: [f64; 8],
    params: [f64; 20],
    burst_duration: f64,
    burst_duration_value: f64,
    freq: f64,
    s_points: usize,
    input_syn: f64,
    input_burst_duration: f64,
}

unsafe extern "C" {
    fn rthybrid_komendantov_kononenko_1996_neuron_c_init(
        state: *mut RthybridKomendantovKononenko1996NeuronCState,
    );
    fn rthybrid_komendantov_kononenko_1996_neuron_c_set_config(
        state: *mut RthybridKomendantovKononenko1996NeuronCState,
        key: *const u8,
        len: usize,
        value: f64,
    );
    fn rthybrid_komendantov_kononenko_1996_neuron_c_set_input(
        state: *mut RthybridKomendantovKononenko1996NeuronCState,
        key: *const u8,
        len: usize,
        value: f64,
    );
    fn rthybrid_komendantov_kononenko_1996_neuron_c_process(
        state: *mut RthybridKomendantovKononenko1996NeuronCState,
        period_seconds: f64,
    );
    fn rthybrid_komendantov_kononenko_1996_neuron_c_get_output(
        state: *const RthybridKomendantovKononenko1996NeuronCState,
        key: *const u8,
        len: usize,
    ) -> f64;
}

struct RthybridKomendantovKononenko1996NeuronC {
    state: RthybridKomendantovKononenko1996NeuronCState,
}

impl Default for RthybridKomendantovKononenko1996NeuronC {
    fn default() -> Self {
        let mut state = MaybeUninit::<RthybridKomendantovKononenko1996NeuronCState>::uninit();
        unsafe {
            rthybrid_komendantov_kononenko_1996_neuron_c_init(state.as_mut_ptr());
            Self {
                state: state.assume_init(),
            }
        }
    }
}

impl PluginDescriptor for RthybridKomendantovKononenko1996NeuronC {
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

impl PluginRuntime for RthybridKomendantovKononenko1996NeuronC {
    fn set_config_value(&mut self, key: &str, value: &Value) {
        if let Some(v) = value.as_f64() {
            unsafe {
                rthybrid_komendantov_kononenko_1996_neuron_c_set_config(
                    &mut self.state,
                    key.as_ptr(),
                    key.len(),
                    v,
                );
            }
        }
    }

    fn set_input_value(&mut self, key: &str, value: f64) {
        unsafe {
            rthybrid_komendantov_kononenko_1996_neuron_c_set_input(
                &mut self.state,
                key.as_ptr(),
                key.len(),
                value,
            );
        }
    }

    fn process_tick(&mut self, _tick: u64, period_seconds: f64) {
        unsafe { rthybrid_komendantov_kononenko_1996_neuron_c_process(&mut self.state, period_seconds) };
    }

    fn get_output_value(&self, key: &str) -> f64 {
        unsafe {
            rthybrid_komendantov_kononenko_1996_neuron_c_get_output(&self.state, key.as_ptr(), key.len())
        }
    }

    fn get_internal_value(&self, key: &str) -> Option<f64> {
        match key {
            "m" => Some(self.state.vars[3]),
            "h" => Some(self.state.vars[4]),
            "n" => Some(self.state.vars[5]),
            "ca" => Some(self.state.vars[7]),
            _ => None,
        }
    }
}

rtsyn_plugin::export_plugin!(RthybridKomendantovKononenko1996NeuronC);
