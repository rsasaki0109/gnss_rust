//! Broadcast GPS/Galileo code SPP and raw Doppler. Never reads reference/IMU.
use gnss_rust::ppp::{PppMotionSeed, PppMotionSeedDiagnostic};
use gnss_rust::spp::SppConfig;
use gnss_rust::spp_qc::{SppQualityConfig, solve_gps_galileo_epoch_with_quality_control};
use gnss_rust::spp_quality::SppVarianceModel;
use gnss_rust::tides::ReceiverPositionReference;
use gnss_rust::velocity::{VelocityClockModel, VelocityEpochConfig, solve_gps_galileo_e1_epoch};
use gnss_rust::{
    GnssTime, PositionSolution, navigation::NavigationData, observation::ObservationEpoch,
};

pub const HEADER: &str = ",seed_status,seed_source,seed_frame_assertion,seed_receiver_reference,seed_x_m,seed_y_m,seed_z_m,seed_gps_clock_s,seed_position_covariance_trace_m2,seed_satellites,seed_code_rms_m,seed_gdop,seed_MAD_rejected,seed_RAIM_rejected,seed_vx_mps,seed_vy_mps,seed_vz_mps,seed_velocity_covariance_trace_m2_s2,seed_velocity_satellites,seed_velocity_residual_rms_mps,seed_velocity_chi_square_per_dof,seed_velocity_clock_systems,seed_error,seed_initialized,seed_restarted,seed_variance_trigger,seed_position_nis,seed_covariance_fallback,seed_velocity_update,seed_committed,service_status,service_x_m,service_y_m,service_z_m,service_vx_mps,service_vy_mps,service_vz_mps";
const SOURCE: &str =
    "GPS-Galileo-broadcast-native-variance-MAD-RAIM;raw-carrier-locked-Doppler-per-system-clock";
fn number(v: f64) -> String {
    format!("{v:.17e}")
}
pub struct Record {
    pub seed: Option<PppMotionSeed>,
    pub fields: Vec<String>,
}
impl Record {
    pub fn solve(e: &ObservationEpoch, nav: &NavigationData) -> Self {
        let mut fields = vec![String::new(); 23];
        fields[0] = "INVALID".into();
        fields[1] = SOURCE.into();
        fields[2] =
            "broadcast-WGS84-approximately-compatible-with-IGS20;caller-assertion-no-transform"
                .into();
        fields[3] = "InstantaneousMarker;receiver-calibration-unknown".into();
        if e.flag != 0 {
            fields[22] = "event-no-seed".into();
            return Self { seed: None, fields };
        }
        let report = match solve_gps_galileo_epoch_with_quality_control(
            e,
            nav,
            SppConfig::default(),
            SppVarianceModel::NativeSnrAtmosphere,
            SppQualityConfig::default(),
        ) {
            Ok(r) => r,
            Err(err) => {
                fields[22] = format!("SPP:{err:?}");
                return Self { seed: None, fields };
            }
        };
        let mut position = report.solution.spp.position;
        let xyz = position.position_ecef.unwrap();
        fields[0] = "SPP".into();
        for (i, v) in [xyz.x, xyz.y, xyz.z].into_iter().enumerate() {
            fields[4 + i] = number(v);
        }
        fields[7] = number(position.receiver_clock_bias);
        let q = position.position_covariance.unwrap();
        fields[8] = number(q[0][0] + q[1][1] + q[2][2]);
        fields[9] = position.num_satellites.to_string();
        fields[10] = number(report.solution.spp.residual_rms_m);
        fields[11] = number(report.evidence.summary.independent_clock_gdop);
        fields[12] = report.evidence.summary.outlier_rejections.to_string();
        fields[13] = report.evidence.summary.raim_rejections.to_string();
        let cfg = VelocityEpochConfig {
            require_carrier_lock: true,
            velocity: gnss_rust::velocity::VelocityConfig {
                clock_model: VelocityClockModel::PerSystem,
                ..Default::default()
            },
            ..Default::default()
        };
        match solve_gps_galileo_e1_epoch(e, nav, xyz, cfg) {
            Ok(r) => match r.solution {
                Ok(v) => {
                    position.velocity_ecef = Some(v.velocity_mps);
                    // Leading principal 3x3 is the marginal; never condition on clocks.
                    position.velocity_covariance = Some(std::array::from_fn(|i| {
                        std::array::from_fn(|j| v.covariance_m2s2[i][j])
                    }));
                    for (i, x) in v.velocity_mps.into_iter().enumerate() {
                        fields[14 + i] = number(x);
                    }
                    fields[17] = number(
                        v.covariance_m2s2[0][0] + v.covariance_m2s2[1][1] + v.covariance_m2s2[2][2],
                    );
                    fields[18] = v.satellites.len().to_string();
                    fields[19] = number(v.residual_rms_mps);
                    fields[20] = number(v.chi_square_per_dof);
                    fields[21] = format!("{:?}", v.clock_systems);
                }
                Err(err) => fields[22] = format!("Doppler:{err:?}"),
            },
            Err(err) => fields[22] = format!("Doppler-preparation:{err:?}"),
        }
        let seed = PppMotionSeed {
            solution: position,
            coordinate_system: "IGS20".into(),
            receiver_reference: ReceiverPositionReference::InstantaneousMarker,
            provenance: SOURCE.into(),
        };
        Self {
            seed: Some(seed),
            fields,
        }
    }
    pub fn append(
        mut self,
        d: Option<&PppMotionSeedDiagnostic>,
        ppp: Option<&PositionSolution>,
        time: GnssTime,
    ) -> Vec<String> {
        let mut diag = vec![String::new(); 7];
        if let Some(d) = d {
            assert_eq!(d.seed_time, time);
            diag[0] = d.initialized.to_string();
            diag[1] = d.restarted.to_string();
            diag[2] = d.variance_trigger.to_string();
            diag[3] = d
                .normalized_position_disagreement
                .map(number)
                .unwrap_or_default();
            diag[4] = d.seed_covariance_fallback.to_string();
            diag[5] = format!("{:?}", d.velocity_update);
            diag[6] = d.committed.to_string();
        }
        self.fields.extend(diag);
        let service = ppp.or(self.seed.as_ref().map(|s| &s.solution));
        let mut out = vec![String::new(); 7];
        out[0] = "INVALID".into();
        if let Some(s) = service {
            out[0] = if ppp.is_some() { "PPP_FLOAT" } else { "SPP" }.into();
            let p = s.position_ecef.unwrap();
            for (i, v) in [p.x, p.y, p.z].into_iter().enumerate() {
                out[1 + i] = number(v);
            }
            if let Some(v) = s.velocity_ecef {
                for (i, x) in v.into_iter().enumerate() {
                    out[4 + i] = number(x);
                }
            }
        }
        self.fields.extend(out);
        self.fields
    }
}
