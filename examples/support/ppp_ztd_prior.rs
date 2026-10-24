//! Explicit parameter-study override; ZTD remains an estimated state.
use gnss_rust::ppp::PppConfig;

pub fn take_sigma_flag(args: &mut Vec<String>) -> Result<Option<f64>, Box<dyn std::error::Error>> {
    if args.len() < 2 || args[args.len() - 2] != "-ztd-sigma" {
        return Ok(None);
    }
    let sigma: f64 = args.pop().unwrap().parse()?;
    args.pop();
    let variance = sigma * sigma;
    if !sigma.is_finite() || sigma <= 0.0 || !variance.is_finite() || variance <= 0.0 {
        return Err("ZTD sigma and its squared variance must be finite and positive".into());
    }
    Ok(Some(sigma))
}
pub fn apply(config: &mut PppConfig, sigma: Option<f64>) {
    if let Some(sigma) = sigma {
        config.filter.initial_head_variance_m2[4] = sigma * sigma;
        eprintln!(
            "ZTD_initial_prior_sigma_m={sigma:.17e} variance_m2={:.17e} initial_mean_m={:.17e} process_noise_m2_per_s={:.17e} ZTD_still_estimated=true other_settings_unchanged=true",
            sigma * sigma,
            config.initial_zenith_delay_m,
            config.filter.head_process_noise_m2_per_s[4]
        );
    }
}
