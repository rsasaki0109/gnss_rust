// Whole pinned PPPProcessor on a shared decoded GPS 1C/2W observation ledger.
// No native solver source is changed. Private access is diagnostic and supports
// an explicit header/zero-clock initialization control, separately labelled.
// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
#include <bits/stdc++.h>
#include <Eigen/Dense>
#define private public
#include <libgnss++/algorithms/ppp.hpp>
#undef private
#include <libgnss++/io/rinex.hpp>
#include <libgnss++/core/constants.hpp>
using namespace libgnss;

static std::vector<std::string> fields(const std::string& line) {
    std::vector<std::string> result;std::size_t start=0;
    for(;;) {const auto end=line.find(',',start);result.push_back(line.substr(start,end==std::string::npos?end:end-start));if(end==std::string::npos)break;start=end+1;}
    return result;
}
static double number(const std::string& value) {
    std::size_t end=0;const double x=std::stod(value,&end);
    if(end!=value.size()||!std::isfinite(x))throw std::runtime_error("nonfinite or malformed ledger number");
    return x;
}
static std::string text(double x) {
    std::ostringstream out;out << std::scientific << std::setprecision(17) << x;return out.str();
}
static std::string optional(double x,bool present) {return present?text(x):"";}
static const char* status_name(SolutionStatus s) {
    switch(s) {
        case SolutionStatus::NONE:return "NONE";case SolutionStatus::SPP:return "SPP";
        case SolutionStatus::DGPS:return "DGPS";case SolutionStatus::FLOAT:return "FLOAT";
        case SolutionStatus::FIXED:return "FIX";case SolutionStatus::PPP_FLOAT:return "PPP_FLOAT";
        case SolutionStatus::PPP_FIXED:return "PPP_FIX";case SolutionStatus::PROPAGATED:return "PROPAGATED";
    }
    throw std::runtime_error("unknown native status");
}
static void row(std::ostream& out,const std::vector<std::string>& values) {
    for(std::size_t i=0;i<values.size();++i) {if(i)out << ',';out << values[i];}out << '\n';
}

int main(int argc,char**argv) {
    try {
        if(argc!=11&&argc!=12) {std::cerr << "usage: native_ppc_ppp LEDGER SP3 CLK COMPAT_CODE_BIA NAV GPS_ANTEX OUTPUT_CSV READBACK_CSV BIAS_READBACK_CSV builtin|header-zero [MAX_EPOCHS]\n";return 2;}
        const bool header_seed=std::string(argv[10])=="header-zero";
        if(!header_seed&&std::string(argv[10])!="builtin")throw std::runtime_error("unknown seed profile");
        const std::size_t maximum=argc==12?std::stoull(argv[11]):0;
        for(int i:{7,8,9})if(std::filesystem::exists(argv[i]))throw std::runtime_error("output already exists");
        std::ifstream input(argv[1]);if(!input)throw std::runtime_error("cannot open ledger");
        std::string line;if(!std::getline(input,line)||line!="kind,week,tow,satellite,tracking,code_m,phase_cycles,doppler_hz,snr_db_hz,lli,ssi,flag,receiver_offset_s,header_x_m,header_y_m,header_z_m")throw std::runtime_error("unknown ledger schema");
        std::ofstream out(argv[7]),readback(argv[8]),bias_out(argv[9]);if(!out||!readback||!bias_out)throw std::runtime_error("cannot create outputs");readback << line << '\n';
        std::cout.rdbuf(std::cerr.rdbuf());
        PPPProcessor::PPPConfig config;config.kinematic_mode=true;config.low_dynamics_mode=false;config.use_dynamics_model=false;config.use_pva_dynamics=false;
        config.use_ionosphere_free=true;config.estimate_ionosphere=false;config.estimate_troposphere=true;
        config.enable_ambiguity_resolution=false;config.apply_solid_earth_tides=false;config.apply_ocean_loading=false;
        config.orbit_file_path=argv[2];config.clock_file_path=argv[3];config.dcb_file_path=argv[4];config.antex_file_path=argv[6];config.receiver_antenna_type="";
        PPPProcessor processor(config);ProcessorConfig general;general.mode=PositioningMode::PPP;general.elevation_mask=10.0;general.orbit_file_path=argv[2];general.clock_file_path=argv[3];
        if(!processor.initialize(general)||!processor.precise_products_loaded_||!processor.dcb_products_loaded_||!processor.satellite_antex_loaded_)throw std::runtime_error("required native precise/bias/ANTEX input unavailable");
        if(processor.dcb_products_.entries.size()!=64)throw std::runtime_error("expected 64 GPS code-only OSB rows");
        bias_out << "satellite,observation,unit,value,sigma,metres\n";
        for(const auto& b:processor.dcb_products_.entries) {
            if(!b.valid||b.bias_type!="OSB"||b.satellite.system!=GNSSSystem::GPS||b.observation_1!=b.observation_2)throw std::runtime_error("unexpected native bias view");
            double metres;if(b.unit=="ns")metres=b.bias*constants::SPEED_OF_LIGHT*1e-9;else if(b.unit=="m")metres=b.bias;else throw std::runtime_error("unsupported native bias unit");
            row(bias_out,{b.satellite.toString(),b.observation_1,b.unit,text(b.bias),text(b.sigma),text(metres)});
        }
        io::RINEXReader navigation_reader;NavigationData navigation;
        if(!navigation_reader.open(argv[5])||!navigation_reader.readNavigationData(navigation))throw std::runtime_error("cannot read navigation");
        out << "week,tow,status,native_status,native_is_valid,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,clock_s,clock_native_raw,clock_native_unit,position_covariance_trace_m2,reported_residual_rms_m,satellites,solution_week,solution_tow,filter_initialized,filter_zenith_m,filter_vx_mps,filter_vy_mps,filter_vz_mps,filter_clock_m,filter_states,filter_ambiguity_states,DCB_applied_count,DCB_applied_abs_sum_m,error\n";
        std::cerr << "whole_native_PPP_processEpoch=true upstream_source_patched=false input=shared_predecoded_GPS_1C_2W phase_OSBs_applied=false receiver_calibration_known=false satellite_PCV_applied=false solid_earth_tide_applied=false AR_enabled=false nominal_code_sigma_m=" << config.pseudorange_sigma << " nominal_phase_sigma_m=" << config.carrier_phase_sigma << " native_variance_model=" << config.use_rtklib_measurement_variance << " native_outlier_threshold=" << config.outlier_threshold << " head_seed_profile=" << argv[10] << " native_SPP_precise_loaded=" << processor.spp_processor_.precise_products_loaded_ << " native_SPP_DCB_loaded=" << processor.spp_processor_.dcb_products_loaded_ << '\n';
        ObservationData epoch;bool pending=false,done=false;std::size_t count=0,groups=0;GNSSTime previous;bool have_previous=false;
        auto finish=[&]() {
            if(!pending)return;
            if(have_previous&&!(epoch.time-previous>0.0))throw std::runtime_error("ledger epochs do not increase");
            previous=epoch.time;have_previous=true;
            if(count==0&&header_seed) {
                PositionSolution seed;seed.status=SolutionStatus::SPP;seed.time=epoch.time;seed.position_ecef=epoch.receiver_position;seed.receiver_clock_bias=0.0;seed.num_satellites=4;
                if(!processor.initializeFilter(epoch,navigation,&seed))throw std::runtime_error("header control initialization failed");
                processor.filter_initialized_=true;processor.convergence_start_time_=epoch.time;
            }
            const auto solution=processor.processEpoch(epoch,navigation);const bool valid=solution.isValid()&&solution.position_ecef.allFinite();const bool spp=solution.status==SolutionStatus::SPP;
            std::vector<std::string> r(30);r[0]=std::to_string(epoch.time.week);r[1]=text(epoch.time.tow);r[2]=valid?status_name(solution.status):"INVALID";r[3]=status_name(solution.status);r[4]=solution.isValid()?"1":"0";
            if(valid) {
                for(int i=0;i<3;++i)r[5+i]=text(solution.position_ecef[i]);
                if(solution.has_velocity)for(int i=0;i<3;++i)r[8+i]=optional(solution.velocity_ecef[i],std::isfinite(solution.velocity_ecef[i]));
                r[11]=text(spp?solution.receiver_clock_bias/constants::SPEED_OF_LIGHT:solution.receiver_clock_bias);r[12]=text(solution.receiver_clock_bias);r[13]=spp?"metres_from_native_SPP":"seconds_from_native_PPP";
                r[14]=optional(solution.position_covariance.trace(),solution.position_covariance.allFinite());r[15]=optional(solution.residual_rms,std::isfinite(solution.residual_rms));r[16]=std::to_string(solution.num_satellites);
            }
            r[17]=std::to_string(solution.time.week);r[18]=text(solution.time.tow);r[19]=processor.filter_initialized_?"1":"0";
            if(processor.filter_initialized_) {
                const auto& f=processor.filter_state_;r[20]=text(f.state[f.trop_index]);for(int i=0;i<3;++i)r[21+i]=text(f.state[f.vel_index+i]);r[24]=text(f.state[f.clock_index]);r[25]=std::to_string(f.total_states);r[26]=std::to_string(f.ambiguity_indices.size());
            }
            r[27]=std::to_string(processor.getLastAppliedDcbCorrections());r[28]=text(processor.getLastAppliedDcbMeters());
            if(!valid)r[29]=solution.isValid()?"native_nonfinite_position":"native_no_valid_solution";
            row(out,r);++count;pending=false;if(maximum&&count>=maximum)done=true;
        };
        while(std::getline(input,line)) {
            const auto f=fields(line);if(f.size()!=16)throw std::runtime_error("wrong ledger field count");
            if(f[0]=="epoch") {
                finish();if(done)break;
                const int week=std::stoi(f[1]);const double tow=number(f[2]);epoch.clear();epoch.time=GNSSTime(week,tow);
                if(f[11]!="0"||(!f[12].empty()&&number(f[12])!=0.0))throw std::runtime_error("native benchmark supports normal epochs with no unapplied clock offset only");
                epoch.receiver_position=Vector3d(number(f[13]),number(f[14]),number(f[15]));epoch.receiver_clock_bias=0.0;epoch.setRinexFrequencySlot(GNSSSystem::GPS,0,"1C");epoch.setRinexFrequencySlot(GNSSSystem::GPS,1,"2W");
                std::vector<std::string> b=f;b[1]=std::to_string(week);b[2]=text(tow);for(int i=0;i<3;++i)b[13+i]=text(epoch.receiver_position[i]);row(readback,b);pending=true;
            } else if(f[0]=="observation") {
                if(!pending||std::stoi(f[1])!=epoch.time.week||number(f[2])!=epoch.time.tow||f[3].size()!=3||f[3][0]!='G'||(f[4]!="1C"&&f[4]!="2W"))throw std::runtime_error("observation outside the selected epoch/system/pair");
                const int prn=std::stoi(f[3].substr(1));if(prn<1||prn>32)throw std::runtime_error("invalid GPS PRN");
                Observation o(SatelliteId(GNSSSystem::GPS,prn),f[4]=="1C"?SignalType::GPS_L1CA:SignalType::GPS_L2C);
                o.has_pseudorange=!f[5].empty();o.pseudorange=o.has_pseudorange?number(f[5]):0.0;o.has_carrier_phase=!f[6].empty();o.carrier_phase=o.has_carrier_phase?number(f[6]):0.0;o.has_doppler=!f[7].empty();o.doppler=o.has_doppler?number(f[7]):0.0;o.snr=f[8].empty()?0.0:number(f[8]);
                const int lli=f[9].empty()?0:std::stoi(f[9]);const int strength=f[10].empty()?0:std::stoi(f[10]);if(lli<0||lli>7||strength<0||strength>9)throw std::runtime_error("invalid raw quality flag");o.lli=static_cast<uint8_t>(lli);o.loss_of_lock=(lli&1)!=0;o.signal_strength=strength;o.valid=true;
                if(o.has_pseudorange)o.pseudorange_observation_type="C"+f[4];
                if(o.has_carrier_phase)o.carrier_phase_observation_type="L"+f[4];
                epoch.observations.push_back(o);epoch.addRinexTrackingObservation(f[4],o);
                std::vector<std::string> b=f;b[1]=std::to_string(epoch.time.week);b[2]=text(epoch.time.tow);b[3]=o.satellite.toString();b[5]=optional(o.pseudorange,o.has_pseudorange);b[6]=optional(o.carrier_phase,o.has_carrier_phase);b[7]=optional(o.doppler,o.has_doppler);b[8]=optional(o.snr,!f[8].empty());b[9]=f[9].empty()?"":std::to_string(o.lli);b[10]=f[10].empty()?"":std::to_string(o.signal_strength);row(readback,b);++groups;
            } else throw std::runtime_error("unknown ledger row kind");
        }
        if(!done)finish();
        out.flush();readback.flush();bias_out.flush();if(!out||!readback||!bias_out)throw std::runtime_error("output failure");
        std::cerr << "native_processEpoch_calls=" << count << " imported_groups=" << groups << " native_RINEX_observation_pipeline_executed=false shared_ledger_readback=true truth_used_by_solver=false practical_parity_accepted=false\n";
        return 0;
    } catch(const std::exception& e) {std::cerr << "native_ppc_ppp: " << e.what() << '\n';return 1;}
}
