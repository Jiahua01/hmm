#include "../include/MuonScaRe.hxx"
#include "../include/Crystal.hxx"
#include "TRandom3.h"
#include <cmath>
#include <boost/math/special_functions/erf.hpp>
#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

class SeedSequence {
public:
    explicit SeedSequence(std::initializer_list<uint32_t> seeds)
        : m_seeds(seeds) {}

    template <typename Iter>
    void generate(Iter begin, Iter end) const {
        const size_t n = std::distance(begin, end);
	if (n == 0) return;

	const uint32_t mult = 0x9e3779b9;
	const uint32_t mix_const = 0x85ebca6b;

	std::vector<uint32_t> buffer(n, 0x8b8b8b8b);

	size_t s = m_seeds.size();
	size_t t = (n >= s) ? n - s : 0;

        size_t i = 0;

	for(; i < std::min(n, s); ++i) {
            buffer[i] = buffer[i] ^ (m_seeds[i] + mult * i);
	}
	for(; i < n; ++i) {
            buffer[i] = buffer[i] ^ (mult * i);
        }

	for (size_t k = 0; k < n; ++k) {
            uint32_t z = buffer[(k + n - 1) % n] ^ (buffer[k] >> 27);
	    buffer[k] = (z * mix_const) ^ (buffer[k] << 13);
        }

	std::copy(buffer.begin(), buffer.end(), begin);
    }

private:
    std::vector<uint32_t> m_seeds;
};

MuonScaRe::MuonScaRe(const std::string &json_file) {
    cset = correction::CorrectionSet::from_file(json_file);
}

double MuonScaRe::get_rndm(double eta, double phi, float nL, uint64_t event, uint32_t lumi) const{
    // obtain parameters from correctionlib
    double mean = cset->at("cb_params")->evaluate({abs(eta), nL, 0});
    double sigma = cset->at("cb_params")->evaluate({abs(eta), nL, 1});
    double n = cset->at("cb_params")->evaluate({abs(eta), nL, 2});
    double alpha = cset->at("cb_params")->evaluate({abs(eta), nL, 3});
   
    // instantiate CB and get random number following the CB
    Crystal cb(mean, sigma, alpha, n);
    int64_t phi_seed = static_cast<int64_t>((phi / M_PI) * ((1LL << 31) - 1)) & 0xFFF;
    SeedSequence seq{static_cast<uint32_t>(event), static_cast<uint32_t>(lumi), static_cast<uint32_t>(phi_seed)};
    uint32_t seed;
    seq.generate(&seed, &seed + 1);

    TRandom3 rnd(seed);
    double rndm = rnd.Rndm();
    return cb.invcdf(rndm);
}

double MuonScaRe::get_std(double pt, double eta, float nL) const {
    auto corr = cset->at("poly_params");
    double p0 = corr->evaluate({std::fabs(eta), nL, 0});
    double p1 = corr->evaluate({std::fabs(eta), nL, 1});
    double p2 = corr->evaluate({std::fabs(eta), nL, 2});
    return std::max(0.0, p0 + p1*pt + p2*pt*pt);
}

double MuonScaRe::get_k(double eta, const std::string &var) const {
    double k_data = cset->at("k_data")->evaluate({std::fabs(eta), var});
    double k_mc   = cset->at("k_mc")->evaluate({std::fabs(eta), var});
    return (k_mc < k_data) ? std::sqrt(k_data*k_data - k_mc*k_mc) : 0.0;
}

double MuonScaRe::pt_scale(bool is_data, double pt, double eta, double phi, int charge) const {
    std::string dtmc = is_data ? "data" : "mc";
    double a = cset->at("a_" + dtmc)->evaluate({eta, phi, "nom"});
    double m = cset->at("m_" + dtmc)->evaluate({eta, phi, "nom"});
    return 1.0 / (m / pt + charge * a);
}

double MuonScaRe::pt_scale_var(double pt, double eta, double phi, int charge, const std::string &updn) const {
    double stat_a = cset->at("a_mc")->evaluate({eta, phi, "stat"});
    double stat_m = cset->at("m_mc")->evaluate({eta, phi, "stat"});
    double stat_rho = cset->at("m_mc")->evaluate({eta, phi, "rho_stat"});
    double unc = pt * pt * std::sqrt(
        stat_m*stat_m / (pt*pt) +
        stat_a*stat_a +
        2 * charge * stat_rho * stat_m / pt * stat_a
    );
    return (updn == "up") ? pt + unc : pt - unc;
}

double MuonScaRe::pt_resol(double pt, double eta, double phi, float nL, uint64_t event, uint32_t lumi) const {
    double low_pt_threshold = 26;
    double rndm = get_rndm(eta, phi, nL, event, lumi);
    double std = get_std(pt, eta, nL);
    double k = get_k(eta, "nom");

    double ptc = pt * (1 + k * std * rndm);
    if (std::isnan(ptc)) ptc = pt;
    if(ptc / pt > 2 || ptc / pt < 0.1 || ptc < 0 || pt < low_pt_threshold || pt > 200){
	ptc = pt;
    }
    // TODO: Understand why for evts with pT < threshold the pt_corr is set to one
    return ptc;

}

double MuonScaRe::pt_resol_var(double pt_woresol, double pt_wresol, double eta, const std::string &updn) const {
    double k = get_k(eta, "nom");
    if (k == 0) return pt_wresol;
    double k_unc = cset->at("k_mc")->evaluate({std::fabs(eta), "stat"});
    double std_x_rndm = (pt_wresol / pt_woresol - 1) / k;
    if (updn == "up") return pt_woresol * (1 + (k + k_unc) * std_x_rndm);
    if (updn == "dn") return pt_woresol * (1 + (k - k_unc) * std_x_rndm);
    std::cerr << "ERROR: updn must be 'up' or 'dn'" << std::endl;
    return pt_wresol;
}

