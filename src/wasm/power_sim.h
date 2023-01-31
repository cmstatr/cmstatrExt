#ifndef POWER_SIM_H
#define POWER_SIM_H

bool power_sim_dual_normal(
    const int n_qual, const int m_equiv,
    const int replicates,
    const double mu_pop_qual, const double sd_pop_qual,
    std::vector<double> mu_pop_equiv, std::vector<double> sd_pop_equiv,
    const double k1, const double k2,
    std::vector<double> &reject_rate);

#endif