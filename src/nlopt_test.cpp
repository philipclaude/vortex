#include "tester.h"
#include "log.h"

#include <nlopt.hpp>
#include <fmt/format.h>

UT_TEST_SUITE(nlopt_test_suite)

struct nlopt_data
{
    double k;
    double p;
    int iter{0};
    std::vector<double> target_mass;
};

double energy(unsigned n, const double *w, double *de_dw, void *data0)
{
    nlopt_data &data = *static_cast<nlopt_data *>(data0);
    data.iter++;

    LOG << fmt::format("iteration {}, nlopt guessed x = {}, y = {}", data.iter, w[0], w[1]);
    // calculate power diagram
    if (de_dw)
    {
        de_dw[0] = data.k * std::pow(w[0], data.p - 1);
        de_dw[1] = data.k * std::pow(w[1], data.p - 1);
    }

    return data.k * (std::pow(w[0], data.p) + std::pow(w[1], data.p));
}

UT_TEST_CASE(test1)
{
    int n = 2;
    nlopt::opt opt(nlopt::LD_LBFGS, n);

    nlopt_data data = {10, 2, 0, {}};
    opt.set_min_objective(&energy, static_cast<void *>(&data));

    // set some optimization parameters
    opt.set_xtol_rel(1e-12);
    opt.set_ftol_rel(1e-12);
    opt.set_maxeval(100);

    // set the lower and upper bounds on the weights
    std::vector<double> lower_bound(n, -10.0);
    std::vector<double> upper_bound(n, 10.0);
    opt.set_lower_bounds(lower_bound);
    opt.set_upper_bounds(upper_bound);

    double f_opt;
    std::vector<double> x = {1, 1};
    auto result = opt.optimize(x, f_opt);
    LOG << fmt::format("x = {}, y = {}, f = {}", x[0], x[1], f_opt);
}
UT_TEST_CASE_END(test1)

UT_TEST_SUITE_END(nlopt_test_suite)