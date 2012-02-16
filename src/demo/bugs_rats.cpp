// Code generated by Stan version alpha.0

#include <cmath>
#include <vector>
#include <sstream>
#include <Eigen/Dense>
#include <stan/agrad/agrad.hpp>
#include <stan/agrad/special_functions.hpp>
#include <stan/agrad/matrix.hpp>
#include <stan/io/dump.hpp>
#include <stan/io/reader.hpp>
#include <stan/io/csv_writer.hpp>
#include <stan/math/matrix.hpp>
#include <stan/math/special_functions.hpp>
#include <stan/mcmc/hmc.hpp>
#include <stan/mcmc/sampler.hpp>
#include <stan/model/prob_grad_ad.hpp>
#include <stan/prob/distributions.hpp>

namespace test_model {

using std::vector;
using std::string;
using std::stringstream;
using stan::agrad::var;
using stan::mcmc::prob_grad_ad;
using stan::io::dump;
using std::istream;

typedef Eigen::Matrix<double,1,Eigen::Dynamic> vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;
typedef Eigen::Matrix<stan::agrad::var,1,Eigen::Dynamic> vector_v;
typedef Eigen::Matrix<stan::agrad::var,Eigen::Dynamic,1> row_vector_v;
typedef Eigen::Matrix<stan::agrad::var,Eigen::Dynamic,Eigen::Dynamic> matrix_v;

class test_model : public prob_grad_ad {
private:
    int N;
    int T;
    vector<double> x;
    double xbar;
public:
    test_model(std::istream& in)
        : prob_grad_ad::prob_grad_ad(0) {
        dump dump(in);
        std::vector<int> vals_i;
        std::vector<double> vals_r;
        assert(dump.contains_i("N"));
        vals_i = dump.var_i_vals("N");
        pos = 0;
        N = vals_i[pos++];
        assert(dump.contains_i("T"));
        vals_i = dump.var_i_vals("T");
        pos = 0;
        T = vals_i[pos++];
        assert(dump.contains_r("x"));
        vals_r = dump.var_r_vals("x");
        pos = 0;
        unsigned int x_limit__0 = T;
        for (unsigned int i_0 = 0; i_0 < x_limit__0; ++i_0) {
            x[i_0] = vals_r[pos++];
        }
        assert(dump.contains_r("xbar"));
        vals_r = dump.var_r_vals("xbar");
        pos = 0;
        xbar = vals_r[pos++];
    } // dump ctor

    var log_prob(vector<var>& params_r,
                 vector<int>& params_i) {

        var lp(0.0);
        stan::io::reader<var> in(params_r,params_i);

        // model parameters
        vector<vector<var> > mu;
        unsigned int dim_mu_0 = N;
        mu.resize(dim_mu_0);
        for (unsigned int k_0 = 0; k_0 < dim_mu_0; ++k_0) {
            unsigned int dim_mu_1 = T;
            for (unsigned int k_1 = 0; k_1 < dim_mu_1; ++k_1) {
                mu[k_0].push_back(in.scalar_constrain(lp));
            }
        }
        var tau_c = in.scalar_lb_constrain(0,lp);
        vector<var> alpha;
        unsigned int dim_alpha_0 = N;
        for (unsigned int k_0 = 0; k_0 < dim_alpha_0; ++k_0) {
            alpha.push_back(in.scalar_constrain(lp));
        }
        vector<var> beta;
        unsigned int dim_beta_0 = N;
        for (unsigned int k_0 = 0; k_0 < dim_beta_0; ++k_0) {
            beta.push_back(in.scalar_constrain(lp));
        }
        var alpha_c = in.scalar_constrain(lp);
        var alpha_tau = in.scalar_constrain(lp);
        var beta_c = in.scalar_constrain(lp);
        var beta_tau = in.scalar_constrain(lp);
        vector<vector<var> > Y;
        unsigned int dim_Y_0 = N;
        Y.resize(dim_Y_0);
        for (unsigned int k_0 = 0; k_0 < dim_Y_0; ++k_0) {
            unsigned int dim_Y_1 = T;
            for (unsigned int k_1 = 0; k_1 < dim_Y_1; ++k_1) {
                Y[k_0].push_back(in.scalar_lb_constrain(0,lp));
            }
        }
        var sigma = in.scalar_constrain(lp);
        var alpha0 = in.scalar_constrain(lp);

        // derived variables

        // model body
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= T; ++j) {
                lp += stan::prob::normal_log(Y[i - 1][j - 1], mu[i - 1][j - 1], tau_c);
                mu[i - 1][j - 1] = (alpha[i - 1] + (beta[i - 1] * (x[j - 1] - xbar)));
            }
            lp += stan::prob::normal_log(alpha[i - 1], alpha_c, alpha_tau);
            lp += stan::prob::normal_log(beta[i - 1], beta_c, beta_tau);
        }
        lp += stan::prob::gamma_log(tau_c, 0.001, 0.001);
        sigma = (1 / sqrt(tau_c));
        lp += stan::prob::normal_log(alpha_c, 0, 1e-06);
        lp += stan::prob::gamma_log(alpha_tau, 0.001, 0.001);
        lp += stan::prob::normal_log(beta_c, 0, 1e-06);
        lp += stan::prob::gamma_log(beta_tau, 0.001, 0.001);
        alpha0 = (alpha_c - (xbar * beta_c));

        return lp;

    } // log_prob()

    void write_csv(const std::vector<double>& params_r__,
                   const std::vector<int>& params_i__,
                   std::ostream& o__) {
        stan::io::reader<double> in__(params_r__,params_i__);
        stan::io::csv_writer writer__(o__);
        //---------------------------------
        vector<vector<var> > mu;
        unsigned int dim_mu_0 = N;
        mu.resize(dim_mu_0);
        for (unsigned int k_0 = 0; k_0 < dim_mu_0; ++k_0) {
            unsigned int dim_mu_1 = T;
            for (unsigned int k_1 = 0; k_1 < dim_mu_1; ++k_1) {
                mu[k_0].push_back(in.scalar_constrain(lp));
            }
        }
        var tau_c = in.scalar_lb_constrain(0,lp);
        vector<var> alpha;
        unsigned int dim_alpha_0 = N;
        for (unsigned int k_0 = 0; k_0 < dim_alpha_0; ++k_0) {
            alpha.push_back(in.scalar_constrain(lp));
        }
        vector<var> beta;
        unsigned int dim_beta_0 = N;
        for (unsigned int k_0 = 0; k_0 < dim_beta_0; ++k_0) {
            beta.push_back(in.scalar_constrain(lp));
        }
        var alpha_c = in.scalar_constrain(lp);
        var alpha_tau = in.scalar_constrain(lp);
        var beta_c = in.scalar_constrain(lp);
        var beta_tau = in.scalar_constrain(lp);
        vector<vector<var> > Y;
        unsigned int dim_Y_0 = N;
        Y.resize(dim_Y_0);
        for (unsigned int k_0 = 0; k_0 < dim_Y_0; ++k_0) {
            unsigned int dim_Y_1 = T;
            for (unsigned int k_1 = 0; k_1 < dim_Y_1; ++k_1) {
                Y[k_0].push_back(in.scalar_lb_constrain(0,lp));
            }
        }
        var sigma = in.scalar_constrain(lp);
        var alpha0 = in.scalar_constrain(lp);
        //---------------------------------
        unsigned int bound_mu_0 = N;
        unsigned int bound_mu_1 = T;
        for (unsigned int i_mu_0 = 0; i_mu_0 < bound_mu_0; ++i_mu_0) {
            for (unsigned int i_mu_1 = 0; i_mu_1 < bound_mu_1; ++i_mu_1) {
                writer__.write(in__.scalar_constrain());
            }
        }
        writer__.write(in__.scalar_lb_constrain(0));
        unsigned int bound_alpha_0 = N;
        for (unsigned int i_alpha_0 = 0; i_alpha_0 < bound_alpha_0; ++i_alpha_0) {
            writer__.write(in__.scalar_constrain());
        }
        unsigned int bound_beta_0 = N;
        for (unsigned int i_beta_0 = 0; i_beta_0 < bound_beta_0; ++i_beta_0) {
            writer__.write(in__.scalar_constrain());
        }
        writer__.write(in__.scalar_constrain());
        writer__.write(in__.scalar_constrain());
        writer__.write(in__.scalar_constrain());
        writer__.write(in__.scalar_constrain());
        unsigned int bound_Y_0 = N;
        unsigned int bound_Y_1 = T;
        for (unsigned int i_Y_0 = 0; i_Y_0 < bound_Y_0; ++i_Y_0) {
            for (unsigned int i_Y_1 = 0; i_Y_1 < bound_Y_1; ++i_Y_1) {
                writer__.write(in__.scalar_lb_constrain(0));
            }
        }
        writer__.write(in__.scalar_constrain());
        writer__.write(in__.scalar_constrain());
    }

    void write_csv_header(std::ostream& o__) {
        stan::io::csv_writer writer__(o__);
        unsigned int bound_mu_0 = N;
        unsigned int bound_mu_1 = T;
        for (unsigned int i_mu_0 = 0; i_mu_0 < bound_mu_0; ++i_mu_0) {
            for (unsigned int i_mu_1 = 0; i_mu_1 < bound_mu_1; ++i_mu_1) {
                stringstream ss_mu__;
                ss_mu__ << "mu[" << i_mu_0 << ',' << i_mu_1 << ',' << ']';
                writer__.write(ss_mu__.str());
            }
        }
        stringstream ss_tau_c__;
        ss_tau_c__ << "tau_c[" << ']';
        writer__.write(ss_tau_c__.str());
        unsigned int bound_alpha_0 = N;
        for (unsigned int i_alpha_0 = 0; i_alpha_0 < bound_alpha_0; ++i_alpha_0) {
            stringstream ss_alpha__;
            ss_alpha__ << "alpha[" << i_alpha_0 << ',' << ']';
            writer__.write(ss_alpha__.str());
        }
        unsigned int bound_beta_0 = N;
        for (unsigned int i_beta_0 = 0; i_beta_0 < bound_beta_0; ++i_beta_0) {
            stringstream ss_beta__;
            ss_beta__ << "beta[" << i_beta_0 << ',' << ']';
            writer__.write(ss_beta__.str());
        }
        stringstream ss_alpha_c__;
        ss_alpha_c__ << "alpha_c[" << ']';
        writer__.write(ss_alpha_c__.str());
        stringstream ss_alpha_tau__;
        ss_alpha_tau__ << "alpha_tau[" << ']';
        writer__.write(ss_alpha_tau__.str());
        stringstream ss_beta_c__;
        ss_beta_c__ << "beta_c[" << ']';
        writer__.write(ss_beta_c__.str());
        stringstream ss_beta_tau__;
        ss_beta_tau__ << "beta_tau[" << ']';
        writer__.write(ss_beta_tau__.str());
        unsigned int bound_Y_0 = N;
        unsigned int bound_Y_1 = T;
        for (unsigned int i_Y_0 = 0; i_Y_0 < bound_Y_0; ++i_Y_0) {
            for (unsigned int i_Y_1 = 0; i_Y_1 < bound_Y_1; ++i_Y_1) {
                stringstream ss_Y__;
                ss_Y__ << "Y[" << i_Y_0 << ',' << i_Y_1 << ',' << ']';
                writer__.write(ss_Y__.str());
            }
        }
        stringstream ss_sigma__;
        ss_sigma__ << "sigma[" << ']';
        writer__.write(ss_sigma__.str());
        stringstream ss_alpha0__;
        ss_alpha0__ << "alpha0[" << ']';
        writer__.write(ss_alpha0__.str());
    }

}; // model

} // namespace

