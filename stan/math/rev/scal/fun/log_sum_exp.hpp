#ifndef STAN_MATH_REV_SCAL_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_REV_SCAL_FUN_LOG_SUM_EXP_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/fun/log_sum_exp.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>

namespace stan {
namespace math {

/**
 * Returns the log sum of exponentials.
 */
inline var log_sum_exp(const var& a, const var& b) {
  operands_and_partials<var, var> res(a, b);
  double f = log_sum_exp(a.vi_->val_, b.vi_->val_);
  res.edge1_.partial_ = std::exp(a.vi_->val_ - f);
  res.edge2_.partial_ = std::exp(b.vi_->val_ - f);
  return res.build(f);
}
/**
 * Returns the log sum of exponentials.
 */
inline var log_sum_exp(const var& a, double b) {
  operands_and_partials<var, double> res(a, b);
  double f = log_sum_exp(a.vi_->val_, b);
  res.edge1_.partial_ = std::exp(a.vi_->val_ - f);
  return res.build(f);
}
/**
 * Returns the log sum of exponentials.
 */
inline var log_sum_exp(double a, const var& b) {
  operands_and_partials<double, var> res(a, b);
  double f = log_sum_exp(a, b.vi_->val_);
  res.edge2_.partial_ = std::exp(b.vi_->val_ - f);
  return res.build(f);
}

}  // namespace math
}  // namespace stan
#endif
