#ifndef STAN_MATH_REV_SCAL_FUN_LOG1P_EXP_HPP
#define STAN_MATH_REV_SCAL_FUN_LOG1P_EXP_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/fun/log1p_exp.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>

namespace stan {
namespace math {

/**
 * Return the log of 1 plus the exponential of the specified
 * variable.
 */
inline var log1p_exp(const var& a) {
  operands_and_partials<var> res(a);
  double f = log1p_exp(a.vi_->val_);
  res.edge1_.partial_ = std::exp(a.vi_->val_ - f);
  return res.build(f);
}

}  // namespace math
}  // namespace stan
#endif
