#ifndef STAN_MATH_PRIM_SCAL_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LOG_SUM_EXP_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/exp.hpp>
#include <stan/math/prim/scal/fun/log1p_exp.hpp>
//#include <boost/math/tools/promotion.hpp>
#include <limits>

namespace stan {
namespace math {

/**
 * Calculates the log sum of exponetials without overflow.
 *
 * \f$\log (\exp(a) + \exp(b)) = m + \log(\exp(a-m) + \exp(b-m))\f$,
 *
 * where \f$m = max(a, b)\f$.
 *
 *
   \f[
   \mbox{log\_sum\_exp}(x, y) =
   \begin{cases}
     \ln(\exp(x)+\exp(y)) & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_sum\_exp}(x, y)}{\partial x} =
   \begin{cases}
     \frac{\exp(x)}{\exp(x)+\exp(y)} & \mbox{if } -\infty\leq x, y \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_sum\_exp}(x, y)}{\partial y} =
   \begin{cases}
     \frac{\exp(y)}{\exp(x)+\exp(y)} & \mbox{if } -\infty\leq x, y \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a the first variable
 * @param b the second variable
 */
template <typename T1, typename T2>
inline typename return_type<T1, T2>::type log_sum_exp(const T1& a,
                                                      const T2& b) {
  typedef typename stan::partials_return_type<T1, T2>::type T_partials_return;

  using std::exp;

  const T_partials_return a_dbl = value_of(a);
  const T_partials_return b_dbl = value_of(b);

  const T_partials_return fab = a_dbl > b_dbl
                                    ? a_dbl + log1p_exp(b_dbl - a_dbl)
                                    : b_dbl + log1p_exp(a_dbl - b_dbl);

  operands_and_partials_2<T1, T2> ops_partials(a, b);

  if (!is_constant_struct<T1>::value)
    ops_partials.edge1_.partials_[0] = exp(a_dbl - fab);

  if (!is_constant_struct<T2>::value)
    ops_partials.edge2_.partials_[0] = exp(b_dbl - fab);

  return ops_partials.build(fab);
}

}  // namespace math
}  // namespace stan

#endif
