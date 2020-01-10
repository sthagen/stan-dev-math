#ifndef STAN_MATH_PRIM_MAT_FUN_INV_SQUARE_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_SQUARE_HPP

#include <stan/math/prim/mat/fun/inv.hpp>
#include <stan/math/prim/mat/fun/square.hpp>
#include <stan/math/prim/fun/inv_square.hpp>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>

namespace stan {
namespace math {

/**
 * Vectorized version of inv_square().
 *
 * @tparam T type of container
 * @param x container
 * @return 1 / the square of each value in x.
 */
template <typename T>
inline auto inv_square(const T& x) {
  return inv(square(x));
}

}  // namespace math
}  // namespace stan

#endif
