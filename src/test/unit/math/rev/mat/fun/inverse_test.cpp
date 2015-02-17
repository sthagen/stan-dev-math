#include <stan/math/prim/mat/fun/inverse.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/prim/mat/fun/sum.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/rev/mat/fun/multiply.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/scal/fun/operator_addition.hpp>
#include <stan/math/rev/scal/fun/operator_divide_equal.hpp>
#include <stan/math/rev/scal/fun/operator_division.hpp>
#include <stan/math/rev/scal/fun/operator_equal.hpp>
#include <stan/math/rev/scal/fun/operator_greater_than.hpp>
#include <stan/math/rev/scal/fun/operator_greater_than_or_equal.hpp>
#include <stan/math/rev/scal/fun/operator_less_than.hpp>
#include <stan/math/rev/scal/fun/operator_less_than_or_equal.hpp>
#include <stan/math/rev/scal/fun/operator_minus_equal.hpp>
#include <stan/math/rev/scal/fun/operator_multiplication.hpp>
#include <stan/math/rev/scal/fun/operator_multiply_equal.hpp>
#include <stan/math/rev/scal/fun/operator_not_equal.hpp>
#include <stan/math/rev/scal/fun/operator_plus_equal.hpp>
#include <stan/math/rev/scal/fun/operator_subtraction.hpp>
#include <stan/math/rev/scal/fun/operator_unary_decrement.hpp>
#include <stan/math/rev/scal/fun/operator_unary_increment.hpp>
#include <stan/math/rev/scal/fun/operator_unary_negative.hpp>
#include <stan/math/rev/scal/fun/operator_unary_not.hpp>
#include <stan/math/rev/scal/fun/operator_unary_plus.hpp>
#include <stan/math/rev/mat/fun/sum.hpp>
#include <stan/math/rev/scal/fun/abs.hpp>

TEST(AgradRevMatrix,inverse_val) {
  using stan::math::inverse;
  using stan::agrad::matrix_v;

  matrix_v a(2,2);
  a << 2.0, 3.0, 
    5.0, 7.0;

  matrix_v a_inv = inverse(a);

  matrix_v I = multiply(a,a_inv);

  EXPECT_NEAR(1.0,I(0,0).val(),1.0E-12);
  EXPECT_NEAR(0.0,I(0,1).val(),1.0E-12);
  EXPECT_NEAR(0.0,I(1,0).val(),1.0E-12);
  EXPECT_NEAR(1.0,I(1,1).val(),1.0e-12);

  EXPECT_THROW(inverse(matrix_v(2,3)), std::invalid_argument);
}
TEST(AgradRevMatrix,inverse_grad) {
  using stan::math::inverse;
  using stan::agrad::matrix_v;
  
  for (size_t k = 0; k < 2; ++k) {
    for (size_t l = 0; l < 2; ++l) {

      matrix_v ad(2,2);
      ad << 2.0, 3.0, 
        5.0, 7.0;

      AVEC x = createAVEC(ad(0,0),ad(0,1),ad(1,0),ad(1,1));

      matrix_v ad_inv = inverse(ad);

      // int k = 0;
      // int l = 1;
      VEC g;
      ad_inv(k,l).grad(x,g);

      int idx = 0;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          EXPECT_FLOAT_EQ(-ad_inv(k,i).val() * ad_inv(j,l).val(), g[idx]);
          ++idx;
        }
      }
    }
  }
}
TEST(AgradRevMatrix,inverse_inverse_sum) {
  using stan::math::sum;
  using stan::math::inverse;
  using stan::agrad::matrix_v;

  matrix_v a(4,4);
  a << 2.0, 3.0, 4.0, 5.0, 
    9.0, -1.0, 2.0, 2.0,
    4.0, 3.0, 7.0, -1.0,
    0.0, 1.0, 19.0, 112.0;

  AVEC x;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      x.push_back(a(i,j));

  AVAR a_inv_inv_sum = sum(inverse(inverse(a)));

  VEC g = cgradvec(a_inv_inv_sum,x);

  for (size_t k = 0; k < x.size(); ++k)
    EXPECT_FLOAT_EQ(1.0,g[k]);
}
