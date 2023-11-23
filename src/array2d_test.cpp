#include "array2d.h"

#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(array2d_test_suite)

UT_TEST_CASE(test_rectangular) {
  array2d<int> a(3);

  UT_ASSERT_EQUALS(a.n(), 0);

  int x[3] = {1, 2, 3};
  a.add(x);
  UT_ASSERT_EQUALS(a.n(), 1);

  UT_ASSERT_EQUALS(a(0, 0), 1);
  UT_ASSERT_EQUALS(a(0, 1), 2);
  UT_ASSERT_EQUALS(a(0, 2), 3);

  UT_ASSERT_EQUALS(a[0][0], 1);
  UT_ASSERT_EQUALS(a[0][1], 2);
  UT_ASSERT_EQUALS(a[0][2], 3);

  const int* a0 = a[0];
  UT_ASSERT_EQUALS(a0[0], 1);
  UT_ASSERT_EQUALS(a0[1], 2);
  UT_ASSERT_EQUALS(a0[2], 3);

  a(0, 0) = 3;
  a(0, 1) = 2;
  a(0, 2) = 1;
  UT_ASSERT_EQUALS(a(0, 0), 3);
  UT_ASSERT_EQUALS(a(0, 1), 2);
  UT_ASSERT_EQUALS(a(0, 2), 1);

  a.clear();
  UT_CATCH_EXCEPTION(a.length(0));
  a.set_stride(4);
  UT_ASSERT_EQUALS(a.n(), 0);
  int y[4] = {4, 5, 6, 7};
  a.add(y);
  UT_ASSERT_EQUALS(a.n(), 1);
  UT_ASSERT_EQUALS(a.length(0), 4);

  UT_ASSERT_EQUALS(a(0, 0), 4);
  UT_ASSERT_EQUALS(a(0, 1), 5);
  UT_ASSERT_EQUALS(a(0, 2), 6);
  UT_ASSERT_EQUALS(a(0, 3), 7);

  array2d<int> b(2, 10, 7);
  UT_ASSERT_EQUALS(b.n(), 10);
  for (int k = 0; k < 10; k++) {
    UT_ASSERT_EQUALS(b.length(k), 2);
    for (int d = 0; d < 2; d++) UT_ASSERT_EQUALS(b(k, d), 7);
  }
  a.print();
  b.print();

  b.clear(false);
  int z[2] = {8, 9};
  b.add(z);
  UT_ASSERT_EQUALS(b.length(0), 2);

  // test removal

  int u[4] = {8, 9, 10, 11};
  int v[4] = {-1, 8, 0, 2};
  a.add(u);
  a.add(v);
  LOG << "a (before)";
  a.print();
  UT_ASSERT_EQUALS(a.n(), 3);
  a.remove(1);
  UT_ASSERT_EQUALS(a.n(), 2);
  UT_ASSERT_EQUALS(a.data().size(), 8);
  LOG << "a (after)";
  a.print();
}
UT_TEST_CASE_END(test_rectangular)

UT_TEST_CASE(test_jagged) {
  array2d<int> a(-1);

  UT_ASSERT_EQUALS(a.n(), 0);

  int x[3] = {1, 2, 3};
  UT_CATCH_EXCEPTION(a.add(x));
  a.add(x, 3);

  UT_ASSERT_EQUALS(a.n(), 1);
  UT_ASSERT_EQUALS(a.length(0), 3);

  UT_ASSERT_EQUALS(a(0, 0), 1);
  UT_ASSERT_EQUALS(a(0, 1), 2);
  UT_ASSERT_EQUALS(a(0, 2), 3);

  UT_ASSERT_EQUALS(a[0][0], 1);
  UT_ASSERT_EQUALS(a[0][1], 2);
  UT_ASSERT_EQUALS(a[0][2], 3);

  const int* a0 = a[0];
  UT_ASSERT_EQUALS(a0[0], 1);
  UT_ASSERT_EQUALS(a0[1], 2);
  UT_ASSERT_EQUALS(a0[2], 3);

  a(0, 0) = 3;
  a(0, 1) = 2;
  a(0, 2) = 1;
  UT_ASSERT_EQUALS(a(0, 0), 3);
  UT_ASSERT_EQUALS(a(0, 1), 2);
  UT_ASSERT_EQUALS(a(0, 2), 1);

  int y[4] = {4, 5, 6, 7};
  a.add(y, 4);
  int z[2] = {8, 9};
  a.add(z, 2);

  UT_ASSERT_EQUALS(a.n(), 3);
  LOG << "a (before)\n";
  a.print();
  a.remove(1);
  UT_ASSERT_EQUALS(a.n(), 2);
  UT_ASSERT_EQUALS(a.data().size(), 5);
  LOG << "a (after 1)\n";
  a.print();

  a.remove(1);  // remove the last element
  UT_ASSERT_EQUALS(a.n(), 1);
  UT_ASSERT_EQUALS(a.data().size(), 3);
  LOG << "a (after 2)\n";
  a.print();

  UT_CATCH_EXCEPTION(a.remove(1));
}
UT_TEST_CASE_END(test_jagged)

UT_TEST_SUITE_END(array2d_test_suite)
