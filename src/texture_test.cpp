#include "texture.h"

#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(texture_test_suite)

UT_TEST_CASE(read_test) {
  TextureOptions options;
  options.format = TextureFormat::kRGBA;
  options.flipy = true;
  Texture texture("earth.jpg", options);

  UT_ASSERT_EQUALS(texture.width(), 2048);
  UT_ASSERT_EQUALS(texture.height(), 1024);
  UT_ASSERT_EQUALS(texture.channels(), 4);
  for (int j = 0; j < texture.height(); j++) {
    for (int i = 0; i < texture.width(); i++) {
      double rgba[4];
      double s = double(i) / texture.width();
      double t = double(j) / texture.height();
      texture.sample(s, t, rgba);
    }
  }
}
UT_TEST_CASE_END(read_test)

UT_TEST_SUITE_END(texture_test_suite)
