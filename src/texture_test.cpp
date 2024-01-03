//
//  vortex: Voronoi mesher and fluid simulator for the Earth's oceans and
//  atmosphere.
//
//  Copyright 2023 - 2024 Philip Claude Caplan
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
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
