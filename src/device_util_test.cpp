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
#include "device_util.h"

#include "tester.h"

using namespace vortex;

UT_TEST_SUITE(device_util_test_suite)

UT_TEST_CASE(device_vector_test) {
  // basic tests
  {
    device_vector<int> x(300);

    UT_ASSERT_EQUALS(x.size(), 0);
    UT_ASSERT(x.empty());

    x.push_back(1);
    UT_ASSERT_EQUALS(x.capacity(), 300);
    UT_ASSERT_EQUALS(x.size(), 1);

    for (int i = 0; i < 200; i++) {
      x.push_back(i + 2);
    }
    UT_ASSERT_EQUALS(x.size(), 201);
    UT_ASSERT_EQUALS(x.capacity(), 300);

    for (size_t i = 0; i < x.size(); i++) {
      UT_ASSERT_EQUALS(size_t(x[i]), i + 1);
    }

    UT_ASSERT_EQUALS(x.back(), 201);

    x.clear();
    UT_ASSERT(x.empty());
    UT_ASSERT_EQUALS(x.capacity(), 300);
  }

  // test resize
  {
    device_vector<int> x(10);
    UT_ASSERT_EQUALS(x.capacity(), 10);
    UT_ASSERT_EQUALS(x.size(), 0);
    for (int i = 0; i < 10; i++) {
      x.push_back(10 * i + 2);
    }
    UT_ASSERT_EQUALS(x.capacity(), 10);
    UT_ASSERT_EQUALS(x.size(), 10);
    x.push_back(234);
    UT_ASSERT_EQUALS(x.capacity(), 20);
    UT_ASSERT_EQUALS(x.size(), 11);
    for (size_t i = 0; i < 10; i++) {
      UT_ASSERT_EQUALS(size_t(x[i]), 10 * i + 2);
    }
    UT_ASSERT_EQUALS(x[10], 234);
  }
}
UT_TEST_CASE_END(device_vector_test)

UT_TEST_CASE(device_hash_set_test) {
  device_hash_set<int32_t> set;
  set.insert(11);

  for (int i = 0; i < 10; i++) {
    set.insert(3 * i);
  }

  UT_ASSERT(set.contains(0));
  UT_ASSERT(set.contains(3));
  UT_ASSERT(set.contains(6));
  UT_ASSERT(set.contains(9));
  UT_ASSERT(set.contains(11));
  UT_ASSERT(!set.contains(13));
  set.display(true);
}
UT_TEST_CASE_END(device_hash_set_test)

UT_TEST_SUITE_END(device_util_test_suite)