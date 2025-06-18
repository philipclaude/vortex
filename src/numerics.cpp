//
//  vortex: Voronoi mesher and fluid simulator for the Earth's oceans and
//  atmosphere.
//
//  Copyright 2023 - 2025 Philip Claude Caplan
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
#include "numerics.h"

#include "predicates.h"

#if ROBUST_PREDICATES == 0

void exactinit() {}
double orient2d(double* a, double* b, double* c) {
  const double acx = a[0] - c[0];
  const double bcx = b[0] - c[0];
  const double acy = a[1] - c[1];
  const double bcy = b[1] - c[1];
  return acx * bcy - acy * bcx;
}

#endif