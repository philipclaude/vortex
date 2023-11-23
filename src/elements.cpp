#include "elements.h"

namespace vortex {

int Quad::faces[8] = {0, 1, 1, 2, 2, 3, 3, 0};
int Triangle::faces[6] = {1, 2, 2, 0, 0, 1};
int Triangle::edges[6] = {0, 1, 1, 2, 2, 0};

}