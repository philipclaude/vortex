#pragma once
#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Initializes the exact geometric predicates. Implemented in
 * external/triangle.c.
 */
void exactinit();

/**
 * \brief Initializes the exact geometric predicates. Implemented in
 * external/predicates.cxx.
 */
// void exactinit(int verbose, int noexact, int nofilter, double maxx, double
// maxy,
//                double maxz);

/**
 * \private 2d orientation test, implemented in external/predicates.cxx
 */
double orient2d(double* pa, double* pb, double* pc);

#ifdef __cplusplus
}
#endif

/**
 * \brief Determines whether pd is above or below the the plane passing through
 * the triangle pa-pb-pc. This is in fact negative six times the volume of the
 * tetrahedron pa-pb-pc-pd (oriented in that way). Returns > 0 if the tet has a
 * positive volume, < 0 if the tet has a negative volume. Returns 0 if the tet
 * has parallel edges. This calculation is **exact** if exactinit() is called
 * first. To calculate the volume of a tetrahedron, you can use
 * -orient3d(pa,pb,pc,pd) / 6.0
 *
 * \param[in] pa - pointer to first vertex coordinates
 * \param[in] pb - pointer to second vertex coordinates
 * \param[in] pc - pointer to third vertex coordinates
 * \param[in] pd - pointer to query vertex coordinates
 *
 * \return -6 * signed volume of tetrahedron (pa,pb,pc,pd)
 */
double orient3d(const double* pa, const double* pb, const double* pc,
                const double* pd);

/**
 * \private 2d incircle test, implemented in external/predicates.cxx
 */
double incircle(double* pa, double* pb, double* pc, double* pd);

/**
 * \private 3d insphere test, implemented in external/predicates.cxx
 */
double insphere(double* pa, double* pb, double* pc, double* pd, double* pe);

namespace vortex {

/**
 * \brief Determines whether pc is above or below the the line defined by pa to
 * pb. This is in fact twice the area of the triangle pa-pb-pc (oriented in that
 * way). Returns > 0 if the triangle has a positive area, < 0 if the triangle
 * has a negative area. Returns 0 if the triangle has parallel edges. This
 * calculation is **exact** if exactinit() is called first. To calculate the
 * area of a triangle, you can use orient2d(pa,pb,pc) / 2.0
 *
 * \param[in] pa - pointer to first vertex coordinates
 * \param[in] pb - pointer to second vertex coordinates
 * \param[in] pc - pointer to query coordinates
 *
 * \return 2 * signed area of triangle (pa,pb,pc)
 */
inline double orient2d(const double* pa, const double* pb, const double* pc) {
  return ::orient2d(const_cast<double*>(pa), const_cast<double*>(pb),
                    const_cast<double*>(pc));
}

/**
 * \brief Determines whether pd is inside or outside the circle passing through
 * pa-pb-pc. Returns > 0 if the pd is inside, < 0 if pd is outside and 0 if it
 * is exactly on the circle. This calculation is **exact** if exactinit() is
 * called first.
 *
 * \param[in] pa - pointer to first vertex coordinates
 * \param[in] pb - pointer to second vertex coordinates
 * \param[in] pc - pointer to third vertex coordinates
 * \param[in] pd - pointer to query vertex coordinates
 *
 * \return whether pd is inside or outside through the sign of the result
 */
inline double incircle(const double* pa, const double* pb, const double* pc,
                       const double* pd) {
  return ::incircle(const_cast<double*>(pa), const_cast<double*>(pb),
                    const_cast<double*>(pc), const_cast<double*>(pd));
}

/**
 * \brief Determines whether pe is inside or outside the sphere passing through
 * pa-pb-pc-pd. Returns > 0 if the pe is inside, and < 0 if pe is outside. This
 * calculation is **exact** if exactinit() is called first.
 *
 * \param[in] pa - pointer to first vertex coordinates
 * \param[in] pb - pointer to second vertex coordinates
 * \param[in] pc - pointer to third vertex coordinates
 * \param[in] pd - pointer to fourth vertex coordinates
 * \param[in] pe - pointer to query vertex coordinates
 *
 * \return whether pe is inside or outside through the sign of the result
 */
inline double insphere(const double* pa, const double* pb, const double* pc,
                       const double* pd, const double* pe) {
  return ::insphere(const_cast<double*>(pa), const_cast<double*>(pb),
                    const_cast<double*>(pc), const_cast<double*>(pd),
                    const_cast<double*>(pe));
}

/**
 * \brief Initializes the predicates.
 */
inline void initialize_predicates() {
  exactinit();
  // exactinit(0, 0, 0, 10, 10, 10);
}

}  // namespace vortex
