#pragma once
#include <array>
#include <vector>
#include "Types.hpp"

// Public: FastEnvelope.h reaches this header through AABB.h, so it must not include
// <implicit_point.h> (or ip_filtered.h). orient_2d / orient_3d, which need it, live in the
// private orient_predicates.h.

namespace fastEnvelope {
namespace algorithms {
int seg_cut_plane(
    const Vector3& seg0,
    const Vector3& seg1,
    const Vector3& t0,
    const Vector3& t1,
    const Vector3& t2);

int is_triangle_degenerated(
    const Vector3& triangle0,
    const Vector3& triangle1,
    const Vector3& triangle2);

Vector3 accurate_normal_vector(const Vector3& p0, const Vector3& p1, const Vector3& p2);

void get_tri_corners(
    const Vector3& triangle0,
    const Vector3& triangle1,
    const Vector3& triangle2,
    Vector3& mint,
    Vector3& maxt);

// Inline: the AABB traversal and the per-prism tests call it in their innermost loops, from
// other translation units.
inline bool box_box_intersection(
    const Vector3& min1,
    const Vector3& max1,
    const Vector3& min2,
    const Vector3& max2)
{
    if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2]) return 0;
    if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2]) return 0;
    return 1;
}

Vector2 to_2d(const Vector3& p, int t);
// resort the facets using Morton's code

void resorting(
    const std::vector<Vector3>& V,
    const std::vector<Vector3i>& F,
    std::vector<Vector3i>& fnew);

void resorting(
    const std::vector<Vector3>& V,
    const std::vector<Vector3i>& F,
    std::vector<Vector3i>& fnew,
    std::vector<int>& new2old);

void resorting(
    const std::vector<Vector3>& V,
    const std::vector<Vector2i>& E,
    std::vector<Vector2i>& enew);

void resorting(
    const std::vector<Vector3>& V,
    const std::vector<Vector2i>& E,
    std::vector<Vector2i>& enew,
    std::vector<int>& new2old);

void seg_cube(
    const Vector3& p1,
    const Vector3& p2,
    const Scalar& width,
    std::array<Vector3, 8>& envbox);

void halfspace_generation(
    const std::vector<Vector3>& m_ver,
    const std::vector<Vector3i>& m_faces,
    std::vector<std::vector<std::array<Vector3, 3>>>& halfspace,
    std::vector<std::array<Vector3, 2>>& cornerlist,
    const Scalar& epsilon);
void halfspace_generation(
    const std::vector<Vector3>& m_ver,
    const std::vector<Vector3i>& m_faces,
    std::vector<std::vector<std::array<Vector3, 3>>>& halfspace,
    std::vector<std::array<Vector3, 2>>& cornerlist,
    const std::vector<Scalar>& epsilons);

void halfspace_generation(
    const std::vector<Vector3>& m_ver,
    const std::vector<Vector2i>& m_edges,
    std::vector<std::vector<std::array<Vector3, 3>>>& halfspace,
    std::vector<std::array<Vector3, 2>>& cornerlist,
    const Scalar& epsilon);
void halfspace_generation(
    const std::vector<Vector3>& m_ver,
    const std::vector<Vector2i>& m_edges,
    std::vector<std::vector<std::array<Vector3, 3>>>& halfspace,
    std::vector<std::array<Vector3, 2>>& cornerlist,
    const std::vector<Scalar>& epsilons);

void get_bb_corners(const std::vector<Vector3>& vertices, Vector3& min, Vector3& max);

} // namespace algorithms
} // namespace fastEnvelope