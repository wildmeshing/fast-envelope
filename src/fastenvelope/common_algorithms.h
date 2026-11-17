#pragma once
#include <array>
#include <vector>
#include "Types.hpp"

// For ::orient2d / ::orient3d, which the inline orient_2d / orient_3d below call. Like
// ip_filtered.h, this header is included only by FastEnvelope's own .cpp files.
#include <implicit_point.h>

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

// Inline, with orient_2d and orient_3d below: the AABB traversal and the per-prism tests call
// them in their innermost loops, from other translation units.
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

// These two do NOT share a handedness, and the asymmetry is not a mistake -- it is what makes
// each of them reproduce the Shewchuk-convention predicate it replaced.
//
// Indirect_Predicates orients from the FIRST argument: orient2D(a,b,c) is sign((b-a) x (c-a))
// and orient3D(a,b,c,d) is sign det[b-a, c-a, d-a]. Shewchuk orients from the LAST. In 3D that
// relabelling is an odd permutation, so the negation below lands exactly on Shewchuk's sign;
// in 2D it is an even permutation, so orient2D already agrees with Shewchuk and the negation
// makes this function its opposite.
//
// Net effect: orient_3d(p0,p1,p2,p3) is positive when p3 is below the plane <p0,p1,p2>, and
// orient_2d(p0,p1,p2) is NEGATIVE when p2 is to the left of p0->p1. Every current caller is
// invariant to a uniform sign flip -- is_triangle_cut_bounding_box compares ori * o_corner,
// is_segment_cut_bounding_box asks only whether all corners agree, is_triangle_degenerated
// tests against zero -- and FastEnvelope2D::contains is written against the convention stated
// here. New code comparing either result to an absolute sign should re-read this comment.
//
// Both call the double-precision predicate that genericPoint::orient2D / orient3D dispatch to
// for explicit points (orient2d_EEE -> ::orient2d, orient3d_EEEE -> ::orient3d), the same answer
// without building the points or checking their types.
inline int orient_2d(const Vector2& p0, const Vector2& p1, const Vector2& p2)
{
    return -::orient2d(p0[0], p0[1], p1[0], p1[1], p2[0], p2[1]);
}
inline int orient_3d(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3)
{
    return -::orient3d(
        p0[0],
        p0[1],
        p0[2],
        p1[0],
        p1[1],
        p1[2],
        p2[0],
        p2[1],
        p2[2],
        p3[0],
        p3[1],
        p3[2]);
}

} // namespace algorithms
} // namespace fastEnvelope