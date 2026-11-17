#pragma once
// orient_2d / orient_3d, inline. PRIVATE: include this only from fast-envelope's own .cpp files
// (and its tests), never from a header that FastEnvelope.h or FastEnvelope2D.h reaches.
//
// They call ::orient2d / ::orient3d, so this header needs <implicit_point.h>, and through it
// NFG's numerics.h, which turns on `#pragma STDC FENV_ACCESS ON` (MSVC: `fenv_access (on)`;
// under AVX2 also `GCC target("fma")`) for the rest of the including translation unit and
// declares orient2d / orient3d at global scope. Reached from a public header, that pragma
// would make clang compile every floating-point operation after the include in every
// consumer translation unit as a constrained (strictfp) intrinsic, which blocks
// vectorization and constant folding in code that has nothing to do with the envelope.
#include <implicit_point.h>

#include "Types.hpp"

namespace fastEnvelope {
namespace algorithms {

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
// without building the points or checking their types. Inline: the AABB traversal and the
// per-prism tests call them in their innermost loops, from several translation units.
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
