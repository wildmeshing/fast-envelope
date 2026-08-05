#pragma once

// Indirect predicates over line-plane (LPI) and three-plane (TPI) intersection points.
//
// This file used to be a vendored copy of an older Indirect_Predicates, carrying ~1500
// lines of hand-written expansion arithmetic. The exact arithmetic is now MarcoAttene's
// Indirect_Predicates, which VolumeRemesher also builds against -- two copies of it
// cannot share a binary, because both put `expansionObject` in the global namespace and
// disagree on whether its members are static, which the Itanium ABI mangles identically.
//
// The signatures below are unchanged, so FastEnvelope's ~93 call sites are untouched.
// What changed underneath:
//
//   - Only the *exact* arithmetic was externalized. The semi-static filter that fronts it
//     is still the vendored one, and has to be: Indirect_Predicates has no
//     double-precision tier -- `orient3d_indirect_IEEE` reaches for interval arithmetic
//     immediately, and on arm64 every interval predicate brackets itself with two FPU
//     rounding-mode switches. Routing every query through it profiled at 48% of samples
//     in `fesetround` alone, for 2.5-3.4x the runtime on wildmeshing-toolkit's
//     envelope-heavy configurations. See the tier notes below.
//
//   - The support objects used to cache the intersection point's lambda as raw doubles.
//     They still do, for the filter. What they gained is the Indirect_Predicates implicit
//     point used to escalate, built only when a query actually needs it -- plus the
//     explicit points it refers to, because implicitPoint3D_* stores references and its
//     operands must outlive it. That makes them non-copyable; FastEnvelope only ever
//     passes them by reference and never reads their members, so this is invisible to it.
//
//   - POSITIVE / NEGATIVE / UNCERTAIN are no longer their own enum. Indirect_Predicates
//     ships the same values (1, -1, 0) as `IP_Sign` / `Filtered_Sign`, and declaring them
//     twice collides, so `Filtered_Orientation` is now an alias for those.

#include <implicit_point.h>

#include <math.h>

#include <optional>

// Indirect_Predicates and this header disagree on the sign of orient3D, and FastEnvelope
// depends on the one this header has always returned.
//
// `genericPoint::orient3D(p, a, b, c)` gives the signed volume of the tetrahedron
// (p, a, b, c) in Indirect_Predicates' own vertex order, which works out to
// -sign(det(p - a, b - a, c - a)): a point one unit above the plane z = 0, tested against
// <(0,0,0), (1,0,0), (0,1,0)> -- whose (b-a) x (c-a) normal is +z -- comes back NEGATIVE.
//
// The predicates this file replaced returned the other sign: the height of the implicit
// point over the oriented plane <a, b, c>, positive above and negative below. FastEnvelope
// reads them that way. A point is inside a prism when it lies below all of that prism's
// outward-oriented faces, which its callers spell as every face returning NEGATIVE -- see
// the `tot == halfspace[prismindex[i]].size()` tests in FastEnvelope.cpp. Handing them the
// opposite sign makes containment fail everywhere: every prism is skipped, every point
// looks like it is outside, and the envelope silently grows far more conservative than
// asked for.
//
// UNCERTAIN is 0, which negates to itself, so the filtered/exact split is unaffected.
inline int ip_orient3D_to_height_sign(int orient3D_sign)
{
    return -orient3D_sign;
}

// ---------------------------------------------------------------------------
// Support objects
// ---------------------------------------------------------------------------
//
// Three tiers, each an order of magnitude dearer than the last:
//
//   1. A semi-static filter over plain doubles, with the error bounds Attene's filter
//      generator derived. Cheap, and it settles the overwhelming majority of queries.
//      Restored here verbatim from the predicates this file used to vendor -- it is
//      self-contained double arithmetic and pulls in none of the expansion code whose
//      duplicate definitions were the reason for externalizing the kernel.
//
//   2. Indirect_Predicates' interval stage, for what tier 1 gives up on. Tier 1 is blunt:
//      on a tetwild_sphere run it cannot decide roughly 1700 queries, of which intervals
//      settle all but 470.
//
//   3. Indirect_Predicates' exact cascade, reached through `_pre_exact` / `_post_exact`
//      when tier 2 also reports UNCERTAIN, exactly as the vendored predicates escalated
//      to their own expansion arithmetic.
//
// Neither of the cheap tiers is optional. Without tier 1, every query pays interval
// arithmetic, which on arm64 brackets itself with two FPU rounding-mode switches --
// profiled at 48% of samples in `fesetround` alone, for 2.5-3.4x the runtime of the
// vendored predicates. Without tier 2, everything tier 1 cannot decide falls through to
// bigfloat, which on near-degenerate inputs (where tier 1 gives up far more often) costs
// more than tier 1 saves.
//
// All three return the same thing: the height of the implicit point over the oriented
// plane, positive above and negative below. Tier 1 already speaks that convention because
// it is the original code; tiers 2 and 3 are converted with ip_orient3D_to_height_sign().

// The values the vendored predicates returned. Indirect_Predicates ships the same three
// under IP_Sign / Filtered_Sign, so this is an alias for the old spelling rather than a
// second set of constants.
namespace Filtered_Orientation {
constexpr int POSITIVE = IP_Sign::POSITIVE;
constexpr int NEGATIVE = IP_Sign::NEGATIVE;
constexpr int UNCERTAIN = Filtered_Sign::UNCERTAIN;
} // namespace Filtered_Orientation

class LPI_filtered_suppvars
{
public:
    // Tier 1 state, filled by orient3D_LPI_prefilter and read by orient3D_LPI_postfilter.
    double a11, a12, a13, d, fa11, fa12, fa13, max1, max2, max5;

    // Tier 2 state. The implicit point is built on demand: constructing one computes an
    // interval lambda, which is wasted work on every query tier 1 settles. It holds
    // references to its operands, so those are members too, which is why this object is
    // non-copyable -- FastEnvelope only ever passes it by reference.
    double line[6], plane[9];
    mutable explicitPoint3D ep, eq, er, es, et;
    mutable std::optional<implicitPoint3D_LPI> pt;

    LPI_filtered_suppvars() = default;
    LPI_filtered_suppvars(const LPI_filtered_suppvars&) = delete;
    LPI_filtered_suppvars& operator=(const LPI_filtered_suppvars&) = delete;

    void set_operands(
        double px, double py, double pz,
        double qx, double qy, double qz,
        double rx, double ry, double rz,
        double sx, double sy, double sz,
        double tx, double ty, double tz)
    {
        line[0] = px; line[1] = py; line[2] = pz;
        line[3] = qx; line[4] = qy; line[5] = qz;
        plane[0] = rx; plane[1] = ry; plane[2] = rz;
        plane[3] = sx; plane[4] = sy; plane[5] = sz;
        plane[6] = tx; plane[7] = ty; plane[8] = tz;
        pt.reset();
    }

    const implicitPoint3D_LPI& implicit() const
    {
        if (!pt) {
            ep = explicitPoint3D(line[0], line[1], line[2]);
            eq = explicitPoint3D(line[3], line[4], line[5]);
            er = explicitPoint3D(plane[0], plane[1], plane[2]);
            es = explicitPoint3D(plane[3], plane[4], plane[5]);
            et = explicitPoint3D(plane[6], plane[7], plane[8]);
            pt.emplace(ep, eq, er, es, et);
        }
        return *pt;
    }

    // Whether the point is defined at all -- the line meets the plane / the three planes
    // meet in a point -- which is exactly a non-zero denominator, and is what the vendored
    // exact pre-steps reported when theirs came out zero.
    //
    // getExactXYZCoordinates() answers this too, but divides out the homogeneous
    // coordinates to do it: three bigrational divisions per call, which profiled at 67% of
    // a simwild_double_sphere_notop_3d run (bignatural::push_bit_back and divide_by). Only
    // the denominator's sign is wanted, and the cached interval settles it without any
    // bignum work in the common case.
    bool exists() const
    {
        interval_number lx, ly, lz, den;
        if (implicit().getIntervalLambda(lx, ly, lz, den)) return true;
        bigfloat bx, by, bz, bden;
        implicit().getBigfloatLambda(bx, by, bz, bden);
        return sgn(bden) != 0;
    }

};

// A distinct class rather than an alias: FastEnvelope.h forward-declares
// `class LPI_exact_suppvars;`, which a typedef cannot satisfy. The exact and filtered
// variants carry the same state, so the derived class adds nothing but the name.
class LPI_exact_suppvars : public LPI_filtered_suppvars
{
};

class TPI_filtered_suppvars
{
public:
    // Tier 1 state.
    double d, n1, n2, n3, max1, max2, max3, max4, max5, max6, max7;

    // Tier 2 state, built on demand -- see LPI_filtered_suppvars.
    double pv[9], pw[9], pu[9];
    mutable explicitPoint3D ev1, ev2, ev3, ew1, ew2, ew3, eu1, eu2, eu3;
    mutable std::optional<implicitPoint3D_TPI> pt;

    TPI_filtered_suppvars() = default;
    TPI_filtered_suppvars(const TPI_filtered_suppvars&) = delete;
    TPI_filtered_suppvars& operator=(const TPI_filtered_suppvars&) = delete;

    void set_operands(
        double v1x, double v1y, double v1z, double v2x, double v2y, double v2z,
        double v3x, double v3y, double v3z,
        double w1x, double w1y, double w1z, double w2x, double w2y, double w2z,
        double w3x, double w3y, double w3z,
        double u1x, double u1y, double u1z, double u2x, double u2y, double u2z,
        double u3x, double u3y, double u3z)
    {
        pv[0] = v1x; pv[1] = v1y; pv[2] = v1z;
        pv[3] = v2x; pv[4] = v2y; pv[5] = v2z;
        pv[6] = v3x; pv[7] = v3y; pv[8] = v3z;
        pw[0] = w1x; pw[1] = w1y; pw[2] = w1z;
        pw[3] = w2x; pw[4] = w2y; pw[5] = w2z;
        pw[6] = w3x; pw[7] = w3y; pw[8] = w3z;
        pu[0] = u1x; pu[1] = u1y; pu[2] = u1z;
        pu[3] = u2x; pu[4] = u2y; pu[5] = u2z;
        pu[6] = u3x; pu[7] = u3y; pu[8] = u3z;
        pt.reset();
    }

    const implicitPoint3D_TPI& implicit() const
    {
        if (!pt) {
            ev1 = explicitPoint3D(pv[0], pv[1], pv[2]);
            ev2 = explicitPoint3D(pv[3], pv[4], pv[5]);
            ev3 = explicitPoint3D(pv[6], pv[7], pv[8]);
            ew1 = explicitPoint3D(pw[0], pw[1], pw[2]);
            ew2 = explicitPoint3D(pw[3], pw[4], pw[5]);
            ew3 = explicitPoint3D(pw[6], pw[7], pw[8]);
            eu1 = explicitPoint3D(pu[0], pu[1], pu[2]);
            eu2 = explicitPoint3D(pu[3], pu[4], pu[5]);
            eu3 = explicitPoint3D(pu[6], pu[7], pu[8]);
            pt.emplace(ev1, ev2, ev3, ew1, ew2, ew3, eu1, eu2, eu3);
        }
        return *pt;
    }

    // Whether the point is defined at all -- the line meets the plane / the three planes
    // meet in a point -- which is exactly a non-zero denominator, and is what the vendored
    // exact pre-steps reported when theirs came out zero.
    //
    // getExactXYZCoordinates() answers this too, but divides out the homogeneous
    // coordinates to do it: three bigrational divisions per call, which profiled at 67% of
    // a simwild_double_sphere_notop_3d run (bignatural::push_bit_back and divide_by). Only
    // the denominator's sign is wanted, and the cached interval settles it without any
    // bignum work in the common case.
    bool exists() const
    {
        interval_number lx, ly, lz, den;
        if (implicit().getIntervalLambda(lx, ly, lz, den)) return true;
        bigfloat bx, by, bz, bden;
        implicit().getBigfloatLambda(bx, by, bz, bden);
        return sgn(bden) != 0;
    }

};

class TPI_exact_suppvars : public TPI_filtered_suppvars
{
};

// ---------------------------------------------------------------------------
// Tier 2: Indirect_Predicates' interval stage
// ---------------------------------------------------------------------------
//
// Sits between the double filter and the exact cascade. The double filter is cheap but
// blunt -- on a tetwild_sphere run it gives up on roughly 1700 queries where interval
// arithmetic settles all but 470 of them. Without this tier each of those becomes a
// bigfloat evaluation in the caller's recompute pass, which is what the vendored
// predicates used their own expansion arithmetic for; on near-degenerate inputs, where
// the double filter gives up far more often, that dominates the run.
//
// Reached rarely, so the FPU rounding-mode switches it carries cost little here. Returns
// UNCERTAIN when the intervals straddle zero, leaving the caller to escalate as before.
inline int ip_interval_tier(
    const genericPoint& ip,
    const double& ax, const double& ay, const double& az,
    const double& bx, const double& by, const double& bz,
    const double& cx, const double& cy, const double& cz)
{
    return ip_orient3D_to_height_sign(
        orient3d_indirect_IEEE_t<interval_number, interval_number>(
            ip, ax, ay, az, bx, by, bz, cx, cy, cz));
}

// ---------------------------------------------------------------------------
// Tier 1: semi-static filter over doubles
// ---------------------------------------------------------------------------
//
// Bodies below are the vendored predicates unchanged, minus the USE_MULTISTAGE_FILTERS
// branches, which this project never compiled in. Do not "simplify" the arithmetic: the
// error bounds are matched to this exact operation order.

inline bool orient3D_LPI_prefilter(
    const double& px, const double& py, const double& pz,
    const double& qx, const double& qy, const double& qz,
    const double& rx, const double& ry, const double& rz,
    const double& sx, const double& sy, const double& sz,
    const double& tx, const double& ty, const double& tz,
    LPI_filtered_suppvars& svs)
{
    svs.set_operands(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz);

    double a21, a22, a23, a31, a32, a33;

    double &a11 = svs.a11;
    double &a12 = svs.a12;
    double &a13 = svs.a13;
    double &d = svs.d;
    double &fa11 = svs.fa11;
    double &fa12 = svs.fa12;
    double &fa13 = svs.fa13;
    double &max1 = svs.max1;
    double &max2 = svs.max2;
    double &max5 = svs.max5;

    a11 = (px - qx);
    a12 = (py - qy);
    a13 = (pz - qz);
    a21 = (sx - rx);
    a22 = (sy - ry);
    a23 = (sz - rz);
    a31 = (tx - rx);
    a32 = (ty - ry);
    a33 = (tz - rz);
    double a2233 = ((a22 * a33) - (a23 * a32));
    double a2133 = ((a21 * a33) - (a23 * a31));
    double a2132 = ((a21 * a32) - (a22 * a31));
    d = (((a11 * a2233) - (a12 * a2133)) + (a13 * a2132));

    fa11 = fabs(a11);
    double fa21 = fabs(a21);
    double fa31 = fabs(a31);
    fa12 = fabs(a12);
    double fa22 = fabs(a22);
    double fa32 = fabs(a32);
    fa13 = fabs(a13);
    double fa23 = fabs(a23);
    double fa33 = fabs(a33);

    max1 = fa23;
    if (max1 < fa13)
        max1 = fa13;
    if (max1 < fa33)
        max1 = fa33;

    max2 = fa12;
    if (max2 < fa22)
        max2 = fa22;
    if (max2 < fa32)
        max2 = fa32;

    max5 = fa11;
    if (max5 < fa21)
        max5 = fa21;
    if (max5 < fa31)
        max5 = fa31;

    double deps = 5.1107127829973299e-015 * max1 * max2 * max5;
    if (!isfinite(d) || (d <= deps && d >= -deps))
        return false;

    double px_rx = px - rx;
    double py_ry = py - ry;
    double pz_rz = pz - rz;

    double n = ((((py_ry)*a2133) - ((px_rx)*a2233)) - ((pz_rz)*a2132));

    a11 *= n;
    a12 *= n;
    a13 *= n;

    if (!isfinite(a11) || !isfinite(a12) || !isfinite(a13))
        return false;

    double fpxrx = fabs(px_rx);
    double fpyry = fabs(py_ry);
    double fpzrz = fabs(pz_rz);

    if (max1 < fpzrz)
        max1 = fpzrz;
    if (max2 < fpyry)
        max2 = fpyry;
    if (max5 < fpxrx)
        max5 = fpxrx;


    return true;
}

// `px, py, pz` is the line's first point, which the prefilter has already folded into the
// support variables; it is passed again because the determinant below is expressed
// relative to it.
inline int orient3D_LPI_postfilter(
    const LPI_filtered_suppvars& svs,
    const double& px, const double& py, const double& pz,
    const double& ax, const double& ay, const double& az,
    const double& bx, const double& by, const double& bz,
    const double& cx, const double& cy, const double& cz)
{
    double px_cx, py_cy, pz_cz;
    double d11, d12, d13, d21, d31, d22, d32, d23, d33;
    double d2233, d2332, d2133, d2331, d2132, d2231;
    double det;

    const double &a11 = svs.a11, &a12 = svs.a12, &a13 = svs.a13;
    const double &d = svs.d;
    const double &fa11 = svs.fa11, &fa12 = svs.fa12, &fa13 = svs.fa13;
    const double &max1 = svs.max1, &max2 = svs.max2, &max5 = svs.max5;

    px_cx = px - cx;
    py_cy = py - cy;
    pz_cz = pz - cz;

    d11 = (d * px_cx) + (a11);
    d21 = (ax - cx);
    d31 = (bx - cx);
    d12 = (d * py_cy) + (a12);
    d22 = (ay - cy);
    d32 = (by - cy);
    d13 = (d * pz_cz) + (a13);
    d23 = (az - cz);
    d33 = (bz - cz);

    d2233 = d22 * d33;
    d2332 = d23 * d32;
    d2133 = d21 * d33;
    d2331 = d23 * d31;
    d2132 = d21 * d32;
    d2231 = d22 * d31;

    det = d11 * (d2233 - d2332) - d12 * (d2133 - d2331) + d13 * (d2132 - d2231);

    if (!isfinite(det))
        return ip_interval_tier(svs.implicit(), ax, ay, az, bx, by, bz, cx, cy, cz);

    double fd11 = fabs(d11);
    double fd21 = fabs(d21);
    double fd31 = fabs(d31);
    double fd12 = fabs(d12);
    double fd22 = fabs(d22);
    double fd32 = fabs(d32);
    double fd13 = fabs(d13);
    double fd23 = fabs(d23);
    double fd33 = fabs(d33);
    double fpxcx = fabs(px_cx);
    double fpycy = fabs(py_cy);
    double fpzcz = fabs(pz_cz);

    double max3, max4, max6;

    max3 = fa12;
    if (max3 < fd32)
        max3 = fd32;
    if (max3 < fd22)
        max3 = fd22;
    if (max3 < fpycy)
        max3 = fpycy;

    max4 = fa13;
    if (max4 < fpzcz)
        max4 = fpzcz;
    if (max4 < fd33)
        max4 = fd33;
    if (max4 < fd23)
        max4 = fd23;

    max6 = fa11;
    if (max6 < fd21)
        max6 = fd21;
    if (max6 < fd31)
        max6 = fd31;
    if (max6 < fpxcx)
        max6 = fpxcx;

    double eps = 1.3865993466947057e-013 * max1 * max2 * max5 * max6 * max3 * max4;

    if ((det > eps))
        return (d > 0) ? (Filtered_Orientation::POSITIVE) : (Filtered_Orientation::NEGATIVE);
    if ((det < -eps))
        return (d > 0) ? (Filtered_Orientation::NEGATIVE) : (Filtered_Orientation::POSITIVE);

    return ip_interval_tier(svs.implicit(), ax, ay, az, bx, by, bz, cx, cy, cz);
}

inline bool orient3D_TPI_prefilter(
    const double& ov1x, const double& ov1y, const double& ov1z,
    const double& ov2x, const double& ov2y, const double& ov2z,
    const double& ov3x, const double& ov3y, const double& ov3z,
    const double& ow1x, const double& ow1y, const double& ow1z,
    const double& ow2x, const double& ow2y, const double& ow2z,
    const double& ow3x, const double& ow3y, const double& ow3z,
    const double& ou1x, const double& ou1y, const double& ou1z,
    const double& ou2x, const double& ou2y, const double& ou2z,
    const double& ou3x, const double& ou3y, const double& ou3z,
    TPI_filtered_suppvars& svs)
{
    svs.set_operands(ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x, ov3y, ov3z,
                     ow1x, ow1y, ow1z, ow2x, ow2y, ow2z, ow3x, ow3y, ow3z,
                     ou1x, ou1y, ou1z, ou2x, ou2y, ou2z, ou3x, ou3y, ou3z);

    double &d = svs.d, &n1 = svs.n1, &n2 = svs.n2, &n3 = svs.n3;
    double &max1 = svs.max1, &max2 = svs.max2, &max3 = svs.max3, &max4 = svs.max4, &max5 = svs.max5, &max6 = svs.max6, &max7 = svs.max7;

    double v3x = ov3x - ov2x;
    double v3y = ov3y - ov2y;
    double v3z = ov3z - ov2z;
    double v2x = ov2x - ov1x;
    double v2y = ov2y - ov1y;
    double v2z = ov2z - ov1z;
    double w3x = ow3x - ow2x;
    double w3y = ow3y - ow2y;
    double w3z = ow3z - ow2z;
    double w2x = ow2x - ow1x;
    double w2y = ow2y - ow1y;
    double w2z = ow2z - ow1z;
    double u3x = ou3x - ou2x;
    double u3y = ou3y - ou2y;
    double u3z = ou3z - ou2z;
    double u2x = ou2x - ou1x;
    double u2y = ou2y - ou1y;
    double u2z = ou2z - ou1z;

    double nvx = v2y * v3z - v2z * v3y;
    double nvy = v3x * v2z - v3z * v2x;
    double nvz = v2x * v3y - v2y * v3x;

    double nwx = w2y * w3z - w2z * w3y;
    double nwy = w3x * w2z - w3z * w2x;
    double nwz = w2x * w3y - w2y * w3x;

    double nux = u2y * u3z - u2z * u3y;
    double nuy = u3x * u2z - u3z * u2x;
    double nuz = u2x * u3y - u2y * u3x;

    double nwyuz = nwy * nuz - nwz * nuy;
    double nwxuz = nwx * nuz - nwz * nux;
    double nwxuy = nwx * nuy - nwy * nux;

    d = nvx * nwyuz - nvy * nwxuz + nvz * nwxuy;

    // Almost static filter for d
    double fv2x = fabs(v2x);
    double fv2y = fabs(v2y);
    double fv2z = fabs(v2z);
    double fv3x = fabs(v3x);
    double fv3y = fabs(v3y);
    double fv3z = fabs(v3z);

    double fw2x = fabs(w2x);
    double fw2y = fabs(w2y);
    double fw2z = fabs(w2z);
    double fw3x = fabs(w3x);
    double fw3y = fabs(w3y);
    double fw3z = fabs(w3z);

    double fu2x = fabs(u2x);
    double fu3x = fabs(u3x);
    double fu2y = fabs(u2y);
    double fu2z = fabs(u2z);
    double fu3y = fabs(u3y);
    double fu3z = fabs(u3z);

    max4 = fv2y;
    if (max4 < fv3y)
        max4 = fv3y;
    if (max4 < fw3y)
        max4 = fw3y;
    if (max4 < fw2y)
        max4 = fw2y;
    max2 = fv3x;
    if (max2 < fv2x)
        max2 = fv2x;
    if (max2 < fw2x)
        max2 = fw2x;
    if (max2 < fw3x)
        max2 = fw3x;
    max5 = fv2z;
    if (max5 < fv3z)
        max5 = fv3z;
    if (max5 < fw3z)
        max5 = fw3z;
    if (max5 < fw2z)
        max5 = fw2z;
    max7 = fu2x;
    if (max7 < fu3x)
        max7 = fu3x;
    if (max7 < fw2x)
        max7 = fw2x;
    if (max7 < fw3x)
        max7 = fw3x;

    double max9 = fu2y;
    if (max9 < fu3y)
        max9 = fu3y;
    if (max9 < fw2y)
        max9 = fw2y;
    if (max9 < fw3y)
        max9 = fw3y;
    double max10 = fu2z;
    if (max10 < fu3z)
        max10 = fu3z;
    if (max10 < fw2z)
        max10 = fw2z;
    if (max10 < fw3z)
        max10 = fw3z;

    double deps = 8.8881169117764924e-014 * (((((max4 * max5) * max2) * max10) * max7) * max9);
    if (!isfinite(d) || (d <= deps && d >= -deps))
        return false;

    double nvyuz = nvy * nuz - nvz * nuy;
    double nvxuz = nvx * nuz - nvz * nux;
    double nvxuy = nvx * nuy - nvy * nux;

    double nvywz = nvy * nwz - nvz * nwy;
    double nvxwz = nvx * nwz - nvz * nwx;
    double nvxwy = nvx * nwy - nvy * nwx;

    double p1 = nvx * ov1x + nvy * ov1y + nvz * ov1z;
    double p2 = nwx * ow1x + nwy * ow1y + nwz * ow1z;
    double p3 = nux * ou1x + nuy * ou1y + nuz * ou1z;

    n1 = p1 * nwyuz - p2 * nvyuz + p3 * nvywz;
    n2 = p2 * nvxuz - p3 * nvxwz - p1 * nwxuz;
    n3 = p3 * nvxwy - p2 * nvxuy + p1 * nwxuy;

    if (!isfinite(n1) || !isfinite(n2) || !isfinite(n3))
        return false;

    if (max4 < fu2y)
        max4 = fu2y;
    if (max4 < fu3y)
        max4 = fu3y;
    if (max2 < fu2x)
        max2 = fu2x;
    if (max2 < fu3x)
        max2 = fu3x;
    if (max5 < fu2z)
        max5 = fu2z;
    if (max5 < fu3z)
        max5 = fu3z;

    max1 = max4;
    if (max1 < max2)
        max1 = max2;
    max3 = max5;
    if (max3 < max4)
        max3 = max4;
    max6 = fu2x;
    if (max6 < fu3x)
        max6 = fu3x;
    if (max6 < fu2z)
        max6 = fu2z;
    if (max6 < fw3y)
        max6 = fw3y;
    if (max6 < fw2x)
        max6 = fw2x;
    if (max6 < fw3z)
        max6 = fw3z;
    if (max6 < fw2y)
        max6 = fw2y;
    if (max6 < fw2z)
        max6 = fw2z;
    if (max6 < fu2y)
        max6 = fu2y;
    if (max6 < fu3z)
        max6 = fu3z;
    if (max6 < fu3y)
        max6 = fu3y;
    if (max6 < fw3x)
        max6 = fw3x;

    return true;
}

inline int orient3D_TPI_postfilter(
    const TPI_filtered_suppvars& svs,
    const double& q1x, const double& q1y, const double& q1z,
    const double& q2x, const double& q2y, const double& q2z,
    const double& q3x, const double& q3y, const double& q3z)
{

    const double &d = svs.d, &n1 = svs.n1, &n2 = svs.n2, &n3 = svs.n3;
    const double &max1 = svs.max1, &max2 = svs.max2, &max3 = svs.max3, &max4 = svs.max4, &max5 = svs.max5, &max6 = svs.max6, &max7 = svs.max7;

    double dq3x = d * q3x;
    double dq3y = d * q3y;
    double dq3z = d * q3z;

    double a11 = n1 - dq3x;
    double a12 = n2 - dq3y;
    double a13 = n3 - dq3z;
    double a21 = q1x - q3x;
    double a22 = q1y - q3y;
    double a23 = q1z - q3z;
    double a31 = q2x - q3x;
    double a32 = q2y - q3y;
    double a33 = q2z - q3z;

    double det = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31);
    bool infinite = true;
    infinite = isfinite(det);
    if (!infinite)
        return ip_interval_tier(svs.implicit(), q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);

    double fa21 = fabs(a21);
    double fa22 = fabs(a22);
    double fa23 = fabs(a23);
    double fa31 = fabs(a31);
    double fa32 = fabs(a32);
    double fa33 = fabs(a33);

    double nmax7 = max7;
    if (nmax7 < fa21)
        nmax7 = fa21;
    if (nmax7 < fa31)
        nmax7 = fa31;

    double nmax6 = max6;
    if (nmax6 < fa22)
        nmax6 = fa22;
    if (nmax6 < fa32)
        nmax6 = fa32;

    double max8 = fa22;
    if (max8 < fa23)
        max8 = fa23;
    if (max8 < fa33)
        max8 = fa33;
    if (max8 < fa32)
        max8 = fa32;

    double eps = 3.4025182954957945e-012 * (((((((max1 * max3) * max2) * max5) * nmax7) * max4) * nmax6) * max8);

    if ((det > eps))
        return (d > 0) ? (Filtered_Orientation::POSITIVE) : (Filtered_Orientation::NEGATIVE);
    if ((det < -eps))
        return (d > 0) ? (Filtered_Orientation::NEGATIVE) : (Filtered_Orientation::POSITIVE);
    return ip_interval_tier(svs.implicit(), q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
}

// ---------------------------------------------------------------------------
// Tier 2: Indirect_Predicates' cascade (interval -> bigfloat)
// ---------------------------------------------------------------------------
//
// The "pre" step reports whether the implicit point exists at all -- the line parallel to
// the plane, or three planes without a common point -- which is what the vendored exact
// pre-steps reported when their denominator came out zero.

inline bool orient3D_LPI_pre_exact(
    double px, double py, double pz,
    double qx, double qy, double qz,
    double rx, double ry, double rz,
    double sx, double sy, double sz,
    double tx, double ty, double tz,
    LPI_exact_suppvars& svs)
{
    svs.set_operands(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz);
    return svs.exists();
}

inline int orient3D_LPI_post_exact(
    LPI_exact_suppvars& svs,
    double, double, double,
    double ax, double ay, double az,
    double bx, double by, double bz,
    double cx, double cy, double cz)
{
    const explicitPoint3D a(ax, ay, az), b(bx, by, bz), c(cx, cy, cz);
    return ip_orient3D_to_height_sign(genericPoint::orient3D(svs.implicit(), a, b, c));
}

inline bool orient3D_TPI_pre_exact(
    double v1x, double v1y, double v1z, double v2x, double v2y, double v2z,
    double v3x, double v3y, double v3z,
    double w1x, double w1y, double w1z, double w2x, double w2y, double w2z,
    double w3x, double w3y, double w3z,
    double u1x, double u1y, double u1z, double u2x, double u2y, double u2z,
    double u3x, double u3y, double u3z,
    TPI_exact_suppvars& svs)
{
    svs.set_operands(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
                     w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
                     u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z);
    return svs.exists();
}

inline int orient3D_TPI_post_exact(
    TPI_exact_suppvars& svs,
    double q1x, double q1y, double q1z,
    double q2x, double q2y, double q2z,
    double q3x, double q3y, double q3z)
{
    const explicitPoint3D q1(q1x, q1y, q1z), q2(q2x, q2y, q2z), q3(q3x, q3y, q3z);
    return ip_orient3D_to_height_sign(genericPoint::orient3D(svs.implicit(), q1, q2, q3));
}

// ---------------------------------------------------------------------------
// One-shot interface: both tiers behind a single call
// ---------------------------------------------------------------------------

inline int orient3D_LPI_filtered(
    double px, double py, double pz, double qx, double qy, double qz,
    double rx, double ry, double rz, double sx, double sy, double sz,
    double tx, double ty, double tz,
    double ax, double ay, double az, double bx, double by, double bz,
    double cx, double cy, double cz)
{
    LPI_filtered_suppvars svs;
    if (!orient3D_LPI_prefilter(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, svs))
        return Filtered_Orientation::UNCERTAIN;
    return orient3D_LPI_postfilter(svs, px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz);
}

inline int orient3D_TPI_filtered(
    double v1x, double v1y, double v1z, double v2x, double v2y, double v2z,
    double v3x, double v3y, double v3z,
    double w1x, double w1y, double w1z, double w2x, double w2y, double w2z,
    double w3x, double w3y, double w3z,
    double u1x, double u1y, double u1z, double u2x, double u2y, double u2z,
    double u3x, double u3y, double u3z,
    double q1x, double q1y, double q1z, double q2x, double q2y, double q2z,
    double q3x, double q3y, double q3z)
{
    TPI_filtered_suppvars svs;
    if (!orient3D_TPI_prefilter(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
                                w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
                                u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, svs))
        return Filtered_Orientation::UNCERTAIN;
    return orient3D_TPI_postfilter(svs, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
}

inline int orient3D_LPI_exact(
    double px, double py, double pz, double qx, double qy, double qz,
    double rx, double ry, double rz, double sx, double sy, double sz,
    double tx, double ty, double tz,
    double ax, double ay, double az, double bx, double by, double bz,
    double cx, double cy, double cz)
{
    LPI_exact_suppvars svs;
    if (!orient3D_LPI_pre_exact(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz, svs))
        return IP_Sign::ZERO;
    return orient3D_LPI_post_exact(svs, px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz);
}

inline int orient3D_TPI_exact(
    double v1x, double v1y, double v1z, double v2x, double v2y, double v2z,
    double v3x, double v3y, double v3z,
    double w1x, double w1y, double w1z, double w2x, double w2y, double w2z,
    double w3x, double w3y, double w3z,
    double u1x, double u1y, double u1z, double u2x, double u2y, double u2z,
    double u3x, double u3y, double u3z,
    double q1x, double q1y, double q1z, double q2x, double q2y, double q2z,
    double q3x, double q3y, double q3z)
{
    TPI_exact_suppvars svs;
    if (!orient3D_TPI_pre_exact(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
                                w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
                                u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z, svs))
        return IP_Sign::ZERO;
    return orient3D_TPI_post_exact(svs, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
}

#include "ip_filtered_ex.h"
