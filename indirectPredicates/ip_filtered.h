#pragma once

// Indirect predicates over line-plane (LPI) and three-plane (TPI) intersection points.
//
// This file used to be a vendored copy of an older Indirect_Predicates, with ~1500 lines
// of hand-written expansion arithmetic behind it. It is now a thin adapter over
// MarcoAttene's Indirect_Predicates, which computes the same predicates with its own
// cascade (semi-static filter -> interval -> expansion -> bigfloat).
//
// The signatures below are unchanged, so FastEnvelope's ~93 call sites are untouched.
// What changed underneath:
//
//   - The support objects used to cache the intersection point's lambda as raw doubles.
//     Indirect_Predicates caches it inside the implicit point itself, so a support object
//     now just owns that point -- plus the explicit points it refers to, because
//     implicitPoint3D_* stores references and its operands must outlive it. That makes
//     them non-copyable; FastEnvelope only ever passes them by reference and never reads
//     their members, so the change is invisible to it.
//
//   - The filtered / exact split is preserved, and it matters. `_postfilter` runs only the
//     interval stage and returns UNCERTAIN (0) when it cannot decide, exactly as before;
//     `_post_exact` runs the full cascade. Routing both through genericPoint::orient3D
//     instead is correct but escalates to bigfloat on every call, and FastEnvelope calls
//     the postfilter in its innermost loops -- measured at more than 8x the whole
//     integration suite's runtime before it was interrupted.
//
//   - POSITIVE / NEGATIVE / UNCERTAIN are no longer declared here. Indirect_Predicates
//     ships them as a global `IP_Sign` / `Filtered_Sign` with the same values (1, -1, 0),
//     and declaring them twice collides.

#include <implicit_point.h>

#include <optional>

// ---------------------------------------------------------------------------
// Support objects
// ---------------------------------------------------------------------------

class LPI_filtered_suppvars
{
public:
    // Operands first: the implicit point below holds references to them.
    explicitPoint3D p, q, r, s, t;
    std::optional<implicitPoint3D_LPI> point;

    LPI_filtered_suppvars() = default;
    LPI_filtered_suppvars(const LPI_filtered_suppvars&) = delete;
    LPI_filtered_suppvars& operator=(const LPI_filtered_suppvars&) = delete;

    // Returns whether the interval filter could determine the lambda, which is what the
    // old prefilter reported.
    bool init(
        double px, double py, double pz,
        double qx, double qy, double qz,
        double rx, double ry, double rz,
        double sx, double sy, double sz,
        double tx, double ty, double tz)
    {
        p = explicitPoint3D(px, py, pz);
        q = explicitPoint3D(qx, qy, qz);
        r = explicitPoint3D(rx, ry, rz);
        s = explicitPoint3D(sx, sy, sz);
        t = explicitPoint3D(tx, ty, tz);
        point.emplace(p, q, r, s, t);

        interval_number lx, ly, lz, d;
        return point->getIntervalLambda(lx, ly, lz, d);
    }

    // Exact counterpart: false when the point is undefined (line parallel to the plane),
    // which is what the old pre_exact reported.
    bool init_exact(
        double px, double py, double pz,
        double qx, double qy, double qz,
        double rx, double ry, double rz,
        double sx, double sy, double sz,
        double tx, double ty, double tz)
    {
        init(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz);
        bigrational ex, ey, ez;
        return point->getExactXYZCoordinates(ex, ey, ez);
    }
};

// A distinct class rather than an alias: FastEnvelope.h forward-declares
// `class TPI_exact_suppvars;`, which a typedef cannot satisfy. The exact and filtered
// variants carry the same state now, so the derived class adds nothing but the name.
class LPI_exact_suppvars : public LPI_filtered_suppvars
{
};

class TPI_filtered_suppvars
{
public:
    explicitPoint3D v1, v2, v3, w1, w2, w3, u1, u2, u3;
    std::optional<implicitPoint3D_TPI> point;

    TPI_filtered_suppvars() = default;
    TPI_filtered_suppvars(const TPI_filtered_suppvars&) = delete;
    TPI_filtered_suppvars& operator=(const TPI_filtered_suppvars&) = delete;

    bool init(
        double v1x, double v1y, double v1z, double v2x, double v2y, double v2z,
        double v3x, double v3y, double v3z,
        double w1x, double w1y, double w1z, double w2x, double w2y, double w2z,
        double w3x, double w3y, double w3z,
        double u1x, double u1y, double u1z, double u2x, double u2y, double u2z,
        double u3x, double u3y, double u3z)
    {
        v1 = explicitPoint3D(v1x, v1y, v1z);
        v2 = explicitPoint3D(v2x, v2y, v2z);
        v3 = explicitPoint3D(v3x, v3y, v3z);
        w1 = explicitPoint3D(w1x, w1y, w1z);
        w2 = explicitPoint3D(w2x, w2y, w2z);
        w3 = explicitPoint3D(w3x, w3y, w3z);
        u1 = explicitPoint3D(u1x, u1y, u1z);
        u2 = explicitPoint3D(u2x, u2y, u2z);
        u3 = explicitPoint3D(u3x, u3y, u3z);
        point.emplace(v1, v2, v3, w1, w2, w3, u1, u2, u3);

        interval_number lx, ly, lz, d;
        return point->getIntervalLambda(lx, ly, lz, d);
    }

    bool init_exact(
        double v1x, double v1y, double v1z, double v2x, double v2y, double v2z,
        double v3x, double v3y, double v3z,
        double w1x, double w1y, double w1z, double w2x, double w2y, double w2z,
        double w3x, double w3y, double w3z,
        double u1x, double u1y, double u1z, double u2x, double u2y, double u2z,
        double u3x, double u3y, double u3z)
    {
        init(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
             w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
             u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z);
        bigrational ex, ey, ez;
        return point->getExactXYZCoordinates(ex, ey, ez);
    }
};

class TPI_exact_suppvars : public TPI_filtered_suppvars
{
};

// ---------------------------------------------------------------------------
// Two-stage interface: build the intersection point once, test it against many planes
// ---------------------------------------------------------------------------

inline bool orient3D_LPI_prefilter(
    const double& px, const double& py, const double& pz,
    const double& qx, const double& qy, const double& qz,
    const double& rx, const double& ry, const double& rz,
    const double& sx, const double& sy, const double& sz,
    const double& tx, const double& ty, const double& tz,
    LPI_filtered_suppvars& s)
{
    return s.init(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz);
}

// `px, py, pz` is the line's first point, already part of the cached intersection point,
// so it is accepted and ignored -- kept only to leave the call sites alone.
inline int orient3D_LPI_postfilter(
    const LPI_filtered_suppvars& s,
    const double&, const double&, const double&,
    const double& ax, const double& ay, const double& az,
    const double& bx, const double& by, const double& bz,
    const double& cx, const double& cy, const double& cz)
{
    // Interval stage only: returns 0 (UNCERTAIN) when it cannot decide and leaves the
    // caller to escalate, which is the contract this function has always had. It restores
    // the rounding mode itself on every exit path, and its INT_MAX overflow return is
    // guarded by `if constexpr (is_same<expansion, T>)`, so 0 is the only undecided value.
    return orient3d_indirect_IEEE_t<interval_number, interval_number>(
        *s.point, ax, ay, az, bx, by, bz, cx, cy, cz);
}

inline bool orient3D_TPI_prefilter(
    const double& v1x, const double& v1y, const double& v1z,
    const double& v2x, const double& v2y, const double& v2z,
    const double& v3x, const double& v3y, const double& v3z,
    const double& w1x, const double& w1y, const double& w1z,
    const double& w2x, const double& w2y, const double& w2z,
    const double& w3x, const double& w3y, const double& w3z,
    const double& u1x, const double& u1y, const double& u1z,
    const double& u2x, const double& u2y, const double& u2z,
    const double& u3x, const double& u3y, const double& u3z,
    TPI_filtered_suppvars& s)
{
    return s.init(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
                  w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
                  u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z);
}

inline int orient3D_TPI_postfilter(
    const TPI_filtered_suppvars& s,
    const double& q1x, const double& q1y, const double& q1z,
    const double& q2x, const double& q2y, const double& q2z,
    const double& q3x, const double& q3y, const double& q3z)
{
    // Interval stage only -- see orient3D_LPI_postfilter.
    return orient3d_indirect_IEEE_t<interval_number, interval_number>(
        *s.point, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
}

// Exact variants. Indirect_Predicates escalates to exact arithmetic by itself, so these
// differ from the filtered ones only in what the "pre" step reports: whether the point
// exists at all, rather than whether the interval filter settled it.

inline bool orient3D_LPI_pre_exact(
    double px, double py, double pz,
    double qx, double qy, double qz,
    double rx, double ry, double rz,
    double sx, double sy, double sz,
    double tx, double ty, double tz,
    LPI_exact_suppvars& s)
{
    return s.init_exact(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz);
}

inline int orient3D_LPI_post_exact(
    LPI_exact_suppvars& s,
    double, double, double,
    double ax, double ay, double az,
    double bx, double by, double bz,
    double cx, double cy, double cz)
{
    const explicitPoint3D a(ax, ay, az), b(bx, by, bz), c(cx, cy, cz);
    return genericPoint::orient3D(*s.point, a, b, c);
}

inline bool orient3D_TPI_pre_exact(
    double v1x, double v1y, double v1z, double v2x, double v2y, double v2z,
    double v3x, double v3y, double v3z,
    double w1x, double w1y, double w1z, double w2x, double w2y, double w2z,
    double w3x, double w3y, double w3z,
    double u1x, double u1y, double u1z, double u2x, double u2y, double u2z,
    double u3x, double u3y, double u3z,
    TPI_exact_suppvars& s)
{
    return s.init_exact(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
                        w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
                        u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z);
}

inline int orient3D_TPI_post_exact(
    TPI_exact_suppvars& s,
    double q1x, double q1y, double q1z,
    double q2x, double q2y, double q2z,
    double q3x, double q3y, double q3z)
{
    const explicitPoint3D q1(q1x, q1y, q1z), q2(q2x, q2y, q2z), q3(q3x, q3y, q3z);
    return genericPoint::orient3D(*s.point, q1, q2, q3);
}

// ---------------------------------------------------------------------------
// One-shot interface
// ---------------------------------------------------------------------------

inline int orient3D_LPI_filtered(
    double px, double py, double pz, double qx, double qy, double qz,
    double rx, double ry, double rz, double sx, double sy, double sz,
    double tx, double ty, double tz,
    double ax, double ay, double az, double bx, double by, double bz,
    double cx, double cy, double cz)
{
    LPI_filtered_suppvars s;
    if (!s.init(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz)) return UNCERTAIN;
    return orient3D_LPI_postfilter(s, px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz);
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
    TPI_filtered_suppvars s;
    if (!s.init(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
                w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
                u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z))
        return UNCERTAIN;
    return orient3D_TPI_postfilter(s, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
}

inline int orient3D_LPI_exact(
    double px, double py, double pz, double qx, double qy, double qz,
    double rx, double ry, double rz, double sx, double sy, double sz,
    double tx, double ty, double tz,
    double ax, double ay, double az, double bx, double by, double bz,
    double cx, double cy, double cz)
{
    LPI_exact_suppvars s;
    if (!s.init_exact(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, tx, ty, tz)) return 0;
    return orient3D_LPI_post_exact(s, px, py, pz, ax, ay, az, bx, by, bz, cx, cy, cz);
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
    TPI_exact_suppvars s;
    if (!s.init_exact(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
                      w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y, w3z,
                      u1x, u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z))
        return 0;
    return orient3D_TPI_post_exact(s, q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z);
}

#include "ip_filtered_ex.h"
