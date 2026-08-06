// Characterization/regression test for the AABB::is_triangle_cut_bounding_box
// hoist (src/AABB.cpp). It replicates the ORIGINAL per-node control flow and the
// REFACTORED control flow (query-triangle-invariant work precomputed once) side
// by side, and fuzzes them over millions of independent (triangle, box) queries.
//
// orient_2d is a pure function of its arguments in both the real code and here,
// so the cached "ori" signs and 2D projections equal the per-node ones bit-for-
// bit; proving the two control-flow structures agree on every query proves the
// hoist changes no envelope answer. The integer-coordinate regime exercises the
// ori==0 / o==0 degenerate branches that the double regime rarely hits.
//
// Self-contained (no library deps); exits non-zero on any mismatch.
#include <array>
#include <catch2/catch_all.hpp>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>


struct V2
{
    double x, y;
};
struct V3
{
    double v[3];
};

static inline int sgn(long double a)
{
    return (a > 0) - (a < 0);
}
static inline int orient_2d(const V2& p1, const V2& p2, const V2& p3)
{
    long double d = ((long double)p2.x - p1.x) * ((long double)p3.y - p1.y) -
                    ((long double)p2.y - p1.y) * ((long double)p3.x - p1.x);
    return sgn(d);
}
static inline V2 to_2d(const V3& p, int t)
{
    return V2{p.v[(t + 1) % 3], p.v[(t + 2) % 3]};
}
static inline void get_tri_corners(const V3& a, const V3& b, const V3& c, V3& mn, V3& mx)
{
    for (int k = 0; k < 3; k++) {
        mn.v[k] = std::min(std::min(a.v[k], b.v[k]), c.v[k]);
        mx.v[k] = std::max(std::max(a.v[k], b.v[k]), c.v[k]);
    }
}
static inline bool box_box(const V3& mn1, const V3& mx1, const V3& mn2, const V3& mx2)
{
    if (mx1.v[0] < mn2.v[0] || mx1.v[1] < mn2.v[1] || mx1.v[2] < mn2.v[2]) return false;
    if (mx2.v[0] < mn1.v[0] || mx2.v[1] < mn1.v[1] || mx2.v[2] < mn1.v[2]) return false;
    return true;
}

// ---- ORIGINAL structure (verbatim from the pre-refactor AABB.cpp) ----
static bool old_cut(const V3& tri0, const V3& tri1, const V3& tri2, const V3& bmin, const V3& bmax)
{
    V3 tmin, tmax;
    get_tri_corners(tri0, tri1, tri2, tmin, tmax);
    bool cut = box_box(tmin, tmax, bmin, bmax);
    if (!cut) return false;
    if (cut) {
        std::array<V2, 3> tri;
        std::array<V2, 4> mp;
        int o0, o1, o2, o3, ori;
        for (int i = 0; i < 3; i++) {
            tri[0] = to_2d(tri0, i);
            tri[1] = to_2d(tri1, i);
            tri[2] = to_2d(tri2, i);
            mp[0] = to_2d(bmin, i);
            mp[1] = to_2d(bmax, i);
            mp[2].x = mp[0].x;
            mp[2].y = mp[1].y;
            mp[3].x = mp[1].x;
            mp[3].y = mp[0].y;
            for (int j = 0; j < 3; j++) {
                o0 = orient_2d(mp[0], tri[j % 3], tri[(j + 1) % 3]);
                o1 = orient_2d(mp[1], tri[j % 3], tri[(j + 1) % 3]);
                o2 = orient_2d(mp[2], tri[j % 3], tri[(j + 1) % 3]);
                o3 = orient_2d(mp[3], tri[j % 3], tri[(j + 1) % 3]);
                ori = orient_2d(tri[(j + 2) % 3], tri[j % 3], tri[(j + 1) % 3]);
                if (ori == 0) continue;
                if (ori * o0 <= 0 && ori * o1 <= 0 && ori * o2 <= 0 && ori * o3 <= 0) return false;
            }
        }
    }
    return cut;
}

// ---- REFACTORED structure (verbatim from the new AABB.h/.cpp) ----
struct Cache
{
    V3 tmin, tmax;
    std::array<std::array<V2, 3>, 3> tri2d;
    std::array<std::array<int, 3>, 3> ori;
};
static Cache build(const V3& t0, const V3& t1, const V3& t2)
{
    Cache c;
    get_tri_corners(t0, t1, t2, c.tmin, c.tmax);
    for (int i = 0; i < 3; i++) {
        c.tri2d[i][0] = to_2d(t0, i);
        c.tri2d[i][1] = to_2d(t1, i);
        c.tri2d[i][2] = to_2d(t2, i);
        for (int j = 0; j < 3; j++)
            c.ori[i][j] =
                orient_2d(c.tri2d[i][(j + 2) % 3], c.tri2d[i][j % 3], c.tri2d[i][(j + 1) % 3]);
    }
    return c;
}
static bool new_cut(const Cache& tc, const V3& bmin, const V3& bmax)
{
    if (!box_box(tc.tmin, tc.tmax, bmin, bmax)) return false;
    std::array<V2, 4> mp;
    for (int i = 0; i < 3; i++) {
        const std::array<V2, 3>& tri = tc.tri2d[i];
        mp[0] = to_2d(bmin, i);
        mp[1] = to_2d(bmax, i);
        mp[2].x = mp[0].x;
        mp[2].y = mp[1].y;
        mp[3].x = mp[1].x;
        mp[3].y = mp[0].y;
        for (int j = 0; j < 3; j++) {
            const int ori = tc.ori[i][j];
            if (ori == 0) continue;
            const int o0 = orient_2d(mp[0], tri[j % 3], tri[(j + 1) % 3]);
            const int o1 = orient_2d(mp[1], tri[j % 3], tri[(j + 1) % 3]);
            const int o2 = orient_2d(mp[2], tri[j % 3], tri[(j + 1) % 3]);
            const int o3 = orient_2d(mp[3], tri[j % 3], tri[(j + 1) % 3]);
            if (ori * o0 <= 0 && ori * o1 <= 0 && ori * o2 <= 0 && ori * o3 <= 0) return false;
        }
    }
    return true;
}

TEST_CASE("box cut equivalence", "[equivalence]")
{
    std::mt19937_64 rng(0xC0FFEE);
    auto run = [&](bool ints, long n) -> long {
        long mismatch = 0, trues = 0;
        std::uniform_int_distribution<int> di(-6, 6);
        std::uniform_real_distribution<double> dr(-3, 3);
        for (long k = 0; k < n; k++) {
            V3 t0, t1, t2, ba, bb, bmin, bmax;
            auto g = [&](V3& p) {
                for (int c = 0; c < 3; c++) p.v[c] = ints ? (double)di(rng) : dr(rng);
            };
            g(t0);
            g(t1);
            g(t2);
            g(ba);
            g(bb);
            for (int c = 0; c < 3; c++) {
                bmin.v[c] = std::min(ba.v[c], bb.v[c]);
                bmax.v[c] = std::max(ba.v[c], bb.v[c]);
            }
            bool a = old_cut(t0, t1, t2, bmin, bmax);
            bool b = new_cut(build(t0, t1, t2), bmin, bmax);
            if (a != b) {
                if (++mismatch <= 5) printf("MISMATCH ints=%d old=%d new=%d\n", ints, a, b);
            }
            if (a) trues++;
        }
        printf(
            "regime ints=%d: n=%ld mismatches=%ld (accept-rate=%.1f%%)\n",
            ints,
            n,
            mismatch,
            100.0 * trues / n);
        return mismatch;
    };
    long m = 0;
    m += run(true, 2000000);
    m += run(false, 2000000);
    printf(m == 0 ? "ALL EQUIVALENT\n" : "FAILED: %ld mismatches\n", m);
    REQUIRE(m == 0);
    // return m == 0 ? 0 : 1;
}
