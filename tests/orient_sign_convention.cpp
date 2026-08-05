// Regression test for the sign convention of the indirect orientation predicates in
// indirectPredicates/ip_filtered.h.
//
// FastEnvelope does not read these predicates as a bare "signed volume of the tetrahedron
// (i, a, b, c)". It reads them as the HEIGHT of the implicit point i over the oriented
// plane <a, b, c>: positive above, negative below, zero on it. Its prism-containment tests
// are written directly against that reading -- a point is inside a prism when it is below
// every one of the prism's outward-oriented faces, which the callers spell as every face
// returning NEGATIVE (see the `tot == halfspace[prismindex[i]].size()` tests in
// FastEnvelope.cpp).
//
// This matters because ip_filtered.h is now a thin adapter over MarcoAttene's
// Indirect_Predicates, and `genericPoint::orient3D` there returns the OPPOSITE sign. If the
// adapter forwards it unchanged, containment fails everywhere: every prism is skipped,
// every point looks like it is outside, and the envelope silently becomes far more
// conservative than the requested epsilon -- meshes grow, Hausdorff distance collapses and
// runtime balloons, with no error anywhere.
//
// The geometry below is chosen so the intersection point is exactly the origin and the
// answers are obvious by inspection, so this test says nothing about the arithmetic; it
// pins only the sign convention, in both the filtered and the exact stage.
//
// Exits non-zero on any mismatch.
#include <cstdio>

#include <indirectpredicates/ip_filtered.h>

static int failures = 0;

static void check(const char* what, int got, int want)
{
    const bool ok = (got == want);
    if (!ok) ++failures;
    std::printf("%-58s got %+d  want %+d  %s\n", what, got, want, ok ? "ok" : "FAILED");
}

int main()
{
    initFPU();

    // A plane at height h whose (b-a) x (c-a) normal points along +z, so "above the plane"
    // means z > h. Passed to the predicates as <a, b, c>.
    struct Plane
    {
        double a[3], b[3], c[3];
    };
    auto plane_up = [](double h) {
        return Plane{{0, 0, h}, {1, 0, h}, {0, 1, h}};
    };
    // Same plane with b and c swapped: the normal now points along -z, so the sign of the
    // height of any point off the plane flips.
    auto plane_down = [](double h) {
        return Plane{{0, 0, h}, {0, 1, h}, {1, 0, h}};
    };

    const Plane below = plane_up(-1); // origin is ABOVE it  -> +1
    const Plane above = plane_up(+1); // origin is BELOW it  -> -1
    const Plane below_flipped = plane_down(-1); //            -> -1

    // ---- LPI: the z axis meets the plane z = 0 at the origin -------------------------
    {
        const double p[3] = {0, 0, -2}, q[3] = {0, 0, 3};              // the line
        const double r[3] = {5, 0, 0}, s[3] = {0, 7, 0}, t[3] = {0, 0, 0}; // the plane z = 0

        LPI_filtered_suppvars sv;
        if (!orient3D_LPI_prefilter(
                p[0], p[1], p[2], q[0], q[1], q[2],
                r[0], r[1], r[2], s[0], s[1], s[2], t[0], t[1], t[2], sv)) {
            std::printf("LPI prefilter unexpectedly failed on exact input\n");
            return 1;
        }
        auto post = [&](const Plane& pl) {
            return orient3D_LPI_postfilter(
                sv, p[0], p[1], p[2],
                pl.a[0], pl.a[1], pl.a[2], pl.b[0], pl.b[1], pl.b[2],
                pl.c[0], pl.c[1], pl.c[2]);
        };
        check("LPI postfilter, point above the plane", post(below), 1);
        check("LPI postfilter, point below the plane", post(above), -1);
        check("LPI postfilter, plane orientation flipped", post(below_flipped), -1);

        LPI_exact_suppvars ex;
        if (!orient3D_LPI_pre_exact(
                p[0], p[1], p[2], q[0], q[1], q[2],
                r[0], r[1], r[2], s[0], s[1], s[2], t[0], t[1], t[2], ex)) {
            std::printf("LPI pre_exact unexpectedly failed on exact input\n");
            return 1;
        }
        auto post_ex = [&](const Plane& pl) {
            return orient3D_LPI_post_exact(
                ex, p[0], p[1], p[2],
                pl.a[0], pl.a[1], pl.a[2], pl.b[0], pl.b[1], pl.b[2],
                pl.c[0], pl.c[1], pl.c[2]);
        };
        check("LPI post_exact, point above the plane", post_ex(below), 1);
        check("LPI post_exact, point below the plane", post_ex(above), -1);
        check("LPI post_exact, plane orientation flipped", post_ex(below_flipped), -1);
    }

    // ---- TPI: the planes x = 0, y = 0 and z = 0 meet at the origin -------------------
    {
        // Each plane is given by three of its points.
        const double v1[3] = {0, 0, 0}, v2[3] = {0, 3, 0}, v3[3] = {0, 0, 5}; // x = 0
        const double w1[3] = {0, 0, 0}, w2[3] = {7, 0, 0}, w3[3] = {0, 0, 2}; // y = 0
        const double u1[3] = {0, 0, 0}, u2[3] = {4, 0, 0}, u3[3] = {0, 6, 0}; // z = 0

        TPI_filtered_suppvars sv;
        if (!orient3D_TPI_prefilter(
                v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2],
                w1[0], w1[1], w1[2], w2[0], w2[1], w2[2], w3[0], w3[1], w3[2],
                u1[0], u1[1], u1[2], u2[0], u2[1], u2[2], u3[0], u3[1], u3[2], sv)) {
            std::printf("TPI prefilter unexpectedly failed on exact input\n");
            return 1;
        }
        auto post = [&](const Plane& pl) {
            return orient3D_TPI_postfilter(
                sv, pl.a[0], pl.a[1], pl.a[2], pl.b[0], pl.b[1], pl.b[2],
                pl.c[0], pl.c[1], pl.c[2]);
        };
        check("TPI postfilter, point above the plane", post(below), 1);
        check("TPI postfilter, point below the plane", post(above), -1);
        check("TPI postfilter, plane orientation flipped", post(below_flipped), -1);

        TPI_exact_suppvars ex;
        if (!orient3D_TPI_pre_exact(
                v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2],
                w1[0], w1[1], w1[2], w2[0], w2[1], w2[2], w3[0], w3[1], w3[2],
                u1[0], u1[1], u1[2], u2[0], u2[1], u2[2], u3[0], u3[1], u3[2], ex)) {
            std::printf("TPI pre_exact unexpectedly failed on exact input\n");
            return 1;
        }
        auto post_ex = [&](const Plane& pl) {
            return orient3D_TPI_post_exact(
                ex, pl.a[0], pl.a[1], pl.a[2], pl.b[0], pl.b[1], pl.b[2],
                pl.c[0], pl.c[1], pl.c[2]);
        };
        check("TPI post_exact, point above the plane", post_ex(below), 1);
        check("TPI post_exact, point below the plane", post_ex(above), -1);
        check("TPI post_exact, plane orientation flipped", post_ex(below_flipped), -1);
    }

    std::printf("\n%s\n", failures ? "FAILED" : "all sign conventions as expected");
    return failures ? 1 : 0;
}
