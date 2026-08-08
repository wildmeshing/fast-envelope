#include <fastenvelope/FastEnvelope2D.h>

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <random>
#include <vector>

using namespace fastEnvelope;

TEST_CASE("2D edge envelopes classify points and segments", "[envelope-2d]")
{
    const std::vector<Vector2> vertices = {
        Vector2(0, 0),
        Vector2(3, 0),
    };
    const std::vector<Vector2i> edges = {Vector2i(0, 1)};
    const Scalar epsilon = std::sqrt(Scalar(2));

    const FastEnvelope2D envelope(vertices, edges, epsilon);

    CHECK_FALSE(envelope.is_outside(Vector2(1.5, 0)));
    CHECK_FALSE(envelope.is_outside(Vector2(-0.5, 0.5)));
    CHECK(envelope.is_outside(Vector2(-1.1, 0)));
    CHECK(envelope.is_outside(Vector2(1.5, 1.1)));

    CHECK_FALSE(envelope.is_outside(Vector2(-0.5, 0.5), Vector2(3.5, -0.5)));
    CHECK(envelope.is_outside(Vector2(-0.5, 0), Vector2(4.1, 0)));
}

TEST_CASE("2D rectangles align with arbitrary edges", "[envelope-2d][geometry]")
{
    std::mt19937_64 random(0x2DE6E);
    std::uniform_real_distribution<Scalar> coordinate(-10, 10);
    std::uniform_real_distribution<Scalar> epsilon_distribution(1e-3, 2);

    for (int sample = 0; sample < 250; ++sample) {
        const Vector2 point0(coordinate(random), coordinate(random));
        Vector2 edge(coordinate(random), coordinate(random));
        if (edge.squaredNorm() < Scalar(1e-6)) edge[0] += 1;
        const Vector2 point1 = point0 + edge;
        const Scalar epsilon = epsilon_distribution(random);
        const Scalar tolerance = epsilon / std::sqrt(Scalar(2));
        const Vector2 direction = edge.normalized();
        const Vector2 normal(-direction[1], direction[0]);
        const Vector2 midpoint = (point0 + point1) / 2;

        const FastEnvelope2D envelope({point0, point1}, {Vector2i(0, 1)}, epsilon);

        CHECK_FALSE(envelope.is_outside(midpoint + Scalar(0.9) * tolerance * normal));
        CHECK(envelope.is_outside(midpoint + Scalar(1.1) * tolerance * normal));
        CHECK_FALSE(envelope.is_outside(point0 - Scalar(0.9) * tolerance * direction));
        CHECK(envelope.is_outside(point0 - Scalar(1.1) * tolerance * direction));
        CHECK_FALSE(envelope.is_outside(point1 + Scalar(0.9) * tolerance * (direction + normal)));
        CHECK(envelope.is_outside(point1 + Scalar(1.1) * tolerance * (direction + normal)));
    }
}

TEST_CASE("2D envelopes support zero-length edges", "[envelope-2d]")
{
    const FastEnvelope2D envelope({Vector2(2, 3)}, {Vector2i(0, 0)}, std::sqrt(Scalar(2)));

    CHECK_FALSE(envelope.is_outside(Vector2(2, 3)));
    CHECK_FALSE(envelope.is_outside(Vector2(2.5, 3.5)));
    CHECK(envelope.is_outside(Vector2(3.1, 3)));
}

TEST_CASE("adaptive 2D envelope widths stay attached to their edges", "[envelope-2d]")
{
    const std::vector<Vector2> vertices = {
        Vector2(10, 0),
        Vector2(13, 0),
        Vector2(0, 0),
        Vector2(3, 0),
    };
    const std::vector<Vector2i> edges = {
        Vector2i(0, 1),
        Vector2i(2, 3),
    };
    const std::vector<Scalar> epsilons = {
        std::sqrt(Scalar(2)),
        Scalar(0.1) * std::sqrt(Scalar(2)),
    };

    FastEnvelope2D envelope;
    envelope.init(vertices, edges, epsilons);

    CHECK_FALSE(envelope.is_outside(Vector2(11.5, 0.5)));
    CHECK(envelope.is_outside(Vector2(1.5, 0.5)));
    CHECK_FALSE(envelope.is_outside(Vector2(1.5, 0)));
}

TEST_CASE("2D segment queries require continuous union coverage", "[envelope-2d]")
{
    SECTION("overlapping edge rectangles cover the whole query")
    {
        const FastEnvelope2D envelope(
            {Vector2(0, 0), Vector2(1, 0), Vector2(2, 0)},
            {Vector2i(0, 1), Vector2i(1, 2)},
            Scalar(0.2) * std::sqrt(Scalar(2)));

        CHECK_FALSE(envelope.is_outside(Vector2(-0.1, 0.1), Vector2(2.1, 0.1)));
    }

    SECTION("a gap between rectangles leaves the query outside")
    {
        const FastEnvelope2D envelope(
            {Vector2(0, 0), Vector2(1, 0), Vector2(3, 0), Vector2(4, 0)},
            {Vector2i(0, 1), Vector2i(2, 3)},
            Scalar(0.2) * std::sqrt(Scalar(2)));

        CHECK(envelope.is_outside(Vector2(0.5, 0), Vector2(3.5, 0)));
    }

    SECTION("coverage of both endpoints is not sufficient")
    {
        const FastEnvelope2D envelope(
            {Vector2(0, 0), Vector2(1, 0), Vector2(1, 0), Vector2(1, 1)},
            {Vector2i(0, 1), Vector2i(2, 3)},
            Scalar(0.1) * std::sqrt(Scalar(2)));

        CHECK(envelope.is_outside(Vector2(0, 0), Vector2(1, 1)));
    }
}

TEST_CASE("2D segment clipping does not bridge a microscopic gap", "[envelope-2d][exact]")
{
    const Scalar tolerance = Scalar(0.25);
    const Scalar epsilon = tolerance * std::sqrt(Scalar(2));

    SECTION("boxes that touch hand coverage over at the shared clipped point")
    {
        const FastEnvelope2D envelope(
            {Vector2(0, 0), Vector2(1, 0), Vector2(1.5, 0), Vector2(2.5, 0)},
            {Vector2i(0, 1), Vector2i(2, 3)},
            epsilon);

        CHECK_FALSE(envelope.is_outside(Vector2(0, 0), Vector2(2.5, 0)));
    }

    SECTION("boxes can hand coverage over at a shared corner")
    {
        const FastEnvelope2D envelope(
            {Vector2(0, 0), Vector2(2, 2)},
            {Vector2i(0, 0), Vector2i(1, 1)},
            std::sqrt(Scalar(2)));

        CHECK_FALSE(envelope.is_outside(Vector2(-0.5, -0.5), Vector2(2.5, 2.5)));
    }

    SECTION("the query may overlap boundary edges and transition at their vertices")
    {
        const FastEnvelope2D envelope(
            {Vector2(0, 0), Vector2(2, 2)},
            {Vector2i(0, 0), Vector2i(1, 1)},
            std::sqrt(Scalar(2)));

        CHECK_FALSE(envelope.is_outside(Vector2(-0.5, 1), Vector2(2.5, 1)));
    }

    SECTION("a box may intersect the query in only the transition point")
    {
        // The middle square touches y = x only at (1, 1). All three closed
        // query-box intervals still belong to the same component there.
        const FastEnvelope2D envelope(
            {Vector2(0, 0), Vector2(0, 2), Vector2(2, 2)},
            {Vector2i(0, 0), Vector2i(1, 1), Vector2i(2, 2)},
            std::sqrt(Scalar(2)));

        CHECK_FALSE(envelope.is_outside(Vector2(-0.5, -0.5), Vector2(2.5, 2.5)));
    }

    SECTION("distinct boxes remain disconnected below the old merge tolerance")
    {
        const Scalar gap = Scalar(1e-14);
        const FastEnvelope2D envelope(
            {Vector2(0, 0), Vector2(1, 0), Vector2(1.5 + gap, 0), Vector2(2.5, 0)},
            {Vector2i(0, 1), Vector2i(2, 3)},
            epsilon);

        CHECK(envelope.is_outside(Vector2(0, 0), Vector2(2.5, 0)));
    }
}

TEST_CASE("2D segment traversal starts from the complete endpoint component", "[envelope-2d]")
{
    SECTION("the first box may only touch the query at its first endpoint")
    {
        // The first square is [-1, 1]^2. The second is [1, 3] x [0, 2].
        // Both contain query0, but only the second one covers the rest of the query.
        const FastEnvelope2D envelope(
            {Vector2(0, 0), Vector2(2, 1)},
            {Vector2i(0, 0), Vector2i(1, 1)},
            std::sqrt(Scalar(2)));

        CHECK_FALSE(envelope.is_outside(Vector2(1, 1), Vector2(2, 1)));
    }

    SECTION("nested starting boxes connect to an arbitrary box chain")
    {
        const Scalar root_two = std::sqrt(Scalar(2));
        FastEnvelope2D envelope;
        envelope.init(
            {
                Vector2(0, 0),
                Vector2(0, 0),
                Vector2(0, 0),
                Vector2(0, 0),
                Vector2(3, 0),
                Vector2(5, 0),
            },
            {
                Vector2i(0, 0),
                Vector2i(1, 1),
                Vector2i(2, 2),
                Vector2i(3, 3),
                Vector2i(4, 4),
                Vector2i(5, 5),
            },
            {
                Scalar(0.25) * root_two,
                Scalar(0.5) * root_two,
                root_two,
                Scalar(2) * root_two,
                root_two,
                root_two,
            });

        CHECK_FALSE(envelope.is_outside(Vector2(0, 0), Vector2(5.5, 0)));
    }
}

TEST_CASE("2D queries handle every collinear boundary configuration", "[envelope-2d][exact]")
{
    const Scalar epsilon = std::sqrt(Scalar(2));

    SECTION("a query spans partial and complete overlaps with a boundary chain")
    {
        const FastEnvelope2D envelope(
            {Vector2(-2, 0), Vector2(0, 0), Vector2(2, 0)},
            {Vector2i(0, 0), Vector2i(1, 1), Vector2i(2, 2)},
            epsilon);

        CHECK_FALSE(envelope.is_outside(Vector2(-2.5, 1), Vector2(2.5, 1)));
        CHECK_FALSE(envelope.is_outside(Vector2(2.5, 1), Vector2(-2.5, 1)));
    }

    SECTION("the query lies entirely inside one boundary side")
    {
        const FastEnvelope2D envelope({Vector2(0, 0)}, {Vector2i(0, 0)}, epsilon);

        CHECK_FALSE(envelope.is_outside(Vector2(-0.5, 1), Vector2(0.5, 1)));
    }

    SECTION("a partial overlap without continued coverage is outside")
    {
        const FastEnvelope2D envelope({Vector2(0, 0)}, {Vector2i(0, 0)}, epsilon);

        CHECK(envelope.is_outside(Vector2(0, 1), Vector2(2, 1)));
    }

    SECTION("touching one boundary endpoint without continued coverage is outside")
    {
        const FastEnvelope2D envelope({Vector2(0, 0)}, {Vector2i(0, 0)}, epsilon);

        CHECK(envelope.is_outside(Vector2(1, 1), Vector2(2, 1)));
    }

    SECTION("disjoint collinear segments remain outside")
    {
        const FastEnvelope2D envelope({Vector2(0, 0)}, {Vector2i(0, 0)}, epsilon);

        CHECK(envelope.is_outside(Vector2(2, 1), Vector2(3, 1)));
    }
}

TEST_CASE(
    "2D AABB queries filter with the segment instead of its bounding box",
    "[envelope-2d][aabb]")
{
    std::vector<std::array<Vector3, 2>> boxes(4);
    boxes[0][0] = Vector3(0, 0, 0);
    boxes[0][1] = Vector3(1, 1, 0);
    boxes[1][0] = Vector3(0, 2, 0);
    boxes[1][1] = Vector3(1, 3, 0);
    boxes[2][0] = Vector3(2, 0, 0);
    boxes[2][1] = Vector3(3, 1, 0);
    boxes[3][0] = Vector3(2, 2, 0);
    boxes[3][1] = Vector3(3, 3, 0);

    AABB tree;
    tree.init(boxes);
    const auto find = [&](const Vector2& point0, const Vector2& point1) {
        std::vector<unsigned int> candidates;
        tree.segment_find_bbox(point0, point1, candidates);
        return candidates;
    };

    CHECK(find(Vector2(-1, -1), Vector2(4, 4)) == std::vector<unsigned int>{0, 3});
    CHECK(find(Vector2(4, 4), Vector2(-1, -1)) == std::vector<unsigned int>{0, 3});
    CHECK(find(Vector2(-1, 1), Vector2(1.5, 1)) == std::vector<unsigned int>{0});
    CHECK(find(Vector2(-1, 1), Vector2(0, 0)) == std::vector<unsigned int>{0});
}

TEST_CASE("2D segment traversal indexes its visited set by candidate", "[envelope-2d]")
{
    // The traversal marks boxes visited by their position in the AABB's candidate list, not
    // by their envelope id, so that a query costs O(candidates) rather than O(#input edges).
    // Morton resorting is not applied to FastEnvelope2D, but the tree still returns candidates
    // in tree order and only for the boxes the query's bounding box touches -- so on a long
    // polyline the two indexings diverge, and confusing them silently breaks the chain.
    //
    // The chain here is a staircase of 64 overlapping unit rectangles, with the query crossing
    // the last 32. Under an id-indexed visited set this passes either way; it fails only if
    // the two index spaces are mixed up.
    const int count = 64;
    const Scalar root_two = std::sqrt(Scalar(2));
    std::vector<Vector2> vertices;
    std::vector<Vector2i> edges;
    for (int i = 0; i < count; ++i) {
        vertices.emplace_back(Scalar(i), Scalar(0));
        vertices.emplace_back(Scalar(i) + Scalar(0.75), Scalar(0));
        edges.emplace_back(2 * i, 2 * i + 1);
    }

    const FastEnvelope2D envelope(vertices, edges, Scalar(0.2) * root_two);

    // Fully covered: consecutive rectangles overlap, so the union is connected end to end.
    CHECK_FALSE(envelope.is_outside(Vector2(32, 0), Vector2(Scalar(count) - 1, 0)));
    CHECK_FALSE(envelope.is_outside(Vector2(0, 0), Vector2(Scalar(count) - 1, 0)));
    // Reversed, so the traversal starts from the other end of the same chain.
    CHECK_FALSE(envelope.is_outside(Vector2(Scalar(count) - 1, 0), Vector2(0, 0)));
    // Past the last rectangle there is nothing left to hand coverage over to.
    CHECK(envelope.is_outside(Vector2(32, 0), Vector2(Scalar(count) + 1, 0)));
    // Off the axis by more than the half-width.
    CHECK(envelope.is_outside(Vector2(32, Scalar(0.5)), Vector2(40, Scalar(0.5))));
}
