#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/common_algorithms.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <random>
#include <vector>

using namespace fastEnvelope;

TEST_CASE("edge halfspaces form an oriented segment box", "[edge-envelope]")
{
    const std::vector<Vector3> vertices = {
        Vector3(0, 0, 0),
        Vector3(3, 0, 0),
    };
    const std::vector<Vector2i> edges = {Vector2i(0, 1)};

    std::vector<std::vector<std::array<Vector3, 3>>> halfspaces;
    std::vector<std::array<Vector3, 2>> corners;
    algorithms::halfspace_generation(
        vertices,
        edges,
        halfspaces,
        corners,
        std::sqrt(Scalar(3)));

    REQUIRE(halfspaces.size() == 1);
    REQUIRE(halfspaces.front().size() == 6);
    REQUIRE(corners.size() == 1);

    const Vector3 midpoint(1.5, 0, 0);
    for (const auto& plane : halfspaces.front()) {
        CHECK(algorithms::orient_3d(plane[0], plane[1], plane[2], midpoint) == 1);
    }
}

TEST_CASE("edge hexes stay aligned with arbitrary edges", "[edge-envelope][geometry]")
{
    std::mt19937_64 random(0xED6E);
    std::uniform_real_distribution<Scalar> coordinate(-10, 10);
    std::uniform_real_distribution<Scalar> epsilon_distribution(1e-3, 2);

    for (int sample = 0; sample < 250; ++sample) {
        const Vector3 p0(coordinate(random), coordinate(random), coordinate(random));
        Vector3 delta(coordinate(random), coordinate(random), coordinate(random));
        if (delta.squaredNorm() < Scalar(1e-6)) delta[0] += 1;
        const Vector3 p1 = p0 + delta;
        const Scalar epsilon = epsilon_distribution(random);
        const Scalar tolerance = epsilon / std::sqrt(Scalar(3));

        const std::vector<Vector3> vertices = {p0, p1};
        const std::vector<Vector2i> edges = {Vector2i(0, 1)};
        std::vector<std::vector<std::array<Vector3, 3>>> halfspaces;
        std::vector<std::array<Vector3, 2>> corners;
        algorithms::halfspace_generation(
            vertices,
            edges,
            halfspaces,
            corners,
            epsilon);

        REQUIRE(halfspaces.size() == 1);
        REQUIRE(halfspaces.front().size() == 6);

        const Vector3 edge_direction = delta.normalized();
        const Vector3 midpoint = (p0 + p1) / 2;
        int end_caps = 0;
        int side_faces = 0;

        for (const auto& plane : halfspaces.front()) {
            Vector3 normal =
                (plane[1] - plane[0]).cross(plane[2] - plane[0]);
            REQUIRE(normal.squaredNorm() > 0);
            normal.normalize();

            const Scalar alignment = std::abs(normal.dot(edge_direction));
            const Scalar distance = std::abs(normal.dot(midpoint - plane[0]));
            if (alignment > Scalar(1) - Scalar(1e-10)) {
                ++end_caps;
                CHECK(
                    distance
                    == Catch::Approx(delta.norm() / 2 + tolerance)
                           .epsilon(1e-10)
                           .margin(1e-10));
            } else {
                ++side_faces;
                CHECK(alignment < Scalar(1e-10));
                CHECK(
                    distance
                    == Catch::Approx(tolerance).epsilon(1e-10).margin(1e-10));
            }

            CHECK(algorithms::orient_3d(plane[0], plane[1], plane[2], midpoint) == 1);
            for (const Vector3& vertex : plane) {
                for (int axis = 0; axis < 3; ++axis) {
                    CHECK(vertex[axis] >= corners.front()[0][axis]);
                    CHECK(vertex[axis] <= corners.front()[1][axis]);
                }
            }
        }

        CHECK(end_caps == 2);
        CHECK(side_faces == 4);
    }
}

TEST_CASE("constant-width edge envelopes classify points", "[edge-envelope]")
{
    const std::vector<Vector3> vertices = {
        Vector3(0, 0, 0),
        Vector3(3, 0, 0),
    };
    const std::vector<Vector2i> edges = {Vector2i(0, 1)};

    FastEnvelope envelope;
    envelope.init(vertices, edges, std::sqrt(Scalar(3)));

    CHECK_FALSE(envelope.is_outside(Vector3(1.5, 0, 0)));
    CHECK_FALSE(envelope.is_outside(Vector3(-0.5, 0, 0)));
    CHECK(envelope.is_outside(Vector3(-1.1, 0, 0)));
    CHECK(envelope.is_outside(Vector3(1.5, 2, 0)));
}

TEST_CASE("edge envelopes support zero-length edges", "[edge-envelope]")
{
    const std::vector<Vector3> vertices = {Vector3(2, 3, 4)};
    const std::vector<Vector2i> edges = {Vector2i(0, 0)};

    FastEnvelope envelope;
    envelope.init(vertices, edges, std::sqrt(Scalar(3)));

    CHECK_FALSE(envelope.is_outside(Vector3(2, 3, 4)));
    CHECK_FALSE(envelope.is_outside(Vector3(2.5, 3.5, 4.5)));
    CHECK(envelope.is_outside(Vector3(3.1, 3, 4)));
}

TEST_CASE("adaptive edge widths follow edges through Morton sorting", "[edge-envelope]")
{
    const std::vector<Vector3> vertices = {
        Vector3(10, 0, 0),
        Vector3(13, 0, 0),
        Vector3(0, 0, 0),
        Vector3(3, 0, 0),
    };
    const std::vector<Vector2i> edges = {
        Vector2i(0, 1),
        Vector2i(2, 3),
    };
    const std::vector<Scalar> epsilons = {
        std::sqrt(Scalar(3)),
        Scalar(0.1) * std::sqrt(Scalar(3)),
    };

    FastEnvelope envelope;
    envelope.init(vertices, edges, epsilons);

    CHECK_FALSE(envelope.is_outside(Vector3(11.5, 0.5, 0)));
    CHECK(envelope.is_outside(Vector3(1.5, 0.5, 0)));
    CHECK_FALSE(envelope.is_outside(Vector3(1.5, 0, 0)));
}
