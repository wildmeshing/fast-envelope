#pragma once

#include "Types.hpp"

#include <array>
#include <vector>

class genericPoint;

namespace fastEnvelope {

/// Conservative two-dimensional envelope around a collection of line segments.
///
/// Each input edge is expanded into an edge-aligned rectangle. A point or segment
/// is inside the envelope when it is covered by the union of these rectangles.
class FastEnvelope2D
{
public:
    using HalfPlane = std::array<Vector2, 2>;

    FastEnvelope2D() = default;
    FastEnvelope2D(
        const std::vector<Vector2>& vertices,
        const std::vector<Vector2i>& edges,
        Scalar epsilon);

    void
    init(const std::vector<Vector2>& vertices, const std::vector<Vector2i>& edges, Scalar epsilon);

    void init(
        const std::vector<Vector2>& vertices,
        const std::vector<Vector2i>& edges,
        const std::vector<Scalar>& epsilons);

    bool is_outside(const Vector2& point) const;
    bool is_outside(const Vector2& point0, const Vector2& point1) const;

private:
    struct EdgeEnvelope
    {
        std::array<HalfPlane, 4> halfplanes;
        std::array<Vector2, 2> bounds;
    };

    static EdgeEnvelope
    make_edge_envelope(const Vector2& point0, const Vector2& point1, Scalar epsilon);
    static bool contains(const EdgeEnvelope& envelope, const Vector2& point);
    static bool contains(const EdgeEnvelope& envelope, const genericPoint& point);

    std::vector<EdgeEnvelope> envelopes_;
};

} // namespace fastEnvelope
