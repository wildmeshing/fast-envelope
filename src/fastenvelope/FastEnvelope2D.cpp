#include "FastEnvelope2D.h"

#include "common_algorithms.h"
#include "indirectPredicates/ip_filtered.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

namespace fastEnvelope {
namespace {

Vector3 lift_to_3d(const Vector2& point)
{
    return Vector3(point[0], point[1], 0);
}

bool point_in_segment_bounds(const Vector2& point, const Vector2& segment0, const Vector2& segment1)
{
    return point[0] >= std::min(segment0[0], segment1[0]) &&
           point[0] <= std::max(segment0[0], segment1[0]) &&
           point[1] >= std::min(segment0[1], segment1[1]) &&
           point[1] <= std::max(segment0[1], segment1[1]);
}

} // namespace

FastEnvelope2D::FastEnvelope2D(
    const std::vector<Vector2>& vertices,
    const std::vector<Vector2i>& edges,
    Scalar epsilon)
{
    init(vertices, edges, epsilon);
}

void FastEnvelope2D::init(
    const std::vector<Vector2>& vertices,
    const std::vector<Vector2i>& edges,
    Scalar epsilon)
{
    init(vertices, edges, std::vector<Scalar>(edges.size(), epsilon));
}

void FastEnvelope2D::init(
    const std::vector<Vector2>& vertices,
    const std::vector<Vector2i>& edges,
    const std::vector<Scalar>& epsilons)
{
    assert(epsilons.size() == edges.size());

    envelopes_.clear();
    envelopes_.reserve(edges.size());
    std::vector<std::array<Vector3, 2>> boxes;
    boxes.reserve(edges.size());
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const int first = edges[i][0];
        const int second = edges[i][1];
        assert(first >= 0 && second >= 0);
        assert(first < vertices.size());
        assert(second < vertices.size());

        envelopes_.push_back(make_edge_envelope(vertices[first], vertices[second], epsilons[i]));
        const EdgeEnvelope& envelope = envelopes_.back();
        std::array<Vector3, 2> box;
        box[0] = lift_to_3d(envelope.bounds[0]);
        box[1] = lift_to_3d(envelope.bounds[1]);
        boxes.push_back(box);
    }

    if (!boxes.empty()) tree_.init(boxes);

    initFPU();
}

FastEnvelope2D::EdgeEnvelope
FastEnvelope2D::make_edge_envelope(const Vector2& point0, const Vector2& point1, Scalar epsilon)
{
    assert(std::isfinite(epsilon) && epsilon > 0);

    const Scalar tolerance = epsilon / std::sqrt(Scalar(2));
    const Vector2 edge = point1 - point0;
    std::array<Vector2, 4> corners;

    if (edge[0] == 0 && edge[1] == 0) {
        corners[0] = point0 + Vector2(-tolerance, -tolerance);
        corners[1] = point0 + Vector2(tolerance, -tolerance);
        corners[2] = point0 + Vector2(tolerance, tolerance);
        corners[3] = point0 + Vector2(-tolerance, tolerance);
    } else {
        const Vector2 direction = edge.normalized();
        const Vector2 normal(-direction[1], direction[0]);
        corners[0] = point0 - tolerance * direction - tolerance * normal;
        corners[1] = point1 + tolerance * direction - tolerance * normal;
        corners[2] = point1 + tolerance * direction + tolerance * normal;
        corners[3] = point0 - tolerance * direction + tolerance * normal;
    }

    EdgeEnvelope envelope;
    envelope.bounds[0] = corners[0];
    envelope.bounds[1] = corners[0];
    for (std::size_t i = 0; i < corners.size(); ++i) {
        envelope.halfplanes[i][0] = corners[i];
        envelope.halfplanes[i][1] = corners[(i + 1) % corners.size()];
        envelope.bounds[0] = envelope.bounds[0].cwiseMin(corners[i]);
        envelope.bounds[1] = envelope.bounds[1].cwiseMax(corners[i]);
    }
    for (int dimension = 0; dimension < 2; ++dimension) {
        envelope.bounds[0][dimension] =
            std::nextafter(envelope.bounds[0][dimension], -std::numeric_limits<Scalar>::infinity());
        envelope.bounds[1][dimension] =
            std::nextafter(envelope.bounds[1][dimension], std::numeric_limits<Scalar>::infinity());
    }
    return envelope;
}

bool FastEnvelope2D::contains(const EdgeEnvelope& envelope, const Vector2& point)
{
    if (point[0] < envelope.bounds[0][0] || point[0] > envelope.bounds[1][0] ||
        point[1] < envelope.bounds[0][1] || point[1] > envelope.bounds[1][1]) {
        return false;
    }

    for (const HalfPlane& halfplane : envelope.halfplanes) {
        // The rectangle is counter-clockwise. The legacy explicit orientation
        // convention returns negative for a point on the left of an oriented line.
        if (algorithms::orient_2d(halfplane[0], halfplane[1], point) > 0) return false;
    }
    return true;
}

bool FastEnvelope2D::contains(const EdgeEnvelope& envelope, const genericPoint& point)
{
    for (const HalfPlane& halfplane : envelope.halfplanes) {
        const explicitPoint2D line0(halfplane[0][0], halfplane[0][1]);
        const explicitPoint2D line1(halfplane[1][0], halfplane[1][1]);
        // The rectangle is counter-clockwise, so its closed interior is on the
        // left of every edge. Indirect_Predicates returns positive on that side.
        if (genericPoint::orient2D(line0, line1, point) < 0) return false;
    }
    return true;
}

bool FastEnvelope2D::is_outside(const Vector2& point) const
{
    if (envelopes_.empty()) return true;

    std::vector<unsigned int> candidates;
    tree_.point_find_bbox(lift_to_3d(point), candidates);
    return std::none_of(candidates.begin(), candidates.end(), [&](unsigned int candidate) {
        return contains(envelopes_[candidate], point);
    });
}

bool FastEnvelope2D::is_outside(const Vector2& point0, const Vector2& point1) const
{
    if (point0[0] == point1[0] && point0[1] == point1[1]) return is_outside(point0);

    const explicitPoint2D query0(point0[0], point0[1]);
    const explicitPoint2D query1(point1[0], point1[1]);
    std::vector<std::size_t> queue;
    std::vector<bool> reached(envelopes_.size(), false);
    if (envelopes_.empty()) return true;

    std::vector<unsigned int> candidates;
    tree_.segment_find_bbox(point0, point1, candidates);

    // Every box containing the first endpoint is initially reachable. Reaching
    // any box containing the second endpoint proves continuous union coverage.
    for (const unsigned int i : candidates) {
        const bool contains0 = contains(envelopes_[i], point0);
        const bool contains1 = contains(envelopes_[i], point1);
        if (contains0 && contains1) return false;

        if (contains0) {
            reached[i] = true;
            queue.push_back(i);
        }
    }
    if (queue.empty()) return true;

    const auto enqueue_boxes_containing = [&](const auto& point) {
        for (const unsigned int candidate_id : candidates) {
            if (reached[candidate_id]) continue;

            const EdgeEnvelope& candidate = envelopes_[candidate_id];
            if (!contains(candidate, point)) continue;

            if (contains(candidate, point1)) return true;
            reached[candidate_id] = true;
            queue.push_back(candidate_id);
        }
        return false;
    };

    for (std::size_t next = 0; next < queue.size(); ++next) {
        const std::size_t current_id = queue[next];
        const EdgeEnvelope& current = envelopes_[current_id];

        for (const HalfPlane& boundary : current.halfplanes) {
            const int side0 = algorithms::orient_2d(boundary[0], boundary[1], point0);
            const int side1 = algorithms::orient_2d(boundary[0], boundary[1], point1);
            const int boundary_side0 = algorithms::orient_2d(point0, point1, boundary[0]);
            const int boundary_side1 = algorithms::orient_2d(point0, point1, boundary[1]);

            // A boundary vertex on the query handles corner tangencies and
            // collinear overlaps without constructing a degenerate SSI point.
            if (boundary_side0 == 0 && point_in_segment_bounds(boundary[0], point0, point1) &&
                enqueue_boxes_containing(boundary[0]))
                return false;
            if (boundary_side1 == 0 && point_in_segment_bounds(boundary[1], point0, point1) &&
                enqueue_boxes_containing(boundary[1]))
                return false;

            // The remaining finite segment intersection is a proper crossing.
            if (side0 * side1 >= 0 || boundary_side0 * boundary_side1 >= 0) continue;

            const explicitPoint2D line0(boundary[0][0], boundary[0][1]);
            const explicitPoint2D line1(boundary[1][0], boundary[1][1]);

            const implicitPoint2D_SSI clipped_point(query0, query1, line0, line1);
            if (enqueue_boxes_containing(clipped_point)) return false;
        }
    }
    return true;
}

} // namespace fastEnvelope
