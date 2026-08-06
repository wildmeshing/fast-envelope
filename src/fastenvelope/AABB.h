#pragma once

#include "Types.hpp"
#include "common_algorithms.h"

#include <array>
#include <cassert>
#include <vector>

namespace fastEnvelope {
// Query-triangle data that is invariant across every AABB node visited during
// a single triangle_find_bbox traversal. Precomputed once and threaded through
// the recursion so the per-node box-cut test does not recompute the triangle's
// 2D projections or the (box-independent) orient_2d edge signs at every node.
struct TriCutCache
{
    Vector3 tmin, tmax; // triangle AABB
    std::array<std::array<Vector2, 3>, 3> tri2d; // [projection axis][triangle vertex]
    std::array<std::array<int, 3>, 3> ori; // [axis][edge j] triangle-only orient_2d sign
};
class AABB
{
private:
    std::vector<std::array<Vector3, 2>> boxlist;
    size_t n_corners = -1;

    void init_envelope_boxes_recursive(
        const std::vector<std::array<Vector3, 2>>& cornerlist,
        int node_index,
        int b,
        int e);

    void triangle_search_bbd_recursive(
        const TriCutCache& tc,
        std::vector<unsigned int>& list,
        int n,
        int b,
        int e) const;
    void point_search_bbd_recursive(
        const Vector3& point,
        std::vector<unsigned int>& list,
        int n,
        int b,
        int e) const;
    void segment_search_bbd_recursive(
        const Vector3& seg0,
        const Vector3& seg1,
        std::vector<unsigned int>& list,
        int n,
        int b,
        int e) const;


    void bbd_searching_recursive(
        const Vector3& bbd0,
        const Vector3& bbd1,
        std::vector<unsigned int>& list,
        int n,
        int b,
        int e) const;

    static int envelope_max_node_index(int node_index, int b, int e);

    bool is_triangle_cut_bounding_box(const TriCutCache& tc, int index) const;
    bool is_point_cut_bounding_box(const Vector3& p, int index) const;
    bool is_segment_cut_bounding_box(const Vector3& seg0, const Vector3& seg1, int index) const;
    bool is_bbd_cut_bounding_box(const Vector3& bbd0, const Vector3& bbd1, int index) const;

public:
    void init(const std::vector<std::array<Vector3, 2>>& cornerlist);
    bool is_initialized = false;
    inline void triangle_find_bbox(
        const Vector3& triangle0,
        const Vector3& triangle1,
        const Vector3& triangle2,
        std::vector<unsigned int>& list) const
    {
        assert(n_corners >= 0);
        assert(boxlist.size() > 0);
        int de = algorithms::is_triangle_degenerated(triangle0, triangle1, triangle2);
        if (de == DEGENERATED_SEGMENT) {
            Vector3 tmin, tmax;
            algorithms::get_tri_corners(triangle0, triangle1, triangle2, tmin, tmax);
            segment_find_bbox(tmin, tmax, list);
            return;
        } else if (de == DEGENERATED_POINT) {
            point_find_bbox(triangle0, list);
            return;
        }

        // Precompute the query-triangle-invariant data once, then reuse it at
        // every AABB node instead of recomputing it in is_triangle_cut_bounding_box.
        TriCutCache tc;
        algorithms::get_tri_corners(triangle0, triangle1, triangle2, tc.tmin, tc.tmax);
        for (int i = 0; i < 3; i++) {
            tc.tri2d[i][0] = algorithms::to_2d(triangle0, i);
            tc.tri2d[i][1] = algorithms::to_2d(triangle1, i);
            tc.tri2d[i][2] = algorithms::to_2d(triangle2, i);
            for (int j = 0; j < 3; j++) {
                tc.ori[i][j] = algorithms::orient_2d(
                    tc.tri2d[i][(j + 2) % 3],
                    tc.tri2d[i][j % 3],
                    tc.tri2d[i][(j + 1) % 3]);
            }
        }
        triangle_search_bbd_recursive(tc, list, 1, 0, n_corners);
    }

    inline void point_find_bbox(const Vector3& p, std::vector<unsigned int>& list) const
    {
        assert(boxlist.size() > 0);
        point_search_bbd_recursive(p, list, 1, 0, n_corners);
    }
    inline void segment_find_bbox(
        const Vector3& seg0,
        const Vector3& seg1,
        std::vector<unsigned int>& list) const
    {
        assert(boxlist.size() > 0);
        segment_search_bbd_recursive(seg0, seg1, list, 1, 0, n_corners);
    }
    inline void
    bbox_find_bbox(const Vector3& bbd0, const Vector3& bbd1, std::vector<unsigned int>& list) const
    {
        list.clear();
        assert(n_corners >= 0);
        bbd_searching_recursive(bbd0, bbd1, list, 1, 0, n_corners);
    }
};
} // namespace fastEnvelope
