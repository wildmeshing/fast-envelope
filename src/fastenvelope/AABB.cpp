#include "AABB.h"

#include "FastEnvelope.h"

#include <cassert>

namespace fastEnvelope {
void AABB::init_envelope_boxes_recursive(
    const std::vector<std::array<Vector3, 2>>& cornerlist,
    int node_index,
    int b,
    int e)
{
    assert(b != e);
    assert(node_index < boxlist.size());

    if (b + 1 == e) {
        boxlist[node_index] = cornerlist[b];
        return;
    }
    int m = b + (e - b) / 2;
    int childl = 2 * node_index;
    int childr = 2 * node_index + 1;

    assert(childl < boxlist.size());
    assert(childr < boxlist.size());

    init_envelope_boxes_recursive(cornerlist, childl, b, m);
    init_envelope_boxes_recursive(cornerlist, childr, m, e);

    assert(childl < boxlist.size());
    assert(childr < boxlist.size());
    for (int c = 0; c < 3; ++c) {
        boxlist[node_index][0][c] = std::min<Scalar>(boxlist[childl][0][c], boxlist[childr][0][c]);
        boxlist[node_index][1][c] = std::max<Scalar>(boxlist[childl][1][c], boxlist[childr][1][c]);
    }
}

void AABB::triangle_search_bbd_recursive(
    const TriCutCache& tc,
    std::vector<unsigned int>& list,
    int n,
    int b,
    int e) const
{
    assert(e != b);

    assert(n < boxlist.size());
    bool cut = is_triangle_cut_bounding_box(tc, n);

    if (cut == false) return;

    // Leaf case
    if (e == b + 1) {
        list.emplace_back(b);
        return;
    }

    int m = b + (e - b) / 2;
    int childl = 2 * n;
    int childr = 2 * n + 1;

    // assert(childl < boxlist.size());
    // assert(childr < boxlist.size());

    // Traverse the "nearest" child first, so that it has more chances
    // to prune the traversal of the other child.
    triangle_search_bbd_recursive(tc, list, childl, b, m);
    triangle_search_bbd_recursive(tc, list, childr, m, e);
}

void AABB::point_search_bbd_recursive(
    const Vector3& point,
    std::vector<unsigned int>& list,
    int n,
    int b,
    int e) const
{
    assert(e != b);

    assert(n < boxlist.size());
    bool cut = is_point_cut_bounding_box(point, n);

    if (cut == false) return;

    // Leaf case
    if (e == b + 1) {
        list.emplace_back(b);
        return;
    }

    int m = b + (e - b) / 2;
    int childl = 2 * n;
    int childr = 2 * n + 1;

    // assert(childl < boxlist.size());
    // assert(childr < boxlist.size());

    // Traverse the "nearest" child first, so that it has more chances
    // to prune the traversal of the other child.
    point_search_bbd_recursive(point, list, childl, b, m);
    point_search_bbd_recursive(point, list, childr, m, e);
}

void AABB::segment_search_bbd_recursive(
    const Vector3& seg0,
    const Vector3& seg1,
    std::vector<unsigned int>& list,
    int n,
    int b,
    int e) const
{
    assert(e != b);

    assert(n < boxlist.size());
    bool cut = is_segment_cut_bounding_box(seg0, seg1, n);

    if (cut == false) return;

    // Leaf case
    if (e == b + 1) {
        list.emplace_back(b);
        return;
    }

    int m = b + (e - b) / 2;
    int childl = 2 * n;
    int childr = 2 * n + 1;

    // assert(childl < boxlist.size());
    // assert(childr < boxlist.size());

    // Traverse the "nearest" child first, so that it has more chances
    // to prune the traversal of the other child.
    segment_search_bbd_recursive(seg0, seg1, list, childl, b, m);
    segment_search_bbd_recursive(seg0, seg1, list, childr, m, e);
}

void AABB::segment_2d_search_bbd_recursive(
    const Vector2& seg0,
    const Vector2& seg1,
    std::vector<unsigned int>& list,
    int n,
    int b,
    int e) const
{
    assert(e != b);
    assert(n < boxlist.size());
    if (!is_segment_cut_bounding_box(seg0, seg1, n)) return;

    if (e == b + 1) {
        list.emplace_back(b);
        return;
    }

    const int m = b + (e - b) / 2;
    segment_2d_search_bbd_recursive(seg0, seg1, list, 2 * n, b, m);
    segment_2d_search_bbd_recursive(seg0, seg1, list, 2 * n + 1, m, e);
}

void AABB::bbd_searching_recursive(
    const Vector3& bbd0,
    const Vector3& bbd1,
    std::vector<unsigned int>& list,
    int n,
    int b,
    int e) const
{
    assert(e != b);

    assert(n < boxlist.size());
    bool cut = is_bbd_cut_bounding_box(bbd0, bbd1, n);

    if (cut == false) return;

    // Leaf case
    if (e == b + 1) {
        list.emplace_back(b);
        return;
    }

    int m = b + (e - b) / 2;
    int childl = 2 * n;
    int childr = 2 * n + 1;

    // assert(childl < boxlist.size());
    // assert(childr < boxlist.size());

    // Traverse the "nearest" child first, so that it has more chances
    // to prune the traversal of the other child.
    bbd_searching_recursive(bbd0, bbd1, list, childl, b, m);
    bbd_searching_recursive(bbd0, bbd1, list, childr, m, e);
}
int AABB::envelope_max_node_index(int node_index, int b, int e)
{
    assert(e > b);
    if (b + 1 == e) {
        return node_index;
    }
    int m = b + (e - b) / 2;
    int childl = 2 * node_index;
    int childr = 2 * node_index + 1;
    return std::max<int>(
        envelope_max_node_index(childl, b, m),
        envelope_max_node_index(childr, m, e));
}

void AABB::init(const std::vector<std::array<Vector3, 2>>& cornerlist)
{
    n_corners = cornerlist.size();

    boxlist.resize(
        envelope_max_node_index(1, 0, n_corners) +
        1 // <-- this is because size == max_index + 1 !!!
    );

    init_envelope_boxes_recursive(cornerlist, 1, 0, n_corners);
}

bool AABB::is_triangle_cut_bounding_box(const TriCutCache& tc, int index) const
{
    const auto& bmin = boxlist[index][0];
    const auto& bmax = boxlist[index][1];

    // Cheap AABB reject using the precomputed triangle bbox.
    if (!algorithms::box_box_intersection(tc.tmin, tc.tmax, bmin, bmax)) return false;

    std::array<Vector2, 4> mp;
    for (int i = 0; i < 3; i++) {
        // Triangle projection for this axis is precomputed (box-independent).
        const std::array<Vector2, 3>& tri = tc.tri2d[i];

        mp[0] = algorithms::to_2d(bmin, i);
        mp[1] = algorithms::to_2d(bmax, i);
        mp[2][0] = mp[0][0];
        mp[2][1] = mp[1][1];
        mp[3][0] = mp[1][0];
        mp[3][1] = mp[0][1];

        for (int j = 0; j < 3; j++) {
            // The triangle-only edge sign is precomputed (box-independent); when it is
            // zero the original code skipped this edge, so skip before touching the box
            // corners (also saves the four box orient_2d calls in that case).
            const int ori = tc.ori[i][j];
            if (ori == 0) continue;
            const int o0 = algorithms::orient_2d(mp[0], tri[j % 3], tri[(j + 1) % 3]);
            const int o1 = algorithms::orient_2d(mp[1], tri[j % 3], tri[(j + 1) % 3]);
            const int o2 = algorithms::orient_2d(mp[2], tri[j % 3], tri[(j + 1) % 3]);
            const int o3 = algorithms::orient_2d(mp[3], tri[j % 3], tri[(j + 1) % 3]);
            if (ori * o0 <= 0 && ori * o1 <= 0 && ori * o2 <= 0 && ori * o3 <= 0) return false;
        }
    }

    return true;
}
bool AABB::is_point_cut_bounding_box(const Vector3& p, int index) const
{
    const auto& bmin = boxlist[index][0];
    const auto& bmax = boxlist[index][1];
    if (p[0] < bmin[0] || p[1] < bmin[1] || p[2] < bmin[2]) return false;
    if (p[0] > bmax[0] || p[1] > bmax[1] || p[2] > bmax[2]) return false;
    return true;
}

bool AABB::is_segment_cut_bounding_box(const Vector3& seg0, const Vector3& seg1, int index) const
{
    const auto& bmin = boxlist[index][0];
    const auto& bmax = boxlist[index][1];
    Scalar min[3], max[3];
    min[0] = std::min<Scalar>(seg0[0], seg1[0]);
    min[1] = std::min<Scalar>(seg0[1], seg1[1]);
    min[2] = std::min<Scalar>(seg0[2], seg1[2]);
    max[0] = std::max<Scalar>(seg0[0], seg1[0]);
    max[1] = std::max<Scalar>(seg0[1], seg1[1]);
    max[2] = std::max<Scalar>(seg0[2], seg1[2]);
    if (max[0] < bmin[0] || max[1] < bmin[1] || max[2] < bmin[2]) return false;
    if (min[0] > bmax[0] || min[1] > bmax[1] || min[2] > bmax[2]) return false;
    return true;
}

bool AABB::is_segment_cut_bounding_box(const Vector2& seg0, const Vector2& seg1, int index) const
{
    const Vector3& bmin = boxlist[index][0];
    const Vector3& bmax = boxlist[index][1];
    if (std::max(seg0[0], seg1[0]) < bmin[0] || std::min(seg0[0], seg1[0]) > bmax[0] ||
        std::max(seg0[1], seg1[1]) < bmin[1] || std::min(seg0[1], seg1[1]) > bmax[1])
        return false;

    std::array<Vector2, 4> corners;
    corners[0] = Vector2(bmin[0], bmin[1]);
    corners[1] = Vector2(bmax[0], bmin[1]);
    corners[2] = Vector2(bmax[0], bmax[1]);
    corners[3] = Vector2(bmin[0], bmax[1]);

    bool all_positive = true;
    bool all_negative = true;
    for (const Vector2& corner : corners) {
        const int side = algorithms::orient_2d(seg0, seg1, corner);
        all_positive = all_positive && side > 0;
        all_negative = all_negative && side < 0;
    }
    return !all_positive && !all_negative;
}
bool AABB::is_bbd_cut_bounding_box(const Vector3& bbd0, const Vector3& bbd1, int index) const
{
    const auto& bmin = boxlist[index][0];
    const auto& bmax = boxlist[index][1];

    return algorithms::box_box_intersection(bbd0, bbd1, bmin, bmax);
}
} // namespace fastEnvelope
