#include "FastEnvelope.h"

#include "indirectPredicates/ip_filtered.h"

#ifdef ENVELOPE_WITH_GMP
#include "Rational.hpp"
#endif
#include "common_algorithms.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>


#ifdef FAST_ENVELOPE_USE_TIMER
#include <igl/Timer.h>
double time1 = 0, time2 = 0, time3 = 0, time4 = 0;
igl::Timer timer, totaltimer, extimer;
int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
#endif

namespace {
// The records of the facets one query visits (FacetRec: the sub-expressions of the tier-1 TPI and
// LPI filters that depend on the facet alone). A triangle query evaluates each facet it visits
// about 24 times, so its record is computed the first time and reused; it depends on the facet
// only, so this is the same arithmetic on the same doubles, bit for bit. Filled lazily, emptied by
// each query, one per thread.
struct FacetRecordCache
{
    uint64_t envelope = 0; // FastEnvelope::m_instance of the query in progress
    unsigned query = 0;
    std::vector<unsigned> stamp; // per prism: == query when `first` points at its records
    std::vector<int> first; // per prism: index of its first record in `records`
    std::vector<FacetRec> records;
};
thread_local FacetRecordCache t_facet_records;
std::atomic<uint64_t> g_next_instance{1};

[[noreturn]] void facet_record_misuse()
{
    // A record taken outside a query begun on the same envelope, or for a prism that is not one
    // of its candidates: either would read another envelope's records or move the ones already
    // handed out. An internal error, never a property of the input.
    std::fprintf(stderr, "fast-envelope: facet record requested outside its query\n");
    std::abort();
}

// The query triangle of the triangle_out_of_envelope in progress: its record as a TPI plane, and
// the off-plane point n and the three planes through n and an edge that is_3_triangle_cut and the
// cut test in is_triangle_cut_envelope_polyhedra evaluate for every pair of facets.
struct QueryTriangle
{
    bool valid = false; // set for the duration of one triangle_out_of_envelope
    bool n_off_plane = false; // n is off the triangle's plane: the callees need no rand() point
    FacetRec plane;
    const double* points = nullptr;
    fastEnvelope::Vector3 n;
    double edge_points[3][9];
    FacetRec edge_plane[3];
};
thread_local QueryTriangle t_query_triangle;
// The filters read a triangle or a facet as nine contiguous doubles, from the first point's data().
static_assert(
    sizeof(std::array<fastEnvelope::Vector3, 3>) == 9 * sizeof(double),
    "three Vector3 must be nine contiguous doubles");
struct QueryTriangleScope
{
    ~QueryTriangleScope() { t_query_triangle.valid = false; }
};

// is_triangle_cut_envelope_polyhedra's vertex-vs-facet orientations from its last call that went
// through the whole facet loop, and triangle_out_of_envelope's copy of them per filtered prism.
constexpr int kMaxFacets = 16;
thread_local int8_t t_last_orientations[3][kMaxFacets];
thread_local std::vector<std::array<int8_t, 3 * kMaxFacets>> t_prism_orientations;

// The facets of one prism, on the stack: what is_triangle_cut_envelope_polyhedra and
// is_seg_cut_polyhedra used a std::vector for.
struct FacetList
{
    int d[kMaxFacets];
    int n = 0;
    void emplace_back(int v) { d[n++] = v; }
    size_t size() const { return n; }
    int& operator[](size_t i) { return d[i]; }
    int operator[](size_t i) const { return d[i]; }
};
} // namespace

namespace fastEnvelope {
void FastEnvelope::assign_instance()
{
    m_instance = g_next_instance++;
    for (const auto& prism : halfspace) {
        if (prism.size() > kMaxFacets) {
            std::fprintf(stderr, "fast-envelope: a prism with %zu facets\n", prism.size());
            std::abort();
        }
    }
}

void FastEnvelope::begin_query(const std::vector<unsigned int>& prismindex) const
{
    auto& c = t_facet_records;
    c.envelope = m_instance;
    if (c.stamp.size() < halfspace.size()) {
        c.stamp.resize(halfspace.size(), 0);
        c.first.resize(halfspace.size(), 0);
    }
    if (++c.query == 0) {
        std::fill(c.stamp.begin(), c.stamp.end(), 0u);
        c.query = 1;
    }
    // Records are handed out as pointers, so they must not move during the query: room for every
    // facet of every candidate, which are the only prisms a query visits.
    size_t n = 0;
    for (const unsigned int p : prismindex) n += halfspace[p].size();
    c.records.clear();
    c.records.reserve(n);
}

const FacetRec* FastEnvelope::facet_records(int prism) const
{
    auto& c = t_facet_records;
    if (c.envelope != m_instance) facet_record_misuse();
    if (c.stamp[prism] != c.query) {
        if (c.records.size() + halfspace[prism].size() > c.records.capacity()) {
            facet_record_misuse();
        }
        c.stamp[prism] = c.query;
        c.first[prism] = int(c.records.size());
        for (const auto& f : halfspace[prism]) {
            c.records.emplace_back();
            make_facet_rec(f[0].data(), c.records.back());
        }
    }
    return c.records.data() + c.first[prism];
}
} // namespace fastEnvelope

// The facet records of the last two prisms a function asked for, held in locals: most callers
// walk one prism's facets, or a pair, many times over.
#define FACET_RECORD_MEMO                     \
    int memo_prism_a = -1, memo_prism_b = -1; \
    const FacetRec* memo_records_a = nullptr; \
    const FacetRec* memo_records_b = nullptr; \
    (void)memo_prism_a;                       \
    (void)memo_prism_b;                       \
    (void)memo_records_a;                     \
    (void)memo_records_b
#define FACET_RECORDS(P)                       \
    (memo_prism_a == (int)(P) ? memo_records_a \
                              : (memo_records_a = facet_records(memo_prism_a = (int)(P))))
#define FACET_RECORDS_2(P)                     \
    (memo_prism_b == (int)(P) ? memo_records_b \
                              : (memo_records_b = facet_records(memo_prism_b = (int)(P))))
// orient3D_TPI_postfilter / orient3D_LPI_postfilter against facet F of prism P, from its record.
#define TPI_POSTFILTER(st, P, F) \
    orient3D_TPI_postfilter_rec(st, halfspace[P][F][0].data(), FACET_RECORDS(P)[F])
#define LPI_POSTFILTER(sv, X, Y, Z, P, F) \
    orient3D_LPI_postfilter_rec(sv, X, Y, Z, halfspace[P][F][0].data(), FACET_RECORDS(P)[F])

namespace fastEnvelope {

void FastEnvelope::printnumber()
{
#ifdef FAST_ENVELOPE_USE_TIMER
    std::cout << "pre calculation float, " << time1 << std::endl;
    std::cout << "post calculation float, " << time2 << std::endl;
    std::cout << "total float, " << time3 << std::endl;
    std::cout << "exact pre, " << time4 << std::endl;
    std::cout << "exact pre nbr, " << count2 << " out of " << count1 << std::endl;
    std::cout << "exact post nbr, " << count4 << " out of " << count3 << std::endl;
#endif
}

#ifdef ENVELOPE_WITH_GMP
static const std::function<int(fastEnvelope::Rational)> check_rational =
    [](fastEnvelope::Rational v) {
        if (v.get_sign() > 0) return 1;
        if (v.get_sign() < 0) return -1;
        return 0;
    };
#endif

FastEnvelope::FastEnvelope(
    const std::vector<Vector3>& m_ver,
    const std::vector<Vector3i>& m_faces,
    const Scalar eps)
{
    cornerlist.clear();
    halfspace.clear();
    init(m_ver, m_faces, eps);
}

void FastEnvelope::init(
    const std::vector<Vector3>& m_ver,
    const std::vector<Vector3i>& m_faces,
    const Scalar eps)
{
    std::vector<Vector3i> faces_new;


    algorithms::resorting(m_ver, m_faces, faces_new); // resort the facets order

    // algorithms::halfspace_init(m_ver, faces_new, halfspace, cornerlist, eps);
    algorithms::halfspace_generation(m_ver, faces_new, halfspace, cornerlist, eps);

    tree.init(cornerlist);
    assign_instance();

    // initializing types
    initFPU();
}

void FastEnvelope::init(
    const std::vector<Vector3>& m_ver,
    const std::vector<Vector2i>& m_edges,
    const std::vector<Scalar> eps)
{
    assert(eps.size() == m_edges.size());

    std::vector<Vector2i> edges_new;
    std::vector<Scalar> epsnew(m_edges.size());
    std::vector<int> new2old;
    algorithms::resorting(m_ver, m_edges, edges_new, new2old); // resort the facets order
    for (std::size_t i = 0; i < m_edges.size(); ++i) {
        epsnew[i] = eps[new2old[i]];
    }
    // algorithms::halfspace_init(m_ver, edges_new, halfspace, cornerlist, eps);
    algorithms::halfspace_generation(m_ver, edges_new, halfspace, cornerlist, epsnew);

    tree.init(cornerlist);
    assign_instance();

    // initializing types
    initFPU();
}


void FastEnvelope::init(
    const std::vector<Vector3>& m_ver,
    const std::vector<Vector2i>& m_edges,
    const Scalar eps)
{
    std::vector<Vector2i> edges_new;

    algorithms::resorting(m_ver, m_edges, edges_new); // resort the facets order

    // algorithms::halfspace_init(m_ver, edges_new, halfspace, cornerlist, eps);
    algorithms::halfspace_generation(m_ver, edges_new, halfspace, cornerlist, eps);

    tree.init(cornerlist);
    assign_instance();

    // initializing types
    initFPU();
}

void FastEnvelope::init(
    const std::vector<Vector3>& m_ver,
    const std::vector<Vector3i>& m_faces,
    const std::vector<Scalar> eps)
{
    std::vector<Vector3i> faces_new;
    std::vector<Scalar> epsnew;
    epsnew.resize(eps.size());
    std::vector<int> new2old;
    algorithms::resorting(m_ver, m_faces, faces_new, new2old); // resort the facets order
    for (int i = 0; i < eps.size(); i++) {
        epsnew[i] = eps[new2old[i]];
    }
    // algorithms::halfspace_init(m_ver, faces_new, halfspace, cornerlist, eps);
    algorithms::halfspace_generation(m_ver, faces_new, halfspace, cornerlist, epsnew);

    tree.init(cornerlist);
    assign_instance();

    // initializing types
    initFPU();
}

void FastEnvelope::init(
    const std::vector<std::vector<std::array<Vector3, 3>>> halfspace_input,
    const std::vector<std::array<Vector3, 2>> cornerlist_input,
    const Scalar eps)
{
    std::vector<Vector3i> faces_new;

    halfspace = halfspace_input;
    cornerlist = cornerlist_input;
    tree.init(cornerlist);
    assign_instance();
    USE_ADJACENT_INFORMATION = false;
    // initializing types
    initFPU();
}


bool FastEnvelope::is_outside(const std::array<Vector3, 3>& triangle) const
{
    std::vector<unsigned int> querylist;
    tree.triangle_find_bbox(triangle[0], triangle[1], triangle[2], querylist);
#ifdef ENVELOPE_WITH_GMP
    const auto res = triangle_out_simple(triangle, querylist);
    // const auto res = triangle_out_of_envelope(triangle, querylist);
#else
    const auto res = triangle_out_of_envelope(triangle, querylist);
#endif

    return res;
}
bool FastEnvelope::is_outside_no_optimazation(const std::array<Vector3, 3>& triangle) const
{
    std::vector<unsigned int> querylist;
    tree.triangle_find_bbox(triangle[0], triangle[1], triangle[2], querylist);
    const auto res = triangle_out_simple(triangle, querylist);
    return res;
}


bool FastEnvelope::is_outside(const Vector3& point) const
{
    std::vector<unsigned int> querylist;
    tree.point_find_bbox(point, querylist);
    if (querylist.size() == 0) return true;
    bool out = point_out_prism(point, querylist, -1);
    if (out == true) {
        return 1;
    }
    return 0;
}
bool FastEnvelope::is_outside(const Vector3& point0, const Vector3& point1) const
{
    if (point0[0] == point1[0] && point0[1] == point1[1] && point0[2] == point1[2])
        return is_outside(point0);
    std::vector<unsigned int> querylist;
    tree.segment_find_bbox(point0, point1, querylist);
    const bool res = segment_out_of_envelope(point0, point1, querylist);
    return res;
}


bool FastEnvelope::triangle_out_of_envelope(
    const std::array<Vector3, 3>& triangle,
    const std::vector<unsigned int>& prismindex) const
{
    FACET_RECORD_MEMO;
    if (prismindex.size() == 0) {
        return true;
    }
    begin_query(prismindex);
    QueryTriangleScope query_triangle_scope;


    int jump1, jump2;
    static const std::array<std::array<int, 2>, 3> triseg = {{{{0, 1}}, {{0, 2}}, {{1, 2}}}};


    std::vector<unsigned int> filted_intersection;
    filted_intersection.reserve(prismindex.size() / 3);
    std::vector<std::vector<int>> intersect_face;
    intersect_face.reserve(prismindex.size() / 3);
    bool out, cut;

    int inter, inter1, record1, record2,

        tti; // triangle-triangle intersection

    jump1 = -1;

    int check_id;


    for (int i = 0; i < 3; i++) {
        out = point_out_prism_return_local_id(triangle[i], prismindex, jump1, check_id);

        if (out) {
            return true;
        }
    }


    if (prismindex.size() == 1) return false;

    ////////////////////degeneration fix


    int degeneration = algorithms::is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);

    if (degeneration == DEGENERATED_POINT) // case 1 degenerate to a point
    {
        return false;
    }


    if (degeneration == DEGENERATED_SEGMENT) {
        std::vector<unsigned int> queue, idlist;
        queue.emplace_back(check_id); // queue contains the id in prismindex
        idlist.emplace_back(prismindex[check_id]);

        std::vector<int> cidl;
        cidl.reserve(8);
        for (int i = 0; i < queue.size(); i++) {
            jump1 = prismindex[queue[i]];
            int seg_inside = 0;
            for (int k = 0; k < 3; k++) {
                cut = is_seg_cut_polyhedra(
                    jump1,
                    triangle[triseg[k][0]],
                    triangle[triseg[k][1]],
                    cidl);
                if (cut && cidl.size() == 0) {
                    seg_inside++;
                    if (seg_inside == 3)
                        return false; // 3 segs are all totally inside of some polyhedrons
                    continue; // means this seg is inside, check next seg
                }
                if (!cut) continue;

                for (int j = 0; j < cidl.size(); j++) {
#ifdef ENVELOPE_WITH_GMP
                    inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_Rational(
                        triangle[triseg[k][0]],
                        triangle[triseg[k][1]],
                        halfspace[prismindex[queue[i]]][cidl[j]][0],
                        halfspace[prismindex[queue[i]]][cidl[j]][1],
                        halfspace[prismindex[queue[i]]][cidl[j]][2],
                        idlist,
                        jump1,
                        check_id);

#else
                    inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(
                        triangle[triseg[k][0]],
                        triangle[triseg[k][1]],
                        halfspace[prismindex[queue[i]]][cidl[j]][0],
                        halfspace[prismindex[queue[i]]][cidl[j]][1],
                        halfspace[prismindex[queue[i]]][cidl[j]][2],
                        idlist,
                        jump1,
                        check_id);
#endif


                    if (inter == 1) {
#ifdef ENVELOPE_WITH_GMP
                        inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_Rational(
                            triangle[triseg[k][0]],
                            triangle[triseg[k][1]],
                            halfspace[prismindex[queue[i]]][cidl[j]][0],
                            halfspace[prismindex[queue[i]]][cidl[j]][1],
                            halfspace[prismindex[queue[i]]][cidl[j]][2],
                            prismindex,
                            jump1,
                            check_id);
#else
                        inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(
                            triangle[triseg[k][0]],
                            triangle[triseg[k][1]],
                            halfspace[prismindex[queue[i]]][cidl[j]][0],
                            halfspace[prismindex[queue[i]]][cidl[j]][1],
                            halfspace[prismindex[queue[i]]][cidl[j]][2],
                            prismindex,
                            jump1,
                            check_id);
#endif


                        if (inter == 1) {
                            return true;
                        }
                        if (inter == 0) {
                            queue.emplace_back(check_id);
                            idlist.emplace_back(prismindex[check_id]);
                        }
                    }
                }
            }
        }

        return false;
    }
    //
    ////////////////////////////////degeneration fix over


    // The query triangle's plane, its off-plane point and its edge planes, which every
    // facet pair below is tested against (see QueryTriangle).
    {
        QueryTriangle& T = t_query_triangle;
        make_facet_rec(triangle[0].data(), T.plane);
        T.points = triangle[0].data();
        T.n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];
        T.n_off_plane = algorithms::orient_3d(T.n, triangle[0], triangle[1], triangle[2]) != 0;
        if (T.n_off_plane) {
            for (int e = 0; e < 3; e++) {
                const Vector3& ta = triangle[e];
                const Vector3& tb = triangle[(e + 1) % 3];
                double* q = T.edge_points[e];
                q[0] = T.n[0];
                q[1] = T.n[1];
                q[2] = T.n[2];
                q[3] = ta[0];
                q[4] = ta[1];
                q[5] = ta[2];
                q[6] = tb[0];
                q[7] = tb[1];
                q[8] = tb[2];
                make_facet_rec(q, T.edge_plane[e]);
            }
        }
        T.valid = true;
    }
    std::vector<int> cidl;
    cidl.reserve(8);
    t_prism_orientations.clear();
    for (int i = 0; i < prismindex.size(); i++) {
        tti = is_triangle_cut_envelope_polyhedra(
            prismindex[i],
            triangle[0],
            triangle[1],
            triangle[2],
            cidl);
        if (tti == 2) {
            return false; // totally inside of this polyhedron
        } else if (tti == 1 && cidl.size() > 0) {
            filted_intersection.emplace_back(prismindex[i]);
            intersect_face.emplace_back(cidl);
            t_prism_orientations.emplace_back();
            for (int v = 0; v < 3; v++)
                for (int f = 0; f < kMaxFacets; f++)
                    t_prism_orientations.back()[kMaxFacets * v + f] = t_last_orientations[v][f];
        }
    }


    if (filted_intersection.size() == 0) {
        return false;
    }


    std::vector<unsigned int> queue, idlist;
    std::vector<bool> coverlist;
    coverlist.resize(filted_intersection.size());
    for (int i = 0; i < coverlist.size(); i++) {
        coverlist[i] = false; // coverlist shows if the element in filtered_intersection is one of
                              // the current covers
    }
    queue.emplace_back(0); // queue contains the id in filted_intersection
    idlist.emplace_back(
        filted_intersection
            [queue
                 [0]]); // idlist contains the id in prismid//it is fine maybe it is not really intersected
    coverlist[queue[0]] =
        true; // when filted_intersection[i] is already in the cover list, coverlist[i]=true

    std::vector<unsigned int> neighbours; // local id
    std::vector<unsigned int> list;
    std::vector<std::vector<int>> neighbour_facets, idlistorder;
    std::vector<bool> neighbour_cover;
    idlistorder.emplace_back(intersect_face[queue[0]]);


    for (int i = 0; i < queue.size(); i++) {
        jump1 = filted_intersection[queue[i]];

        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < intersect_face[queue[i]].size(); j++) {
                {
                    // seg_cut_plane(t_a, t_b, facet) is FE_CUT_FACE iff orient_3d(t_a, facet) *
                    // orient_3d(t_b, facet) < 0 -- the same predicate on the same points that
                    // is_triangle_cut_envelope_polyhedra already evaluated for this prism.
                    const auto& O = t_prism_orientations[queue[i]];
                    const int f = intersect_face[queue[i]][j];
                    tti = (O[kMaxFacets * triseg[k][0] + f] * O[kMaxFacets * triseg[k][1] + f] < 0)
                              ? FE_CUT_FACE
                              : FE_CUT_COPLANAR;
                }
                if (tti != FE_CUT_FACE) continue;
#ifdef ENVELOPE_WITH_GMP
                inter =
                    Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_Rational(
                        triangle[triseg[k][0]],
                        triangle[triseg[k][1]],
                        halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][0],
                        halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][1],
                        halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][2],
                        idlist,
                        idlistorder,
                        jump1,
                        check_id);
#else
                inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
                    triangle[triseg[k][0]],
                    triangle[triseg[k][1]],
                    halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][0],
                    halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][1],
                    halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]][2],
                    idlist,
                    idlistorder,
                    jump1,
                    check_id);
#endif


                if (inter == 1) {
#ifdef ENVELOPE_WITH_GMP
                    inter =
                        Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over_Rational(
                            triangle[triseg[k][0]],
                            triangle[triseg[k][1]],
                            halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]]
                                     [0],
                            halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]]
                                     [1],
                            halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]]
                                     [2],
                            filted_intersection,
                            intersect_face,
                            coverlist,
                            jump1,
                            check_id);
#else
                    inter =
                        Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
                            triangle[triseg[k][0]],
                            triangle[triseg[k][1]],
                            halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]]
                                     [0],
                            halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]]
                                     [1],
                            halfspace[filted_intersection[queue[i]]][intersect_face[queue[i]][j]]
                                     [2],
                            filted_intersection,
                            intersect_face,
                            coverlist,
                            jump1,
                            check_id);
#endif

                    assert(
                        inter !=
                        2); // the point must exist because it is a seg-halfplane intersection
                    if (inter == 1) {
                        return true;
                    }
                    if (inter == 0) {
                        idlistorder.emplace_back(intersect_face[check_id]);
                        queue.emplace_back(check_id);
                        idlist.emplace_back(filted_intersection[check_id]);
                        coverlist[check_id] = true;
                    }
                }
            }
        }
    }


    // tpi part

    // tree

    const Vector3* hoisted_n = t_query_triangle.n_off_plane ? &t_query_triangle.n : nullptr;
    AABB localtree;
    std::vector<std::array<Vector3, 2>> localcorner;
    localcorner.resize(filted_intersection.size());

    for (int i = 0; i < filted_intersection.size(); i++) {
        localcorner[i] = cornerlist[filted_intersection[i]];
    }

    localtree.init(localcorner);

    // tree end

    for (int i = 1; i < queue.size(); i++) {
        jump1 = filted_intersection[queue[i]];

        localtree.bbox_find_bbox(cornerlist[jump1][0], cornerlist[jump1][1], list);
        neighbours.clear();
        neighbour_cover.clear();
        neighbour_facets.clear();

        neighbours.resize(list.size());
        neighbour_facets.resize(list.size());
        neighbour_cover.resize(list.size());
        for (int j = 0; j < list.size(); j++) {
            neighbours[j] = filted_intersection[list[j]];
            neighbour_facets[j] = intersect_face[list[j]];
            if (coverlist[list[j]] == true)
                neighbour_cover[j] = true;
            else
                neighbour_cover[j] = false;
        }


        for (int j = 0; j < i; j++) {
            jump2 = filted_intersection[queue[j]];
            if (!algorithms::box_box_intersection(
                    cornerlist[jump1][0],
                    cornerlist[jump1][1],
                    cornerlist[jump2][0],
                    cornerlist[jump2][1]))
                continue;
            for (int k = 0; k < intersect_face[queue[i]].size(); k++) {
                for (int h = 0; h < intersect_face[queue[j]].size(); h++) {
#ifdef ENVELOPE_WITH_GMP
                    cut = is_3_triangle_cut_Rational(
                        triangle,
                        halfspace[jump1][intersect_face[queue[i]][k]][0],
                        halfspace[jump1][intersect_face[queue[i]][k]][1],
                        halfspace[jump1][intersect_face[queue[i]][k]][2],

                        halfspace[jump2][intersect_face[queue[j]][h]][0],
                        halfspace[jump2][intersect_face[queue[j]][h]][1],
                        halfspace[jump2][intersect_face[queue[j]][h]][2]);

                    // timettt += timer.getElapsedTimeInSec();
                    if (!cut) continue;


                    cut = is_tpp_on_polyhedra_Rational(
                        triangle,
                        halfspace[jump1][intersect_face[queue[i]][k]][0],
                        halfspace[jump1][intersect_face[queue[i]][k]][1],
                        halfspace[jump1][intersect_face[queue[i]][k]][2],

                        halfspace[jump2][intersect_face[queue[j]][h]][0],
                        halfspace[jump2][intersect_face[queue[j]][h]][1],
                        halfspace[jump2][intersect_face[queue[j]][h]][2],
                        jump1,
                        intersect_face[queue[i]][k]);

                    if (!cut) continue;


                    cut = is_tpp_on_polyhedra_Rational(
                        triangle,
                        halfspace[jump1][intersect_face[queue[i]][k]][0],
                        halfspace[jump1][intersect_face[queue[i]][k]][1],
                        halfspace[jump1][intersect_face[queue[i]][k]][2],

                        halfspace[jump2][intersect_face[queue[j]][h]][0],
                        halfspace[jump2][intersect_face[queue[j]][h]][1],
                        halfspace[jump2][intersect_face[queue[j]][h]][2],
                        jump2,
                        intersect_face[queue[j]][h]);

                    if (!cut) continue;
#else
#ifdef FAST_ENVELOPE_USE_TIMER
                    totaltimer.start();
#endif
                    // One prefilter per (triangle, facet k, facet h) point, shared by the five
                    // functions below that used to recompute it from the same 27 doubles.
                    TPI_filtered_suppvars st_shared;
                    const int pre_shared =
                        orient3D_TPI_prefilter_rec(
                            t_query_triangle.plane,
                            FACET_RECORDS(jump1)[intersect_face[queue[i]][k]],
                            FACET_RECORDS_2(jump2)[intersect_face[queue[j]][h]],
                            t_query_triangle.points,
                            halfspace[jump1][intersect_face[queue[i]][k]][0].data(),
                            halfspace[jump2][intersect_face[queue[j]][h]][0].data(),
                            st_shared)
                            ? 1
                            : 0;
                    const TPI_filtered_suppvars* pst_shared = &st_shared;
                    cut = is_3_triangle_cut(
                        triangle,
                        halfspace[jump1][intersect_face[queue[i]][k]][0],
                        halfspace[jump1][intersect_face[queue[i]][k]][1],
                        halfspace[jump1][intersect_face[queue[i]][k]][2],

                        halfspace[jump2][intersect_face[queue[j]][h]][0],
                        halfspace[jump2][intersect_face[queue[j]][h]][1],
                        halfspace[jump2][intersect_face[queue[j]][h]][2],
                        pst_shared,
                        pre_shared,
                        hoisted_n);
#ifdef FAST_ENVELOPE_USE_TIMER
                    time3 += totaltimer.getElapsedTimeInSec();
#endif
                    // timettt += timer.getElapsedTimeInSec();
                    if (!cut) continue;


                    cut = is_tpp_on_polyhedra(
                        triangle,
                        halfspace[jump1][intersect_face[queue[i]][k]][0],
                        halfspace[jump1][intersect_face[queue[i]][k]][1],
                        halfspace[jump1][intersect_face[queue[i]][k]][2],

                        halfspace[jump2][intersect_face[queue[j]][h]][0],
                        halfspace[jump2][intersect_face[queue[j]][h]][1],
                        halfspace[jump2][intersect_face[queue[j]][h]][2],
                        jump1,
                        intersect_face[queue[i]][k],
                        pst_shared,
                        pre_shared);

                    if (!cut) continue;


                    cut = is_tpp_on_polyhedra(
                        triangle,
                        halfspace[jump1][intersect_face[queue[i]][k]][0],
                        halfspace[jump1][intersect_face[queue[i]][k]][1],
                        halfspace[jump1][intersect_face[queue[i]][k]][2],

                        halfspace[jump2][intersect_face[queue[j]][h]][0],
                        halfspace[jump2][intersect_face[queue[j]][h]][1],
                        halfspace[jump2][intersect_face[queue[j]][h]][2],
                        jump2,
                        intersect_face[queue[j]][h],
                        pst_shared,
                        pre_shared);

                    if (!cut) continue;
#endif


#ifdef ENVELOPE_WITH_GMP
                    inter =
                        Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_Rational(
                            triangle,
                            halfspace[jump1][intersect_face[queue[i]][k]][0],
                            halfspace[jump1][intersect_face[queue[i]][k]][1],
                            halfspace[jump1][intersect_face[queue[i]][k]][2],

                            halfspace[jump2][intersect_face[queue[j]][h]][0],
                            halfspace[jump2][intersect_face[queue[j]][h]][1],
                            halfspace[jump2][intersect_face[queue[j]][h]][2],
                            idlist,
                            idlistorder,
                            jump1,
                            jump2,
                            check_id);
#else
                    inter =
                        Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
                            triangle,
                            halfspace[jump1][intersect_face[queue[i]][k]][0],
                            halfspace[jump1][intersect_face[queue[i]][k]][1],
                            halfspace[jump1][intersect_face[queue[i]][k]][2],

                            halfspace[jump2][intersect_face[queue[j]][h]][0],
                            halfspace[jump2][intersect_face[queue[j]][h]][1],
                            halfspace[jump2][intersect_face[queue[j]][h]][2],
                            idlist,
                            idlistorder,
                            jump1,
                            jump2,
                            check_id,
                            pst_shared,
                            pre_shared);
#endif


                    if (inter == 1) {
#ifdef ENVELOPE_WITH_GMP
                        inter =
                            Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_Rational(
                                triangle,
                                halfspace[jump1][intersect_face[queue[i]][k]][0],
                                halfspace[jump1][intersect_face[queue[i]][k]][1],
                                halfspace[jump1][intersect_face[queue[i]][k]][2],

                                halfspace[jump2][intersect_face[queue[j]][h]][0],
                                halfspace[jump2][intersect_face[queue[j]][h]][1],
                                halfspace[jump2][intersect_face[queue[j]][h]][2],
                                neighbours,
                                neighbour_facets,
                                jump1,
                                jump2,
                                check_id);
#else
                        inter =
                            Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
                                triangle,
                                halfspace[jump1][intersect_face[queue[i]][k]][0],
                                halfspace[jump1][intersect_face[queue[i]][k]][1],
                                halfspace[jump1][intersect_face[queue[i]][k]][2],

                                halfspace[jump2][intersect_face[queue[j]][h]][0],
                                halfspace[jump2][intersect_face[queue[j]][h]][1],
                                halfspace[jump2][intersect_face[queue[j]][h]][2],
                                neighbours,
                                neighbour_facets,
                                neighbour_cover,
                                jump1,
                                jump2,
                                check_id,
                                pst_shared,
                                pre_shared);
#endif


                        if (inter == 1) {
                            return true;
                        }
                        if (inter == 0) {
                            idlistorder.emplace_back(intersect_face[list[check_id]]);
                            queue.emplace_back(list[check_id]);
                            idlist.emplace_back(filted_intersection[list[check_id]]);
                            coverlist[list[check_id]] = true;
                        }
                    }
                }
            }
        }
    }


    return false;
}

bool FastEnvelope::triangle_out_simple(
    const std::array<Vector3, 3>& triangle,
    const std::vector<unsigned int>& prismindex) const
{
    if (prismindex.size() == 0) {
        return true;
    }
    begin_query(prismindex);


    int jump1, jump2;
    static const std::array<std::array<int, 2>, 3> triseg = {{{{0, 1}}, {{0, 2}}, {{1, 2}}}};
    bool out, cut;
    int inter, check_id;

    jump1 = -1;
    for (int i = 0; i < 3; i++) {
        out = point_out_prism(triangle[i], prismindex, jump1);

        if (out) {
            return true;
        }
    }

    if (prismindex.size() == 1) return false;

    ////////////////////degeneration fix


    int degeneration = algorithms::is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);

    if (degeneration == DEGENERATED_POINT) // case 1 degenerate to a point
    {
        return false;
    }


    if (degeneration == DEGENERATED_SEGMENT) {
        for (int i = 0; i < prismindex.size(); i++) { // for each prisms
            jump1 = prismindex[i];
            for (int k = 0; k < 3; k++) { // for each edge
                for (int j = 0; j < halfspace[jump1].size(); j++) { // for each plane
                    int tti = algorithms::seg_cut_plane(
                        triangle[triseg[k][0]],
                        triangle[triseg[k][1]],
                        halfspace[jump1][j][0],
                        halfspace[jump1][j][1],
                        halfspace[jump1][j][2]); // make sure the intersection is on the edge
                    if (tti != FE_CUT_FACE) continue;
#ifdef ENVELOPE_WITH_GMP
                    inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_Rational(
                        triangle[triseg[k][0]],
                        triangle[triseg[k][1]],
                        halfspace[jump1][j][0],
                        halfspace[jump1][j][1],
                        halfspace[jump1][j][2],
                        prismindex,
                        jump1,
                        check_id);
#else
                    inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(
                        triangle[triseg[k][0]],
                        triangle[triseg[k][1]],
                        halfspace[jump1][j][0],
                        halfspace[jump1][j][1],
                        halfspace[jump1][j][2],
                        prismindex,
                        jump1,
                        check_id);
#endif
                    if (inter == 1) {
                        return true;
                    }
                }
            }
        }

        return false;
    }
    //
    ////////////////////////////////degeneration fix over


    // lpi

    for (int i = 0; i < prismindex.size(); i++) { // for each prisms
        jump1 = prismindex[i];
        for (int k = 0; k < 3; k++) { // for each edge
            for (int j = 0; j < halfspace[jump1].size(); j++) { // for each plane
                int tti = algorithms::seg_cut_plane(
                    triangle[triseg[k][0]],
                    triangle[triseg[k][1]],
                    halfspace[jump1][j][0],
                    halfspace[jump1][j][1],
                    halfspace[jump1][j][2]); // make sure the intersection is on the edge
                if (tti != FE_CUT_FACE) continue;

#ifdef ENVELOPE_WITH_GMP
                inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_Rational(
                    triangle[triseg[k][0]],
                    triangle[triseg[k][1]],
                    halfspace[jump1][j][0],
                    halfspace[jump1][j][1],
                    halfspace[jump1][j][2],
                    prismindex,
                    jump1,
                    check_id);
#else
                inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(
                    triangle[triseg[k][0]],
                    triangle[triseg[k][1]],
                    halfspace[jump1][j][0],
                    halfspace[jump1][j][1],
                    halfspace[jump1][j][2],
                    prismindex,
                    jump1,
                    check_id);
#endif
                if (inter == 1) {
                    return true;
                }
            }
        }
    }

    for (int i = 0; i < prismindex.size(); i++) { // for each prisms
        jump1 = prismindex[i];

        for (int j = 0; j < i; j++) {
            jump2 = prismindex[j];

            if (!algorithms::box_box_intersection(
                    cornerlist[jump1][0],
                    cornerlist[jump1][1],
                    cornerlist[jump2][0],
                    cornerlist[jump2][1]))
                continue;


            for (int k = 0; k < halfspace[jump1].size(); k++) {
                for (int h = 0; h < halfspace[jump2].size(); h++) {
#ifdef ENVELOPE_WITH_GMP
                    cut = is_3_triangle_cut_Rational(
                        triangle,
                        halfspace[jump1][k][0],
                        halfspace[jump1][k][1],
                        halfspace[jump1][k][2],

                        halfspace[jump2][h][0],
                        halfspace[jump2][h][1],
                        halfspace[jump2][h][2]);
#else
                    cut = is_3_triangle_cut(
                        triangle,
                        halfspace[jump1][k][0],
                        halfspace[jump1][k][1],
                        halfspace[jump1][k][2],

                        halfspace[jump2][h][0],
                        halfspace[jump2][h][1],
                        halfspace[jump2][h][2]);

#endif
                    // timettt += timer.getElapsedTimeInSec();
                    if (!cut) continue;

#ifdef ENVELOPE_WITH_GMP
                    inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_Rational(
                        triangle, // can not have any bounding box checking inside
                        halfspace[jump1][k][0],
                        halfspace[jump1][k][1],
                        halfspace[jump1][k][2],

                        halfspace[jump2][h][0],
                        halfspace[jump2][h][1],
                        halfspace[jump2][h][2],
                        prismindex,
                        jump1,
                        jump2);
#else
                    inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism(
                        triangle, // can not have any bounding box checking inside
                        halfspace[jump1][k][0],
                        halfspace[jump1][k][1],
                        halfspace[jump1][k][2],

                        halfspace[jump2][h][0],
                        halfspace[jump2][h][1],
                        halfspace[jump2][h][2],
                        prismindex,
                        jump1,
                        jump2);
#endif
                    if (inter == 1) {
                        return true;
                    }
                }
            }
        }
    }

    return false;
}


bool FastEnvelope::segment_out_of_envelope(
    const Vector3& seg0,
    const Vector3& seg1,
    const std::vector<unsigned int>& prismindex) const
{
    if (prismindex.size() == 0) {
        return true;
    }
    begin_query(prismindex);


    int jump1;


    bool out, cut;

    int inter,

        tti; // triangle-triangle intersection

    jump1 = -1;

    int check_id, check_id1;

    out = point_out_prism_return_local_id(seg0, prismindex, jump1, check_id);

    if (out) {
        return true;
    }
    out = point_out_prism_return_local_id(seg1, prismindex, jump1, check_id1);

    if (out) {
        return true;
    }
    if (check_id == check_id1) return false;

    if (prismindex.size() == 1) return false;


    std::vector<unsigned int> queue, idlist;
    queue.emplace_back(check_id); // queue contains the id in prismindex
    idlist.emplace_back(prismindex[check_id]);

    std::vector<int> cidl;
    cidl.reserve(8);
    for (int i = 0; i < queue.size(); i++) {
        jump1 = prismindex[queue[i]];

        cut = is_seg_cut_polyhedra(jump1, seg0, seg1, cidl);
        if (cut && cidl.size() == 0) return false;
        if (!cut) continue;
        for (int j = 0; j < cidl.size(); j++) {
            inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(
                seg0,
                seg1,
                halfspace[prismindex[queue[i]]][cidl[j]][0],
                halfspace[prismindex[queue[i]]][cidl[j]][1],
                halfspace[prismindex[queue[i]]][cidl[j]][2],
                idlist,
                jump1,
                check_id);

            if (inter == 1) {
                inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(
                    seg0,
                    seg1,
                    halfspace[prismindex[queue[i]]][cidl[j]][0],
                    halfspace[prismindex[queue[i]]][cidl[j]][1],
                    halfspace[prismindex[queue[i]]][cidl[j]][2],
                    prismindex,
                    jump1,
                    check_id);

                if (inter == 1) {
                    return true;
                }
                if (inter == 0) {
                    queue.emplace_back(check_id);
                    idlist.emplace_back(prismindex[check_id]);
                }
            }
        }
    }

    return false;
}

bool FastEnvelope::debugcode(
    const std::array<Vector3, 3>& triangle,
    const std::vector<unsigned int>& prismindex) const
{
    return false;
}


int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(
    const Vector3& segpoint0,
    const Vector3& segpoint1,
    const Vector3& triangle0,
    const Vector3& triangle1,
    const Vector3& triangle2,
    const std::vector<unsigned int>& prismindex,
    const int& jump,
    int& id) const
{
    FACET_RECORD_MEMO;
    LPI_filtered_suppvars sf;
    bool precom = orient3D_LPI_prefilter( //
            segpoint0[0], segpoint0[1], segpoint0[2],
            segpoint1[0], segpoint1[1], segpoint1[2],
            triangle0[0], triangle0[1], triangle0[2],
            triangle1[0], triangle1[1], triangle1[2],
            triangle2[0], triangle2[1], triangle2[2],
            sf);
    if (precom == false) {
        int ori;


        LPI_exact_suppvars s;
        bool premulti = orient3D_LPI_pre_exact(
            segpoint0[0],
            segpoint0[1],
            segpoint0[2],
            segpoint1[0],
            segpoint1[1],
            segpoint1[2],

            triangle0[0],
            triangle0[1],
            triangle0[2],
            triangle1[0],
            triangle1[1],
            triangle1[2],
            triangle2[0],
            triangle2[1],
            triangle2[2],
            s);


        if (premulti == false) return 2;
        for (int i = 0; i < prismindex.size(); i++) {
            if (prismindex[i] == jump) {
                continue;
            }


            for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
                ori = orient3D_LPI_post_exact(
                    s,
                    segpoint0[0],
                    segpoint0[1],
                    segpoint0[2],

                    halfspace[prismindex[i]][j][0][0],
                    halfspace[prismindex[i]][j][0][1],
                    halfspace[prismindex[i]][j][0][2],

                    halfspace[prismindex[i]][j][1][0],
                    halfspace[prismindex[i]][j][1][1],
                    halfspace[prismindex[i]][j][1][2],

                    halfspace[prismindex[i]][j][2][0],
                    halfspace[prismindex[i]][j][2][1],
                    halfspace[prismindex[i]][j][2][2]);
                if (ori == 1 || ori == 0) {
                    break;
                }
            }
            if (ori == -1) {
                id = i;
                return IN_PRISM;
            }
        }


        return OUT_PRISM;
    }
    int tot;
    int ori;
    INDEX index;
    std::vector<INDEX> recompute;

    recompute.clear();
    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump) continue;

        index.FACES.clear();
        tot = 0;


        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            ori = LPI_POSTFILTER(sf, segpoint0[0], segpoint0[1], segpoint0[2], prismindex[i], j);

            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                index.FACES.emplace_back(j);
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (tot == halfspace[prismindex[i]].size()) {
            id = i;
            return IN_PRISM;
        }

        if (ori != 1) {
            assert(!index.FACES.empty());
            index.Pi = i; // local id
            recompute.emplace_back(index);
        }
    }

    if (!recompute.empty()) {
        LPI_exact_suppvars s;
        bool premulti = orient3D_LPI_pre_exact(
            segpoint0[0],
            segpoint0[1],
            segpoint0[2],
            segpoint1[0],
            segpoint1[1],
            segpoint1[2],
            triangle0[0],
            triangle0[1],
            triangle0[2],
            triangle1[0],
            triangle1[1],
            triangle1[2],
            triangle2[0],
            triangle2[1],
            triangle2[2],
            s);


        for (int k = 0; k < recompute.size(); k++) {
            int in1 = prismindex[recompute[k].Pi];

            for (int j = 0; j < recompute[k].FACES.size(); j++) {
                int in2 = recompute[k].FACES[j];

                ori = orient3D_LPI_post_exact(
                    s,
                    segpoint0[0],
                    segpoint0[1],
                    segpoint0[2],
                    halfspace[in1][in2][0][0],
                    halfspace[in1][in2][0][1],
                    halfspace[in1][in2][0][2],
                    halfspace[in1][in2][1][0],
                    halfspace[in1][in2][1][1],
                    halfspace[in1][in2][1][2],
                    halfspace[in1][in2][2][0],
                    halfspace[in1][in2][2][1],
                    halfspace[in1][in2][2][2]);

                if (ori == 1 || ori == 0) break;
            }

            if (ori == -1) {
                id = recompute[k].Pi;
                return IN_PRISM;
            }
        }
    }

    return OUT_PRISM;
}


#ifdef ENVELOPE_WITH_GMP
int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_Rational(
    const Vector3& segpoint0,
    const Vector3& segpoint1,
    const Vector3& triangle0,
    const Vector3& triangle1,
    const Vector3& triangle2,
    const std::vector<unsigned int>& prismindex,
    const int& jump,
    int& id) const
{
    Rational s00(segpoint0[0]), s01(segpoint0[1]), s02(segpoint0[2]), s10(segpoint1[0]),
        s11(segpoint1[1]), s12(segpoint1[2]), t00(triangle0[0]), t01(triangle0[1]),
        t02(triangle0[2]), t10(triangle1[0]), t11(triangle1[1]), t12(triangle1[2]),
        t20(triangle2[0]), t21(triangle2[1]), t22(triangle2[2]), a11, a12, a13, d;
    bool precom = orient3D_LPI_prefilter_multiprecision(
        s00,
        s01,
        s02,
        s10,
        s11,
        s12,
        t00,
        t01,
        t02,
        t10,
        t11,
        t12,
        t20,
        t21,
        t22,
        a11,
        a12,
        a13,
        d,
        check_rational);

    if (precom == false) {
        return 2;
    }
    int tot;
    int ori;

    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump) continue;


        tot = 0;

        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            Rational h00(halfspace[prismindex[i]][j][0][0]), h01(halfspace[prismindex[i]][j][0][1]),
                h02(halfspace[prismindex[i]][j][0][2]), h10(halfspace[prismindex[i]][j][1][0]),
                h11(halfspace[prismindex[i]][j][1][1]), h12(halfspace[prismindex[i]][j][1][2]),
                h20(halfspace[prismindex[i]][j][2][0]), h21(halfspace[prismindex[i]][j][2][1]),
                h22(halfspace[prismindex[i]][j][2][2]);
            ori = orient3D_LPI_postfilter_multiprecision(
                a11,
                a12,
                a13,
                d,
                s00,
                s01,
                s02,
                h00,
                h01,
                h02,
                h10,
                h11,
                h12,
                h20,
                h21,
                h22,
                check_rational);

            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                break;
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (tot == halfspace[prismindex[i]].size()) {
            id = i;
            return IN_PRISM;
        }
    }
    return OUT_PRISM;
}
#endif

int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
    const Vector3& segpoint0,
    const Vector3& segpoint1,
    const Vector3& triangle0,
    const Vector3& triangle1,
    const Vector3& triangle2,
    const std::vector<unsigned int>& prismindex,
    const std::vector<std::vector<int>>& intersect_face,
    const int& jump,
    int& id) const
{
    FACET_RECORD_MEMO;
    LPI_filtered_suppvars sl;
    bool precom = orient3D_LPI_prefilter( //
            segpoint0[0], segpoint0[1], segpoint0[2],
            segpoint1[0], segpoint1[1], segpoint1[2],
            triangle0[0], triangle0[1], triangle0[2],
            triangle1[0], triangle1[1], triangle1[2],
            triangle2[0], triangle2[1], triangle2[2],
            sl);
    if (precom == false) {
        int tot, ori, fid;

        LPI_exact_suppvars s;
        bool premulti = orient3D_LPI_pre_exact(
            segpoint0[0],
            segpoint0[1],
            segpoint0[2],
            segpoint1[0],
            segpoint1[1],
            segpoint1[2],

            triangle0[0],
            triangle0[1],
            triangle0[2],
            triangle1[0],
            triangle1[1],
            triangle1[2],
            triangle2[0],
            triangle2[1],
            triangle2[2],
            s);


        if (premulti == false) return 2;
        for (int i = 0; i < prismindex.size(); i++) {
            if (prismindex[i] == jump) {
                continue;
            }
            tot = 0;
            fid = 0;
            ori = -1;
            for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
                if (intersect_face[i][fid] == j) {
                    if (fid + 1 < intersect_face[i].size()) fid++;
                } else
                    continue;
                ori = orient3D_LPI_post_exact(
                    s,
                    segpoint0[0],
                    segpoint0[1],
                    segpoint0[2],

                    halfspace[prismindex[i]][j][0][0],
                    halfspace[prismindex[i]][j][0][1],
                    halfspace[prismindex[i]][j][0][2],

                    halfspace[prismindex[i]][j][1][0],
                    halfspace[prismindex[i]][j][1][1],
                    halfspace[prismindex[i]][j][1][2],

                    halfspace[prismindex[i]][j][2][0],
                    halfspace[prismindex[i]][j][2][1],
                    halfspace[prismindex[i]][j][2][2]);
                if (ori == 1 || ori == 0) {
                    break;
                }

                if (ori == -1) {
                    tot++;
                }
            }
            if (ori == 1 || ori == 0) continue;
            fid = 0;
            ori = -1;
            for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
                if (intersect_face[i][fid] == j) {
                    if (fid + 1 < intersect_face[i].size()) fid++;
                    continue;
                }

                ori = orient3D_LPI_post_exact(
                    s,
                    segpoint0[0],
                    segpoint0[1],
                    segpoint0[2],

                    halfspace[prismindex[i]][j][0][0],
                    halfspace[prismindex[i]][j][0][1],
                    halfspace[prismindex[i]][j][0][2],

                    halfspace[prismindex[i]][j][1][0],
                    halfspace[prismindex[i]][j][1][1],
                    halfspace[prismindex[i]][j][1][2],

                    halfspace[prismindex[i]][j][2][0],
                    halfspace[prismindex[i]][j][2][1],
                    halfspace[prismindex[i]][j][2][2]);
                if (ori == 1 || ori == 0) {
                    break;
                }

                if (ori == -1) {
                    tot++;
                }
            }
            if (ori == 1 || ori == 0) continue;
            if (tot == halfspace[prismindex[i]].size()) {
                id = i;
                return IN_PRISM;
            }
        }


        return OUT_PRISM;
    }
    int tot, fid;
    int ori;
    INDEX index;
    std::vector<INDEX> recompute;

    recompute.clear();
    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump) continue;

        index.FACES.clear();
        tot = 0;
        fid = 0;
        ori = -1;

        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
            } else
                continue;
            ori = LPI_POSTFILTER(sl, segpoint0[0], segpoint0[1], segpoint0[2], prismindex[i], j);

            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                index.FACES.emplace_back(j);
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1) continue;
        fid = 0;
        ori = -1;
        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
                continue;
            }

            ori = LPI_POSTFILTER(sl, segpoint0[0], segpoint0[1], segpoint0[2], prismindex[i], j);

            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                index.FACES.emplace_back(j);
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1) continue;
        if (tot == halfspace[prismindex[i]].size()) {
            id = i;
            return IN_PRISM;
        }

        if (ori != 1) {
            assert(!index.FACES.empty());
            index.Pi = i; // local id
            recompute.emplace_back(index);
        }
    }

    if (!recompute.empty()) {
        LPI_exact_suppvars s;
        bool premulti = orient3D_LPI_pre_exact(
            segpoint0[0],
            segpoint0[1],
            segpoint0[2],
            segpoint1[0],
            segpoint1[1],
            segpoint1[2],
            triangle0[0],
            triangle0[1],
            triangle0[2],
            triangle1[0],
            triangle1[1],
            triangle1[2],
            triangle2[0],
            triangle2[1],
            triangle2[2],
            s);


        for (int k = 0; k < recompute.size(); k++) {
            int in1 = prismindex[recompute[k].Pi];

            for (int j = 0; j < recompute[k].FACES.size(); j++) {
                int in2 = recompute[k].FACES[j];

                ori = orient3D_LPI_post_exact(
                    s,
                    segpoint0[0],
                    segpoint0[1],
                    segpoint0[2],
                    halfspace[in1][in2][0][0],
                    halfspace[in1][in2][0][1],
                    halfspace[in1][in2][0][2],
                    halfspace[in1][in2][1][0],
                    halfspace[in1][in2][1][1],
                    halfspace[in1][in2][1][2],
                    halfspace[in1][in2][2][0],
                    halfspace[in1][in2][2][1],
                    halfspace[in1][in2][2][2]);

                if (ori == 1 || ori == 0) break;
            }

            if (ori == -1) {
                id = recompute[k].Pi;
                return IN_PRISM;
            }
        }
    }

    return OUT_PRISM;
}

#ifdef ENVELOPE_WITH_GMP
int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_Rational(
    const Vector3& segpoint0,
    const Vector3& segpoint1,
    const Vector3& triangle0,
    const Vector3& triangle1,
    const Vector3& triangle2,
    const std::vector<unsigned int>& prismindex,
    const std::vector<std::vector<int>>& intersect_face,
    const int& jump,
    int& id) const
{
    Rational s00(segpoint0[0]), s01(segpoint0[1]), s02(segpoint0[2]), s10(segpoint1[0]),
        s11(segpoint1[1]), s12(segpoint1[2]), t00(triangle0[0]), t01(triangle0[1]),
        t02(triangle0[2]), t10(triangle1[0]), t11(triangle1[1]), t12(triangle1[2]),
        t20(triangle2[0]), t21(triangle2[1]), t22(triangle2[2]), a11, a12, a13, d;
    bool precom = orient3D_LPI_prefilter_multiprecision(
        s00,
        s01,
        s02,
        s10,
        s11,
        s12,
        t00,
        t01,
        t02,
        t10,
        t11,
        t12,
        t20,
        t21,
        t22,
        a11,
        a12,
        a13,
        d,
        check_rational);

    if (precom == false) {
        return 2;
    }

    int tot, fid;
    int ori;

    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump) continue;

        tot = 0;
        fid = 0;
        ori = -1;

        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
            } else
                continue;

            Rational h00(halfspace[prismindex[i]][j][0][0]), h01(halfspace[prismindex[i]][j][0][1]),
                h02(halfspace[prismindex[i]][j][0][2]), h10(halfspace[prismindex[i]][j][1][0]),
                h11(halfspace[prismindex[i]][j][1][1]), h12(halfspace[prismindex[i]][j][1][2]),
                h20(halfspace[prismindex[i]][j][2][0]), h21(halfspace[prismindex[i]][j][2][1]),
                h22(halfspace[prismindex[i]][j][2][2]);
            ori = orient3D_LPI_postfilter_multiprecision(
                a11,
                a12,
                a13,
                d,
                s00,
                s01,
                s02,
                h00,
                h01,
                h02,
                h10,
                h11,
                h12,
                h20,
                h21,
                h22,
                check_rational);


            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                break;
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1 || ori == 0) continue;
        fid = 0;
        ori = -1;
        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
                continue;
            }
            Rational h00(halfspace[prismindex[i]][j][0][0]), h01(halfspace[prismindex[i]][j][0][1]),
                h02(halfspace[prismindex[i]][j][0][2]), h10(halfspace[prismindex[i]][j][1][0]),
                h11(halfspace[prismindex[i]][j][1][1]), h12(halfspace[prismindex[i]][j][1][2]),
                h20(halfspace[prismindex[i]][j][2][0]), h21(halfspace[prismindex[i]][j][2][1]),
                h22(halfspace[prismindex[i]][j][2][2]);
            ori = orient3D_LPI_postfilter_multiprecision(
                a11,
                a12,
                a13,
                d,
                s00,
                s01,
                s02,
                h00,
                h01,
                h02,
                h10,
                h11,
                h12,
                h20,
                h21,
                h22,
                check_rational);

            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                break;
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1 || ori == 0) continue;
        if (tot == halfspace[prismindex[i]].size()) {
            id = i;
            return IN_PRISM;
        }
    }

    return OUT_PRISM;
}
#endif
int FastEnvelope::Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
    const Vector3& segpoint0,
    const Vector3& segpoint1,
    const Vector3& triangle0,
    const Vector3& triangle1,
    const Vector3& triangle2,
    const std::vector<unsigned int>& prismindex,
    const std::vector<std::vector<int>>& intersect_face,
    const std::vector<bool>& coverlist,
    const int& jump,
    int& id) const
{
    FACET_RECORD_MEMO;
    LPI_filtered_suppvars sl;
    bool precom = orient3D_LPI_prefilter( //
            segpoint0[0], segpoint0[1], segpoint0[2],
            segpoint1[0], segpoint1[1], segpoint1[2],
            triangle0[0], triangle0[1], triangle0[2],
            triangle1[0], triangle1[1], triangle1[2],
            triangle2[0], triangle2[1], triangle2[2],
            sl);
    if (precom == false) {
        int tot, ori, fid;

        LPI_exact_suppvars s;
        bool premulti = orient3D_LPI_pre_exact(
            segpoint0[0],
            segpoint0[1],
            segpoint0[2],
            segpoint1[0],
            segpoint1[1],
            segpoint1[2],

            triangle0[0],
            triangle0[1],
            triangle0[2],
            triangle1[0],
            triangle1[1],
            triangle1[2],
            triangle2[0],
            triangle2[1],
            triangle2[2],
            s);


        if (premulti == false) return 2;
        for (int i = 0; i < prismindex.size(); i++) {
            if (prismindex[i] == jump) {
                continue;
            }
            if (coverlist[i] == true) continue;
            tot = 0;
            fid = 0;
            ori = -1;
            for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
                if (intersect_face[i][fid] == j) {
                    if (fid + 1 < intersect_face[i].size()) fid++;
                } else
                    continue;
                ori = orient3D_LPI_post_exact(
                    s,
                    segpoint0[0],
                    segpoint0[1],
                    segpoint0[2],

                    halfspace[prismindex[i]][j][0][0],
                    halfspace[prismindex[i]][j][0][1],
                    halfspace[prismindex[i]][j][0][2],

                    halfspace[prismindex[i]][j][1][0],
                    halfspace[prismindex[i]][j][1][1],
                    halfspace[prismindex[i]][j][1][2],

                    halfspace[prismindex[i]][j][2][0],
                    halfspace[prismindex[i]][j][2][1],
                    halfspace[prismindex[i]][j][2][2]);
                if (ori == 1 || ori == 0) {
                    break;
                }

                if (ori == -1) {
                    tot++;
                }
            }
            if (ori == 1 || ori == 0) continue;
            fid = 0;
            ori = -1;
            for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
                if (intersect_face[i][fid] == j) {
                    if (fid + 1 < intersect_face[i].size()) fid++;
                    continue;
                }

                ori = orient3D_LPI_post_exact(
                    s,
                    segpoint0[0],
                    segpoint0[1],
                    segpoint0[2],

                    halfspace[prismindex[i]][j][0][0],
                    halfspace[prismindex[i]][j][0][1],
                    halfspace[prismindex[i]][j][0][2],

                    halfspace[prismindex[i]][j][1][0],
                    halfspace[prismindex[i]][j][1][1],
                    halfspace[prismindex[i]][j][1][2],

                    halfspace[prismindex[i]][j][2][0],
                    halfspace[prismindex[i]][j][2][1],
                    halfspace[prismindex[i]][j][2][2]);
                if (ori == 1 || ori == 0) {
                    break;
                }

                if (ori == -1) {
                    tot++;
                }
            }
            if (ori == 1 || ori == 0) continue;
            if (tot == halfspace[prismindex[i]].size()) {
                id = i;
                return IN_PRISM;
            }
        }


        return OUT_PRISM;
    }
    int tot, fid;
    int ori;
    INDEX index;
    std::vector<INDEX> recompute;

    recompute.clear();
    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump) continue;
        if (coverlist[i] == true) continue;
        index.FACES.clear();
        tot = 0;
        fid = 0;
        ori = -1;

        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
            } else
                continue;
            ori = LPI_POSTFILTER(sl, segpoint0[0], segpoint0[1], segpoint0[2], prismindex[i], j);

            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                index.FACES.emplace_back(j);
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1) continue;
        fid = 0;
        ori = -1;
        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
                continue;
            }

            ori = LPI_POSTFILTER(sl, segpoint0[0], segpoint0[1], segpoint0[2], prismindex[i], j);

            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                index.FACES.emplace_back(j);
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1) continue;
        if (tot == halfspace[prismindex[i]].size()) {
            id = i;
            return IN_PRISM;
        }

        if (ori != 1) {
            assert(!index.FACES.empty());
            index.Pi = i; // local id
            recompute.emplace_back(index);
        }
    }

    if (!recompute.empty()) {
        LPI_exact_suppvars s;
        bool premulti = orient3D_LPI_pre_exact(
            segpoint0[0],
            segpoint0[1],
            segpoint0[2],
            segpoint1[0],
            segpoint1[1],
            segpoint1[2],
            triangle0[0],
            triangle0[1],
            triangle0[2],
            triangle1[0],
            triangle1[1],
            triangle1[2],
            triangle2[0],
            triangle2[1],
            triangle2[2],
            s);


        for (int k = 0; k < recompute.size(); k++) {
            int in1 = prismindex[recompute[k].Pi];

            for (int j = 0; j < recompute[k].FACES.size(); j++) {
                int in2 = recompute[k].FACES[j];

                ori = orient3D_LPI_post_exact(
                    s,
                    segpoint0[0],
                    segpoint0[1],
                    segpoint0[2],
                    halfspace[in1][in2][0][0],
                    halfspace[in1][in2][0][1],
                    halfspace[in1][in2][0][2],
                    halfspace[in1][in2][1][0],
                    halfspace[in1][in2][1][1],
                    halfspace[in1][in2][1][2],
                    halfspace[in1][in2][2][0],
                    halfspace[in1][in2][2][1],
                    halfspace[in1][in2][2][2]);

                if (ori == 1 || ori == 0) break;
            }

            if (ori == -1) {
                id = recompute[k].Pi;
                return IN_PRISM;
            }
        }
    }

    return OUT_PRISM;
}

#ifdef ENVELOPE_WITH_GMP
int FastEnvelope::
    Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over_Rational(
        const Vector3& segpoint0,
        const Vector3& segpoint1,
        const Vector3& triangle0,
        const Vector3& triangle1,
        const Vector3& triangle2,
        const std::vector<unsigned int>& prismindex,
        const std::vector<std::vector<int>>& intersect_face,
        const std::vector<bool>& coverlist,
        const int& jump,
        int& id) const
{
    Rational s00(segpoint0[0]), s01(segpoint0[1]), s02(segpoint0[2]), s10(segpoint1[0]),
        s11(segpoint1[1]), s12(segpoint1[2]), t00(triangle0[0]), t01(triangle0[1]),
        t02(triangle0[2]), t10(triangle1[0]), t11(triangle1[1]), t12(triangle1[2]),
        t20(triangle2[0]), t21(triangle2[1]), t22(triangle2[2]), a11, a12, a13, d;
    bool precom = orient3D_LPI_prefilter_multiprecision(
        s00,
        s01,
        s02,
        s10,
        s11,
        s12,
        t00,
        t01,
        t02,
        t10,
        t11,
        t12,
        t20,
        t21,
        t22,
        a11,
        a12,
        a13,
        d,
        check_rational);

    if (precom == false) {
        return 2;
    }

    int tot, fid;
    int ori;

    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump) continue;
        if (coverlist[i] == true) continue;
        tot = 0;
        fid = 0;
        ori = -1;

        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
            } else
                continue;

            Rational h00(halfspace[prismindex[i]][j][0][0]), h01(halfspace[prismindex[i]][j][0][1]),
                h02(halfspace[prismindex[i]][j][0][2]), h10(halfspace[prismindex[i]][j][1][0]),
                h11(halfspace[prismindex[i]][j][1][1]), h12(halfspace[prismindex[i]][j][1][2]),
                h20(halfspace[prismindex[i]][j][2][0]), h21(halfspace[prismindex[i]][j][2][1]),
                h22(halfspace[prismindex[i]][j][2][2]);
            ori = orient3D_LPI_postfilter_multiprecision(
                a11,
                a12,
                a13,
                d,
                s00,
                s01,
                s02,
                h00,
                h01,
                h02,
                h10,
                h11,
                h12,
                h20,
                h21,
                h22,
                check_rational);


            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                break;
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1 || ori == 0) continue;
        fid = 0;
        ori = -1;
        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
                continue;
            }
            Rational h00(halfspace[prismindex[i]][j][0][0]), h01(halfspace[prismindex[i]][j][0][1]),
                h02(halfspace[prismindex[i]][j][0][2]), h10(halfspace[prismindex[i]][j][1][0]),
                h11(halfspace[prismindex[i]][j][1][1]), h12(halfspace[prismindex[i]][j][1][2]),
                h20(halfspace[prismindex[i]][j][2][0]), h21(halfspace[prismindex[i]][j][2][1]),
                h22(halfspace[prismindex[i]][j][2][2]);
            ori = orient3D_LPI_postfilter_multiprecision(
                a11,
                a12,
                a13,
                d,
                s00,
                s01,
                s02,
                h00,
                h01,
                h02,
                h10,
                h11,
                h12,
                h20,
                h21,
                h22,
                check_rational);

            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                break;
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1 || ori == 0) continue;
        if (tot == halfspace[prismindex[i]].size()) {
            id = i;
            return IN_PRISM;
        }
    }

    return OUT_PRISM;
}
#endif

int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
    const std::array<Vector3, 3>& triangle,
    const Vector3& facet10,
    const Vector3& facet11,
    const Vector3& facet12,
    const Vector3& facet20,
    const Vector3& facet21,
    const Vector3& facet22,
    const std::vector<unsigned int>& prismindex,
    const std::vector<std::vector<int>>& intersect_face,
    const int& jump1,
    const int& jump2,
    int& id,
    const TPI_filtered_suppvars* shared_st,
    int shared_pre) const
{
    FACET_RECORD_MEMO;
    TPI_exact_suppvars s;
    TPI_filtered_suppvars st_own;

    int tot, ori, fid;
    bool pre = shared_st ? (shared_pre != 0)
                         : orient3D_TPI_prefilter(
                               triangle[0][0],
                               triangle[0][1],
                               triangle[0][2],
                               triangle[1][0],
                               triangle[1][1],
                               triangle[1][2],
                               triangle[2][0],
                               triangle[2][1],
                               triangle[2][2],
                               facet10[0],
                               facet10[1],
                               facet10[2],
                               facet11[0],
                               facet11[1],
                               facet11[2],
                               facet12[0],
                               facet12[1],
                               facet12[2],
                               facet20[0],
                               facet20[1],
                               facet20[2],
                               facet21[0],
                               facet21[1],
                               facet21[2],
                               facet22[0],
                               facet22[1],
                               facet22[2],
                               st_own);
    const TPI_filtered_suppvars& st = shared_st ? *shared_st : st_own;

    if (pre == false) {
        bool premulti = orient3D_TPI_pre_exact(
            triangle[0][0],
            triangle[0][1],
            triangle[0][2],
            triangle[1][0],
            triangle[1][1],
            triangle[1][2],
            triangle[2][0],
            triangle[2][1],
            triangle[2][2],
            facet10[0],
            facet10[1],
            facet10[2],
            facet11[0],
            facet11[1],
            facet11[2],
            facet12[0],
            facet12[1],
            facet12[2],
            facet20[0],
            facet20[1],
            facet20[2],
            facet21[0],
            facet21[1],
            facet21[2],
            facet22[0],
            facet22[1],
            facet22[2],
            s);
        if (premulti == false) return 2;
        for (int i = 0; i < prismindex.size(); i++) {
            if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
            if (!algorithms::box_box_intersection(
                    cornerlist[prismindex[i]][0],
                    cornerlist[prismindex[i]][1],
                    cornerlist[jump1][0],
                    cornerlist[jump1][1]))
                continue;
            if (!algorithms::box_box_intersection(
                    cornerlist[prismindex[i]][0],
                    cornerlist[prismindex[i]][1],
                    cornerlist[jump2][0],
                    cornerlist[jump2][1]))
                continue;

            tot = 0;
            fid = 0;
            ori = -1;
            for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
                if (intersect_face[i][fid] == j) {
                    if (fid + 1 < intersect_face[i].size()) fid++;
                } else
                    continue;

                ori = orient3D_TPI_post_exact(
                    s,
                    halfspace[prismindex[i]][j][0][0],
                    halfspace[prismindex[i]][j][0][1],
                    halfspace[prismindex[i]][j][0][2],
                    halfspace[prismindex[i]][j][1][0],
                    halfspace[prismindex[i]][j][1][1],
                    halfspace[prismindex[i]][j][1][2],
                    halfspace[prismindex[i]][j][2][0],
                    halfspace[prismindex[i]][j][2][1],
                    halfspace[prismindex[i]][j][2][2]);

                if (ori == 1 || ori == 0) {
                    break;
                }

                if (ori == -1) {
                    tot++;
                }
            }
            if (ori == 1 || ori == 0) continue;
            fid = 0;
            ori = -1;
            for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
                if (intersect_face[i][fid] == j) {
                    if (fid + 1 < intersect_face[i].size()) fid++;
                    continue;
                }


                ori = orient3D_TPI_post_exact(
                    s,
                    halfspace[prismindex[i]][j][0][0],
                    halfspace[prismindex[i]][j][0][1],
                    halfspace[prismindex[i]][j][0][2],
                    halfspace[prismindex[i]][j][1][0],
                    halfspace[prismindex[i]][j][1][1],
                    halfspace[prismindex[i]][j][1][2],
                    halfspace[prismindex[i]][j][2][0],
                    halfspace[prismindex[i]][j][2][1],
                    halfspace[prismindex[i]][j][2][2]);

                if (ori == 1 || ori == 0) {
                    break;
                }

                if (ori == -1) {
                    tot++;
                }
            }
            if (ori == 1 || ori == 0) continue;
            if (tot == halfspace[prismindex[i]].size()) {
                id = i;
                return IN_PRISM;
            }
        }

        return OUT_PRISM;
    }

    INDEX index;
    std::vector<INDEX> recompute;
    recompute.clear();


    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
        if (!algorithms::box_box_intersection(
                cornerlist[prismindex[i]][0],
                cornerlist[prismindex[i]][1],
                cornerlist[jump1][0],
                cornerlist[jump1][1]))
            continue;
        if (!algorithms::box_box_intersection(
                cornerlist[prismindex[i]][0],
                cornerlist[prismindex[i]][1],
                cornerlist[jump2][0],
                cornerlist[jump2][1]))
            continue;

        index.FACES.clear();
        tot = 0;
        fid = 0;
        ori = -1;
        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
            } else
                continue;

            ori = TPI_POSTFILTER(st, prismindex[i], j);


            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                index.FACES.emplace_back(j);
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1) continue;
        fid = 0;
        ori = -1;
        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
                continue;
            }


            ori = TPI_POSTFILTER(st, prismindex[i], j);


            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                index.FACES.emplace_back(j);
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1) continue;
        if (tot == halfspace[prismindex[i]].size()) {
            id = i;
            return IN_PRISM;
        }

        if (ori != 1) {
            index.Pi = i;
            recompute.emplace_back(index);
        }
    }

    if (recompute.size() > 0) {
        bool premulti = orient3D_TPI_pre_exact(
            triangle[0][0],
            triangle[0][1],
            triangle[0][2],
            triangle[1][0],
            triangle[1][1],
            triangle[1][2],
            triangle[2][0],
            triangle[2][1],
            triangle[2][2],

            facet10[0],
            facet10[1],
            facet10[2],
            facet11[0],
            facet11[1],
            facet11[2],
            facet12[0],
            facet12[1],
            facet12[2],

            facet20[0],
            facet20[1],
            facet20[2],
            facet21[0],
            facet21[1],
            facet21[2],
            facet22[0],
            facet22[1],
            facet22[2],
            s);
        if (premulti == false) return 2;
        for (int k = 0; k < recompute.size(); k++) {
            int in1 = prismindex[recompute[k].Pi];

            for (int j = 0; j < recompute[k].FACES.size(); j++) {
                int in2 = recompute[k].FACES[j];

                ori = orient3D_TPI_post_exact(
                    s,
                    halfspace[in1][in2][0][0],
                    halfspace[in1][in2][0][1],
                    halfspace[in1][in2][0][2],
                    halfspace[in1][in2][1][0],
                    halfspace[in1][in2][1][1],
                    halfspace[in1][in2][1][2],
                    halfspace[in1][in2][2][0],
                    halfspace[in1][in2][2][1],
                    halfspace[in1][in2][2][2]);

                if (ori == 1 || ori == 0) break;
            }
            if (ori == -1) {
                id = recompute[k].Pi;
                return IN_PRISM;
            }
        }
    }
    return OUT_PRISM;
}

#ifdef ENVELOPE_WITH_GMP
int FastEnvelope::
    Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_Rational(
        const std::array<Vector3, 3>& triangle,
        const Vector3& facet10,
        const Vector3& facet11,
        const Vector3& facet12,
        const Vector3& facet20,
        const Vector3& facet21,
        const Vector3& facet22,
        const std::vector<unsigned int>& prismindex,
        const std::vector<std::vector<int>>& intersect_face,
        const int& jump1,
        const int& jump2,
        int& id) const
{
    int tot, ori, fid;

    Rational t00, t01, t02, t10, t11, t12, t20, t21, t22,

        f100, f101, f102, f110, f111, f112, f120, f121, f122,

        f200, f201, f202, f210, f211, f212, f220, f221, f222,

        dr, n1r, n2r, n3r;

    t00 = (triangle[0][0]);
    t01 = (triangle[0][1]);
    t02 = (triangle[0][2]);
    t10 = (triangle[1][0]);
    t11 = (triangle[1][1]);
    t12 = (triangle[1][2]);
    t20 = (triangle[2][0]);
    t21 = (triangle[2][1]);
    t22 = (triangle[2][2]);

    f100 = (facet10[0]);
    f101 = (facet10[1]);
    f102 = (facet10[2]);
    f110 = (facet11[0]);
    f111 = (facet11[1]);
    f112 = (facet11[2]);
    f120 = (facet12[0]);
    f121 = (facet12[1]);
    f122 = (facet12[2]);

    f200 = (facet20[0]);
    f201 = (facet20[1]);
    f202 = (facet20[2]);
    f210 = (facet21[0]);
    f211 = (facet21[1]);
    f212 = (facet21[2]);
    f220 = (facet22[0]);
    f221 = (facet22[1]);
    f222 = (facet22[2]);

    bool premulti = orient3D_TPI_prefilter_multiprecision(
        t00,
        t01,
        t02,
        t10,
        t11,
        t12,
        t20,
        t21,
        t22,
        f100,
        f101,
        f102,
        f110,
        f111,
        f112,
        f120,
        f121,
        f122,
        f200,
        f201,
        f202,
        f210,
        f211,
        f212,
        f220,
        f221,
        f222,
        dr,
        n1r,
        n2r,
        n3r,
        check_rational);
    if (premulti == false) return 2;

    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
        if (!algorithms::box_box_intersection(
                cornerlist[prismindex[i]][0],
                cornerlist[prismindex[i]][1],
                cornerlist[jump1][0],
                cornerlist[jump1][1]))
            continue;
        if (!algorithms::box_box_intersection(
                cornerlist[prismindex[i]][0],
                cornerlist[prismindex[i]][1],
                cornerlist[jump2][0],
                cornerlist[jump2][1]))
            continue;

        tot = 0;
        fid = 0;
        ori = -1;
        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
            } else
                continue;
            Rational h00(halfspace[prismindex[i]][j][0][0]), h01(halfspace[prismindex[i]][j][0][1]),
                h02(halfspace[prismindex[i]][j][0][2]), h10(halfspace[prismindex[i]][j][1][0]),
                h11(halfspace[prismindex[i]][j][1][1]), h12(halfspace[prismindex[i]][j][1][2]),
                h20(halfspace[prismindex[i]][j][2][0]), h21(halfspace[prismindex[i]][j][2][1]),
                h22(halfspace[prismindex[i]][j][2][2]);
            ori = orient3D_TPI_postfilter_multiprecision(
                dr,
                n1r,
                n2r,
                n3r,
                h00,
                h01,
                h02,
                h10,
                h11,
                h12,
                h20,
                h21,
                h22,
                check_rational);

            if (ori == 1 || ori == 0) {
                break;
            }

            if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1 || ori == 0) continue;
        fid = 0;
        ori = -1;
        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
                continue;
            }

            Rational h00(halfspace[prismindex[i]][j][0][0]), h01(halfspace[prismindex[i]][j][0][1]),
                h02(halfspace[prismindex[i]][j][0][2]), h10(halfspace[prismindex[i]][j][1][0]),
                h11(halfspace[prismindex[i]][j][1][1]), h12(halfspace[prismindex[i]][j][1][2]),
                h20(halfspace[prismindex[i]][j][2][0]), h21(halfspace[prismindex[i]][j][2][1]),
                h22(halfspace[prismindex[i]][j][2][2]);
            ori = orient3D_TPI_postfilter_multiprecision(
                dr,
                n1r,
                n2r,
                n3r,
                h00,
                h01,
                h02,
                h10,
                h11,
                h12,
                h20,
                h21,
                h22,
                check_rational);


            if (ori == 1 || ori == 0) {
                break;
            }

            if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1 || ori == 0) continue;
        if (tot == halfspace[prismindex[i]].size()) {
            id = i;
            return IN_PRISM;
        }
    }

    return OUT_PRISM;
}

int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism_Rational(
    const std::array<Vector3, 3>& triangle,
    const Vector3& facet10,
    const Vector3& facet11,
    const Vector3& facet12,
    const Vector3& facet20,
    const Vector3& facet21,
    const Vector3& facet22,
    const std::vector<unsigned int>& prismindex,
    const int& jump1,
    const int& jump2) const
{
    int tot, ori;

    Rational t00, t01, t02, t10, t11, t12, t20, t21, t22,

        f100, f101, f102, f110, f111, f112, f120, f121, f122,

        f200, f201, f202, f210, f211, f212, f220, f221, f222,

        dr, n1r, n2r, n3r;

    t00 = (triangle[0][0]);
    t01 = (triangle[0][1]);
    t02 = (triangle[0][2]);
    t10 = (triangle[1][0]);
    t11 = (triangle[1][1]);
    t12 = (triangle[1][2]);
    t20 = (triangle[2][0]);
    t21 = (triangle[2][1]);
    t22 = (triangle[2][2]);

    f100 = (facet10[0]);
    f101 = (facet10[1]);
    f102 = (facet10[2]);
    f110 = (facet11[0]);
    f111 = (facet11[1]);
    f112 = (facet11[2]);
    f120 = (facet12[0]);
    f121 = (facet12[1]);
    f122 = (facet12[2]);

    f200 = (facet20[0]);
    f201 = (facet20[1]);
    f202 = (facet20[2]);
    f210 = (facet21[0]);
    f211 = (facet21[1]);
    f212 = (facet21[2]);
    f220 = (facet22[0]);
    f221 = (facet22[1]);
    f222 = (facet22[2]);

    bool premulti = orient3D_TPI_prefilter_multiprecision(
        t00,
        t01,
        t02,
        t10,
        t11,
        t12,
        t20,
        t21,
        t22,
        f100,
        f101,
        f102,
        f110,
        f111,
        f112,
        f120,
        f121,
        f122,
        f200,
        f201,
        f202,
        f210,
        f211,
        f212,
        f220,
        f221,
        f222,
        dr,
        n1r,
        n2r,
        n3r,
        check_rational);
    if (premulti == false) return 2;

    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
        // if (!algorithms::box_box_intersection(cornerlist[prismindex[i]][0],
        // cornerlist[prismindex[i]][1], cornerlist[jump1][0], cornerlist[jump1][1])) continue; if
        // (!algorithms::box_box_intersection(cornerlist[prismindex[i]][0],
        // cornerlist[prismindex[i]][1], cornerlist[jump2][0], cornerlist[jump2][1])) continue;

        tot = 0;
        ori = -1;

        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            Rational h00(halfspace[prismindex[i]][j][0][0]), h01(halfspace[prismindex[i]][j][0][1]),
                h02(halfspace[prismindex[i]][j][0][2]), h10(halfspace[prismindex[i]][j][1][0]),
                h11(halfspace[prismindex[i]][j][1][1]), h12(halfspace[prismindex[i]][j][1][2]),
                h20(halfspace[prismindex[i]][j][2][0]), h21(halfspace[prismindex[i]][j][2][1]),
                h22(halfspace[prismindex[i]][j][2][2]);
            ori = orient3D_TPI_postfilter_multiprecision(
                dr,
                n1r,
                n2r,
                n3r,
                h00,
                h01,
                h02,
                h10,
                h11,
                h12,
                h20,
                h21,
                h22,
                check_rational);


            if (ori == 1 || ori == 0) {
                break;
            }

            if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1 || ori == 0) continue;
        if (tot == halfspace[prismindex[i]].size()) {
            return IN_PRISM;
        }
    }

    return OUT_PRISM;
}
#endif
int FastEnvelope::Implicit_Tri_Facet_Facet_interpoint_Out_Prism(
    const std::array<Vector3, 3>& triangle,
    const Vector3& facet10,
    const Vector3& facet11,
    const Vector3& facet12,
    const Vector3& facet20,
    const Vector3& facet21,
    const Vector3& facet22,
    const std::vector<unsigned int>& prismindex,
    const int& jump1,
    const int& jump2) const
{
    FACET_RECORD_MEMO;
    TPI_exact_suppvars s;
    TPI_filtered_suppvars st;

    int tot, ori;
    bool pre = orient3D_TPI_prefilter(
        triangle[0][0],
        triangle[0][1],
        triangle[0][2],
        triangle[1][0],
        triangle[1][1],
        triangle[1][2],
        triangle[2][0],
        triangle[2][1],
        triangle[2][2],
        facet10[0],
        facet10[1],
        facet10[2],
        facet11[0],
        facet11[1],
        facet11[2],
        facet12[0],
        facet12[1],
        facet12[2],
        facet20[0],
        facet20[1],
        facet20[2],
        facet21[0],
        facet21[1],
        facet21[2],
        facet22[0],
        facet22[1],
        facet22[2],
        st);

    if (pre == false) {
        bool premulti = orient3D_TPI_pre_exact(
            triangle[0][0],
            triangle[0][1],
            triangle[0][2],
            triangle[1][0],
            triangle[1][1],
            triangle[1][2],
            triangle[2][0],
            triangle[2][1],
            triangle[2][2],
            facet10[0],
            facet10[1],
            facet10[2],
            facet11[0],
            facet11[1],
            facet11[2],
            facet12[0],
            facet12[1],
            facet12[2],
            facet20[0],
            facet20[1],
            facet20[2],
            facet21[0],
            facet21[1],
            facet21[2],
            facet22[0],
            facet22[1],
            facet22[2],
            s);
        if (premulti == false) return 2;
        for (int i = 0; i < prismindex.size(); i++) {
            if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
            // if (!algorithms::box_box_intersection(cornerlist[prismindex[i]][0],
            // cornerlist[prismindex[i]][1], cornerlist[jump1][0], cornerlist[jump1][1])) continue;
            // if (!algorithms::box_box_intersection(cornerlist[prismindex[i]][0],
            // cornerlist[prismindex[i]][1], cornerlist[jump2][0], cornerlist[jump2][1])) continue;

            tot = 0;
            ori = -1;

            for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
                ori = orient3D_TPI_post_exact(
                    s,
                    halfspace[prismindex[i]][j][0][0],
                    halfspace[prismindex[i]][j][0][1],
                    halfspace[prismindex[i]][j][0][2],
                    halfspace[prismindex[i]][j][1][0],
                    halfspace[prismindex[i]][j][1][1],
                    halfspace[prismindex[i]][j][1][2],
                    halfspace[prismindex[i]][j][2][0],
                    halfspace[prismindex[i]][j][2][1],
                    halfspace[prismindex[i]][j][2][2]);

                if (ori == 1 || ori == 0) {
                    break;
                }

                if (ori == -1) {
                    tot++;
                }
            }
            if (ori == 1 || ori == 0) continue;
            if (tot == halfspace[prismindex[i]].size()) {
                return IN_PRISM;
            }
        }

        return OUT_PRISM;
    }

    INDEX index;
    std::vector<INDEX> recompute;
    recompute.clear();


    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
        // if (!algorithms::box_box_intersection(cornerlist[prismindex[i]][0],
        // cornerlist[prismindex[i]][1], cornerlist[jump1][0], cornerlist[jump1][1])) continue; if
        // (!algorithms::box_box_intersection(cornerlist[prismindex[i]][0],
        // cornerlist[prismindex[i]][1], cornerlist[jump2][0], cornerlist[jump2][1])) continue;

        index.FACES.clear();
        tot = 0;
        ori = -1;

        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            ori = TPI_POSTFILTER(st, prismindex[i], j);
            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                index.FACES.emplace_back(j);
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1) continue;
        if (tot == halfspace[prismindex[i]].size()) {
            return IN_PRISM;
        }

        if (ori != 1) {
            index.Pi = i;
            recompute.emplace_back(index);
        }
    }

    if (recompute.size() > 0) {
        bool premulti = orient3D_TPI_pre_exact(
            triangle[0][0],
            triangle[0][1],
            triangle[0][2],
            triangle[1][0],
            triangle[1][1],
            triangle[1][2],
            triangle[2][0],
            triangle[2][1],
            triangle[2][2],

            facet10[0],
            facet10[1],
            facet10[2],
            facet11[0],
            facet11[1],
            facet11[2],
            facet12[0],
            facet12[1],
            facet12[2],

            facet20[0],
            facet20[1],
            facet20[2],
            facet21[0],
            facet21[1],
            facet21[2],
            facet22[0],
            facet22[1],
            facet22[2],
            s);
        if (premulti == false) return 2;
        for (int k = 0; k < recompute.size(); k++) {
            int in1 = prismindex[recompute[k].Pi];

            for (int j = 0; j < recompute[k].FACES.size(); j++) {
                int in2 = recompute[k].FACES[j];

                ori = orient3D_TPI_post_exact(
                    s,
                    halfspace[in1][in2][0][0],
                    halfspace[in1][in2][0][1],
                    halfspace[in1][in2][0][2],
                    halfspace[in1][in2][1][0],
                    halfspace[in1][in2][1][1],
                    halfspace[in1][in2][1][2],
                    halfspace[in1][in2][2][0],
                    halfspace[in1][in2][2][1],
                    halfspace[in1][in2][2][2]);

                if (ori == 1 || ori == 0) break;
            }
            if (ori == -1) {
                return IN_PRISM;
            }
        }
    }
    return OUT_PRISM;
}


int FastEnvelope::
    Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
        const std::array<Vector3, 3>& triangle,
        const Vector3& facet10,
        const Vector3& facet11,
        const Vector3& facet12,
        const Vector3& facet20,
        const Vector3& facet21,
        const Vector3& facet22,
        const std::vector<unsigned int>& prismindex,
        const std::vector<std::vector<int>>& intersect_face,
        const std::vector<bool>& coverlist,
        const int& jump1,
        const int& jump2,
        int& id,
        const TPI_filtered_suppvars* shared_st,
        int shared_pre) const
{
    FACET_RECORD_MEMO;
    TPI_exact_suppvars s;
    TPI_filtered_suppvars st_own;

    int tot, ori, fid;
    bool pre = shared_st ? (shared_pre != 0)
                         : orient3D_TPI_prefilter(
                               triangle[0][0],
                               triangle[0][1],
                               triangle[0][2],
                               triangle[1][0],
                               triangle[1][1],
                               triangle[1][2],
                               triangle[2][0],
                               triangle[2][1],
                               triangle[2][2],
                               facet10[0],
                               facet10[1],
                               facet10[2],
                               facet11[0],
                               facet11[1],
                               facet11[2],
                               facet12[0],
                               facet12[1],
                               facet12[2],
                               facet20[0],
                               facet20[1],
                               facet20[2],
                               facet21[0],
                               facet21[1],
                               facet21[2],
                               facet22[0],
                               facet22[1],
                               facet22[2],
                               st_own);
    const TPI_filtered_suppvars& st = shared_st ? *shared_st : st_own;

    if (pre == false) {
        bool premulti = orient3D_TPI_pre_exact(
            triangle[0][0],
            triangle[0][1],
            triangle[0][2],
            triangle[1][0],
            triangle[1][1],
            triangle[1][2],
            triangle[2][0],
            triangle[2][1],
            triangle[2][2],
            facet10[0],
            facet10[1],
            facet10[2],
            facet11[0],
            facet11[1],
            facet11[2],
            facet12[0],
            facet12[1],
            facet12[2],
            facet20[0],
            facet20[1],
            facet20[2],
            facet21[0],
            facet21[1],
            facet21[2],
            facet22[0],
            facet22[1],
            facet22[2],
            s);
        if (premulti == false) return 2;
        for (int i = 0; i < prismindex.size(); i++) {
            if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
            if (!algorithms::box_box_intersection(
                    cornerlist[prismindex[i]][0],
                    cornerlist[prismindex[i]][1],
                    cornerlist[jump1][0],
                    cornerlist[jump1][1]))
                continue;
            if (!algorithms::box_box_intersection(
                    cornerlist[prismindex[i]][0],
                    cornerlist[prismindex[i]][1],
                    cornerlist[jump2][0],
                    cornerlist[jump2][1]))
                continue;
            if (coverlist[i] == true) continue;
            tot = 0;
            fid = 0;
            ori = -1;
            for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
                if (intersect_face[i][fid] == j) {
                    if (fid + 1 < intersect_face[i].size()) fid++;
                } else
                    continue;

                ori = orient3D_TPI_post_exact(
                    s,
                    halfspace[prismindex[i]][j][0][0],
                    halfspace[prismindex[i]][j][0][1],
                    halfspace[prismindex[i]][j][0][2],
                    halfspace[prismindex[i]][j][1][0],
                    halfspace[prismindex[i]][j][1][1],
                    halfspace[prismindex[i]][j][1][2],
                    halfspace[prismindex[i]][j][2][0],
                    halfspace[prismindex[i]][j][2][1],
                    halfspace[prismindex[i]][j][2][2]);

                if (ori == 1 || ori == 0) {
                    break;
                }

                if (ori == -1) {
                    tot++;
                }
            }
            if (ori == 1 || ori == 0) continue;
            fid = 0;
            ori = -1;
            for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
                if (intersect_face[i][fid] == j) {
                    if (fid + 1 < intersect_face[i].size()) fid++;
                    continue;
                }


                ori = orient3D_TPI_post_exact(
                    s,
                    halfspace[prismindex[i]][j][0][0],
                    halfspace[prismindex[i]][j][0][1],
                    halfspace[prismindex[i]][j][0][2],
                    halfspace[prismindex[i]][j][1][0],
                    halfspace[prismindex[i]][j][1][1],
                    halfspace[prismindex[i]][j][1][2],
                    halfspace[prismindex[i]][j][2][0],
                    halfspace[prismindex[i]][j][2][1],
                    halfspace[prismindex[i]][j][2][2]);

                if (ori == 1 || ori == 0) {
                    break;
                }

                if (ori == -1) {
                    tot++;
                }
            }
            if (ori == 1 || ori == 0) continue;
            if (tot == halfspace[prismindex[i]].size()) {
                id = i;
                return IN_PRISM;
            }
        }

        return OUT_PRISM;
    }

    INDEX index;
    std::vector<INDEX> recompute;
    recompute.clear();


    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump1 || prismindex[i] == jump2) continue;
        if (!algorithms::box_box_intersection(
                cornerlist[prismindex[i]][0],
                cornerlist[prismindex[i]][1],
                cornerlist[jump1][0],
                cornerlist[jump1][1]))
            continue;
        if (!algorithms::box_box_intersection(
                cornerlist[prismindex[i]][0],
                cornerlist[prismindex[i]][1],
                cornerlist[jump2][0],
                cornerlist[jump2][1]))
            continue;
        if (coverlist[i] == true) continue;
        index.FACES.clear();
        tot = 0;
        fid = 0;
        ori = -1;
        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
            } else
                continue;

            ori = TPI_POSTFILTER(st, prismindex[i], j);


            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                index.FACES.emplace_back(j);
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1) continue;
        fid = 0;
        ori = -1;
        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j) {
                if (fid + 1 < intersect_face[i].size()) fid++;
                continue;
            }


            ori = TPI_POSTFILTER(st, prismindex[i], j);


            if (ori == 1) {
                break;
            }
            if (ori == 0) {
                index.FACES.emplace_back(j);
            }

            else if (ori == -1) {
                tot++;
            }
        }
        if (ori == 1) continue;
        if (tot == halfspace[prismindex[i]].size()) {
            id = i;
            return IN_PRISM;
        }

        if (ori != 1) {
            index.Pi = i;
            recompute.emplace_back(index);
        }
    }

    if (recompute.size() > 0) {
        bool premulti = orient3D_TPI_pre_exact(
            triangle[0][0],
            triangle[0][1],
            triangle[0][2],
            triangle[1][0],
            triangle[1][1],
            triangle[1][2],
            triangle[2][0],
            triangle[2][1],
            triangle[2][2],

            facet10[0],
            facet10[1],
            facet10[2],
            facet11[0],
            facet11[1],
            facet11[2],
            facet12[0],
            facet12[1],
            facet12[2],

            facet20[0],
            facet20[1],
            facet20[2],
            facet21[0],
            facet21[1],
            facet21[2],
            facet22[0],
            facet22[1],
            facet22[2],
            s);
        if (premulti == false) return 2;
        for (int k = 0; k < recompute.size(); k++) {
            int in1 = prismindex[recompute[k].Pi];

            for (int j = 0; j < recompute[k].FACES.size(); j++) {
                int in2 = recompute[k].FACES[j];

                ori = orient3D_TPI_post_exact(
                    s,
                    halfspace[in1][in2][0][0],
                    halfspace[in1][in2][0][1],
                    halfspace[in1][in2][0][2],
                    halfspace[in1][in2][1][0],
                    halfspace[in1][in2][1][1],
                    halfspace[in1][in2][1][2],
                    halfspace[in1][in2][2][0],
                    halfspace[in1][in2][2][1],
                    halfspace[in1][in2][2][2]);

                if (ori == 1 || ori == 0) break;
            }
            if (ori == -1) {
                id = recompute[k].Pi;
                return IN_PRISM;
            }
        }
    }

    return OUT_PRISM;
}

bool FastEnvelope::is_3_triangle_cut_pure_multiprecision(
    const std::array<Vector3, 3>& triangle,
    TPI_exact_suppvars& s)
{
    int o1, o2, o3;
    Vector3 n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];

    if (algorithms::orient_3d(n, triangle[0], triangle[1], triangle[2]) == 0) {
        // logger().debug("Degeneration happens");
        n = {{Vector3(rand(), rand(), rand())}};
    }
    Scalar t00, t01, t02, t10, t11, t12, t20, t21, t22;

    t00 = triangle[0][0];
    t01 = triangle[0][1];
    t02 = triangle[0][2];
    t10 = triangle[1][0];
    t11 = triangle[1][1];
    t12 = triangle[1][2];
    t20 = triangle[2][0];
    t21 = triangle[2][1];
    t22 = triangle[2][2];

    o1 = orient3D_TPI_post_exact(s, n[0], n[1], n[2], t00, t01, t02, t10, t11, t12);

    if (o1 == 0) return false;

    o2 = orient3D_TPI_post_exact(s, n[0], n[1], n[2], t10, t11, t12, t20, t21, t22);

    if (o2 == 0 || o1 + o2 == 0) return false;

    o3 = orient3D_TPI_post_exact(s, n[0], n[1], n[2], t20, t21, t22, t00, t01, t02);

    if (o3 == 0 || o1 + o3 == 0 || o2 + o3 == 0) return false;
    return true;
}


bool FastEnvelope::is_3_triangle_cut(
    const std::array<Vector3, 3>& triangle,
    const Vector3& facet10,
    const Vector3& facet11,
    const Vector3& facet12,
    const Vector3& facet20,
    const Vector3& facet21,
    const Vector3& facet22,
    const TPI_filtered_suppvars* shared_st,
    int shared_pre,
    const Vector3* hoisted_n)
{
#ifdef FAST_ENVELOPE_USE_TIMER
    count1++;
#endif
    TPI_filtered_suppvars st_own;
    TPI_exact_suppvars s;
    Scalar t00, t01, t02, t10, t11, t12, t20, t21, t22,

        f100, f101, f102, f110, f111, f112, f120, f121, f122,

        f200, f201, f202, f210, f211, f212, f220, f221, f222;
#ifdef FAST_ENVELOPE_USE_TIMER
    timer.start();
#endif
    bool pre = shared_st ? (shared_pre != 0)
                         : orient3D_TPI_prefilter(
                               triangle[0][0],
                               triangle[0][1],
                               triangle[0][2],
                               triangle[1][0],
                               triangle[1][1],
                               triangle[1][2],
                               triangle[2][0],
                               triangle[2][1],
                               triangle[2][2],
                               facet10[0],
                               facet10[1],
                               facet10[2],
                               facet11[0],
                               facet11[1],
                               facet11[2],
                               facet12[0],
                               facet12[1],
                               facet12[2],
                               facet20[0],
                               facet20[1],
                               facet20[2],
                               facet21[0],
                               facet21[1],
                               facet21[2],
                               facet22[0],
                               facet22[1],
                               facet22[2],
                               st_own);
    const TPI_filtered_suppvars& st = shared_st ? *shared_st : st_own;
#ifdef FAST_ENVELOPE_USE_TIMER
    time1 += timer.getElapsedTimeInSec();
#endif
    if (pre == false) {
#ifdef FAST_ENVELOPE_USE_TIMER
        extimer.start();
        count2++;
#endif
        bool premulti = orient3D_TPI_pre_exact(
            triangle[0][0],
            triangle[0][1],
            triangle[0][2],
            triangle[1][0],
            triangle[1][1],
            triangle[1][2],
            triangle[2][0],
            triangle[2][1],
            triangle[2][2],
            facet10[0],
            facet10[1],
            facet10[2],
            facet11[0],
            facet11[1],
            facet11[2],
            facet12[0],
            facet12[1],
            facet12[2],
            facet20[0],
            facet20[1],
            facet20[2],
            facet21[0],
            facet21[1],
            facet21[2],
            facet22[0],
            facet22[1],
            facet22[2],
            s);
#ifdef FAST_ENVELOPE_USE_TIMER
        time4 += extimer.getElapsedTimeInSec();
#endif
        if (premulti == false) return false;
        return is_3_triangle_cut_pure_multiprecision(triangle, s);
    }

    Vector3 n;
    if (hoisted_n) {
        // hoisted by the caller: same expression, and the caller only passes it when the
        // orient_3d check below came out non-zero, so no rand() call is skipped.
        n = *hoisted_n;
    } else {
        n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];
        if (algorithms::orient_3d(n, triangle[0], triangle[1], triangle[2]) == 0) {
            // logger().debug("Degeneration happens");
            n = {{Vector3(rand(), rand(), rand())}};
        }
    }


    bool premulti = false;
#ifdef FAST_ENVELOPE_USE_TIMER
    timer.start();
#endif
    int o1 = hoisted_n ? orient3D_TPI_postfilter_rec(
                             st,
                             t_query_triangle.edge_points[0],
                             t_query_triangle.edge_plane[0])
                       : orient3D_TPI_postfilter(
                             st,
                             n[0],
                             n[1],
                             n[2],
                             triangle[0][0],
                             triangle[0][1],
                             triangle[0][2],
                             triangle[1][0],
                             triangle[1][1],
                             triangle[1][2]);
#ifdef FAST_ENVELOPE_USE_TIMER
    count3++;
    time2 += timer.getElapsedTimeInSec();
#endif
    if (o1 == 0) {
        t00 = (triangle[0][0]);
        t01 = (triangle[0][1]);
        t02 = (triangle[0][2]);
        t10 = (triangle[1][0]);
        t11 = (triangle[1][1]);
        t12 = (triangle[1][2]);
        t20 = (triangle[2][0]);
        t21 = (triangle[2][1]);
        t22 = (triangle[2][2]);

        f100 = (facet10[0]);
        f101 = (facet10[1]);
        f102 = (facet10[2]);
        f110 = (facet11[0]);
        f111 = (facet11[1]);
        f112 = (facet11[2]);
        f120 = (facet12[0]);
        f121 = (facet12[1]);
        f122 = (facet12[2]);

        f200 = (facet20[0]);
        f201 = (facet20[1]);
        f202 = (facet20[2]);
        f210 = (facet21[0]);
        f211 = (facet21[1]);
        f212 = (facet21[2]);
        f220 = (facet22[0]);
        f221 = (facet22[1]);
        f222 = (facet22[2]);
#ifdef FAST_ENVELOPE_USE_TIMER
        extimer.start();
#endif
        premulti = orient3D_TPI_pre_exact(
            t00,
            t01,
            t02,
            t10,
            t11,
            t12,
            t20,
            t21,
            t22,
            f100,
            f101,
            f102,
            f110,
            f111,
            f112,
            f120,
            f121,
            f122,
            f200,
            f201,
            f202,
            f210,
            f211,
            f212,
            f220,
            f221,
            f222,
            s);
#ifdef FAST_ENVELOPE_USE_TIMER
        time4 += extimer.getElapsedTimeInSec();
#endif
        o1 = orient3D_TPI_post_exact(s, n[0], n[1], n[2], t00, t01, t02, t10, t11, t12);
#ifdef FAST_ENVELOPE_USE_TIMER
        count4++;
#endif
    }

    if (o1 == 0) return false;
#ifdef FAST_ENVELOPE_USE_TIMER
    timer.start();
#endif
    int o2 = hoisted_n ? orient3D_TPI_postfilter_rec(
                             st,
                             t_query_triangle.edge_points[1],
                             t_query_triangle.edge_plane[1])
                       : orient3D_TPI_postfilter(
                             st,
                             n[0],
                             n[1],
                             n[2],
                             triangle[1][0],
                             triangle[1][1],
                             triangle[1][2],
                             triangle[2][0],
                             triangle[2][1],
                             triangle[2][2]);
#ifdef FAST_ENVELOPE_USE_TIMER
    count3++;
    time2 += timer.getElapsedTimeInSec();
#endif
    if (o2 == 0) {
        if (premulti == false) {
            t00 = (triangle[0][0]);
            t01 = (triangle[0][1]);
            t02 = (triangle[0][2]);
            t10 = (triangle[1][0]);
            t11 = (triangle[1][1]);
            t12 = (triangle[1][2]);
            t20 = (triangle[2][0]);
            t21 = (triangle[2][1]);
            t22 = (triangle[2][2]);

            f100 = (facet10[0]);
            f101 = (facet10[1]);
            f102 = (facet10[2]);
            f110 = (facet11[0]);
            f111 = (facet11[1]);
            f112 = (facet11[2]);
            f120 = (facet12[0]);
            f121 = (facet12[1]);
            f122 = (facet12[2]);

            f200 = (facet20[0]);
            f201 = (facet20[1]);
            f202 = (facet20[2]);
            f210 = (facet21[0]);
            f211 = (facet21[1]);
            f212 = (facet21[2]);
            f220 = (facet22[0]);
            f221 = (facet22[1]);
            f222 = (facet22[2]);
#ifdef FAST_ENVELOPE_USE_TIMER
            extimer.start();
#endif
            premulti = orient3D_TPI_pre_exact(
                t00,
                t01,
                t02,
                t10,
                t11,
                t12,
                t20,
                t21,
                t22,
                f100,
                f101,
                f102,
                f110,
                f111,
                f112,
                f120,
                f121,
                f122,
                f200,
                f201,
                f202,
                f210,
                f211,
                f212,
                f220,
                f221,
                f222,
                s);
#ifdef FAST_ENVELOPE_USE_TIMER
            time4 += extimer.getElapsedTimeInSec();
#endif
        }
        o2 = orient3D_TPI_post_exact(s, n[0], n[1], n[2], t10, t11, t12, t20, t21, t22);
#ifdef FAST_ENVELOPE_USE_TIMER
        count4++;
#endif
    }
    if (o2 == 0 || o1 + o2 == 0) return false;
#ifdef FAST_ENVELOPE_USE_TIMER
    timer.start();
#endif
    int o3 = hoisted_n ? orient3D_TPI_postfilter_rec(
                             st,
                             t_query_triangle.edge_points[2],
                             t_query_triangle.edge_plane[2])
                       : orient3D_TPI_postfilter(
                             st,
                             n[0],
                             n[1],
                             n[2],
                             triangle[2][0],
                             triangle[2][1],
                             triangle[2][2],
                             triangle[0][0],
                             triangle[0][1],
                             triangle[0][2]);
#ifdef FAST_ENVELOPE_USE_TIMER
    count3++;
    time2 += timer.getElapsedTimeInSec();
#endif
    if (o3 == 0) {
        if (premulti == false) {
            t00 = (triangle[0][0]);
            t01 = (triangle[0][1]);
            t02 = (triangle[0][2]);
            t10 = (triangle[1][0]);
            t11 = (triangle[1][1]);
            t12 = (triangle[1][2]);
            t20 = (triangle[2][0]);
            t21 = (triangle[2][1]);
            t22 = (triangle[2][2]);

            f100 = (facet10[0]);
            f101 = (facet10[1]);
            f102 = (facet10[2]);
            f110 = (facet11[0]);
            f111 = (facet11[1]);
            f112 = (facet11[2]);
            f120 = (facet12[0]);
            f121 = (facet12[1]);
            f122 = (facet12[2]);

            f200 = (facet20[0]);
            f201 = (facet20[1]);
            f202 = (facet20[2]);
            f210 = (facet21[0]);
            f211 = (facet21[1]);
            f212 = (facet21[2]);
            f220 = (facet22[0]);
            f221 = (facet22[1]);
            f222 = (facet22[2]);
#ifdef FAST_ENVELOPE_USE_TIMER
            extimer.start();
#endif
            premulti = orient3D_TPI_pre_exact(
                t00,
                t01,
                t02,
                t10,
                t11,
                t12,
                t20,
                t21,
                t22,
                f100,
                f101,
                f102,
                f110,
                f111,
                f112,
                f120,
                f121,
                f122,
                f200,
                f201,
                f202,
                f210,
                f211,
                f212,
                f220,
                f221,
                f222,
                s);
#ifdef FAST_ENVELOPE_USE_TIMER
            time4 += extimer.getElapsedTimeInSec();
#endif
        }
        o3 = orient3D_TPI_post_exact(s, n[0], n[1], n[2], t20, t21, t22, t00, t01, t02);
#ifdef FAST_ENVELOPE_USE_TIMER
        count4++;
#endif
    }
    if (o3 == 0 || o1 + o3 == 0 || o2 + o3 == 0) return false;

    return true;
}

#ifdef ENVELOPE_WITH_GMP
bool FastEnvelope::is_3_triangle_cut_Rational(
    const std::array<Vector3, 3>& triangle,
    const Vector3& facet10,
    const Vector3& facet11,
    const Vector3& facet12,
    const Vector3& facet20,
    const Vector3& facet21,
    const Vector3& facet22)
{
    Vector3 n = (triangle[0] - triangle[1]).cross(triangle[0] - triangle[2]) + triangle[0];

    if (algorithms::orient_3d(n, triangle[0], triangle[1], triangle[2]) == 0) {
        // std::cout << "Degeneration happens" << std::endl;
        n = {{Vector3(rand(), rand(), rand())}};
    }
    Rational t00, t01, t02, t10, t11, t12, t20, t21, t22,

        f100, f101, f102, f110, f111, f112, f120, f121, f122,

        f200, f201, f202, f210, f211, f212, f220, f221, f222,

        nr0, nr1, nr2,

        dr, n1r, n2r, n3r;

    t00 = (triangle[0][0]);
    t01 = (triangle[0][1]);
    t02 = (triangle[0][2]);
    t10 = (triangle[1][0]);
    t11 = (triangle[1][1]);
    t12 = (triangle[1][2]);
    t20 = (triangle[2][0]);
    t21 = (triangle[2][1]);
    t22 = (triangle[2][2]);

    f100 = (facet10[0]);
    f101 = (facet10[1]);
    f102 = (facet10[2]);
    f110 = (facet11[0]);
    f111 = (facet11[1]);
    f112 = (facet11[2]);
    f120 = (facet12[0]);
    f121 = (facet12[1]);
    f122 = (facet12[2]);

    f200 = (facet20[0]);
    f201 = (facet20[1]);
    f202 = (facet20[2]);
    f210 = (facet21[0]);
    f211 = (facet21[1]);
    f212 = (facet21[2]);
    f220 = (facet22[0]);
    f221 = (facet22[1]);
    f222 = (facet22[2]);

    nr0 = (n[0]);
    nr1 = (n[1]);
    nr2 = (n[2]);
    bool premulti = orient3D_TPI_prefilter_multiprecision(
        t00,
        t01,
        t02,
        t10,
        t11,
        t12,
        t20,
        t21,
        t22,
        f100,
        f101,
        f102,
        f110,
        f111,
        f112,
        f120,
        f121,
        f122,
        f200,
        f201,
        f202,
        f210,
        f211,
        f212,
        f220,
        f221,
        f222,
        dr,
        n1r,
        n2r,
        n3r,
        check_rational);
    if (premulti == false) return false; // point not exist

    int o1 = orient3D_TPI_postfilter_multiprecision(
        dr,
        n1r,
        n2r,
        n3r,
        nr0,
        nr1,
        nr2,
        t00,
        t01,
        t02,
        t10,
        t11,
        t12,
        check_rational);

    if (o1 == 0) return false;

    int o2 = orient3D_TPI_postfilter_multiprecision(
        dr,
        n1r,
        n2r,
        n3r,
        nr0,
        nr1,
        nr2,
        t10,
        t11,
        t12,
        t20,
        t21,
        t22,
        check_rational);

    if (o2 == 0 || o1 + o2 == 0) return false;

    int o3 = orient3D_TPI_postfilter_multiprecision(
        dr,
        n1r,
        n2r,
        n3r,
        nr0,
        nr1,
        nr2,
        t20,
        t21,
        t22,
        t00,
        t01,
        t02,
        check_rational);

    if (o3 == 0 || o1 + o3 == 0 || o2 + o3 == 0) {
        return false;
    }

    return true;
}
#endif

int FastEnvelope::is_3_triangle_cut_float_fast(
    const Vector3& tri0,
    const Vector3& tri1,
    const Vector3& tri2,
    const Vector3& facet10,
    const Vector3& facet11,
    const Vector3& facet12,
    const Vector3& facet20,
    const Vector3& facet21,
    const Vector3& facet22)
{
    Vector3 n = (tri0 - tri1).cross(tri0 - tri2) + tri0;

    if (algorithms::orient_3d(n, tri0, tri1, tri2) == 0) {
        // logger().debug("Degeneration happens !");

        n = {{Vector3(rand(), rand(), rand())}};
    }
    TPI_filtered_suppvars st;
#ifdef FAST_ENVELOPE_USE_TIMER
    timer.start();
#endif
    bool pre = orient3D_TPI_prefilter(
        tri0[0],
        tri0[1],
        tri0[2],
        tri1[0],
        tri1[1],
        tri1[2],
        tri2[0],
        tri2[1],
        tri2[2],
        facet10[0],
        facet10[1],
        facet10[2],
        facet11[0],
        facet11[1],
        facet11[2],
        facet12[0],
        facet12[1],
        facet12[2],
        facet20[0],
        facet20[1],
        facet20[2],
        facet21[0],
        facet21[1],
        facet21[2],
        facet22[0],
        facet22[1],
        facet22[2],
        st);
#ifdef FAST_ENVELOPE_USE_TIMER
    time1 += timer.getElapsedTimeInSec();
#endif
    if (pre == false) return 2; // means we dont know
#ifdef FAST_ENVELOPE_USE_TIMER
    timer.start();
#endif
    int o1 = orient3D_TPI_postfilter(
        st,
        n[0],
        n[1],
        n[2],
        tri0[0],
        tri0[1],
        tri0[2],
        tri1[0],
        tri1[1],
        tri1[2]);
    int o2 = orient3D_TPI_postfilter(
        st,
        n[0],
        n[1],
        n[2],
        tri1[0],
        tri1[1],
        tri1[2],
        tri2[0],
        tri2[1],
        tri2[2]);
#ifdef FAST_ENVELOPE_USE_TIMER
    time2 += timer.getElapsedTimeInSec();
#endif
    if (o1 * o2 == -1) return 0;
#ifdef FAST_ENVELOPE_USE_TIMER
    timer.start();
#endif
    int o3 = orient3D_TPI_postfilter(
        st,
        n[0],
        n[1],
        n[2],
        tri2[0],
        tri2[1],
        tri2[2],
        tri0[0],
        tri0[1],
        tri0[2]);
#ifdef FAST_ENVELOPE_USE_TIMER
    time2 += timer.getElapsedTimeInSec();
#endif
    if (o1 * o3 == -1 || o2 * o3 == -1) return 0;
    if (o1 * o2 * o3 == 0) return 2; // means we dont know
    return 1;
}


int FastEnvelope::is_triangle_cut_envelope_polyhedra(
    const int& cindex, // the triangle is not degenerated
    const Vector3& tri0,
    const Vector3& tri1,
    const Vector3& tri2,
    std::vector<int>& cid) const
{
    FACET_RECORD_MEMO;
    cid.clear();
    cid.reserve(3);
    // Prisms have at most 8 facets (6 for the degenerate box; assign_instance checks the
    // bound): fixed-size scratch instead of five heap allocations per candidate prism.
    assert(halfspace[cindex].size() <= kMaxFacets);
    std::array<bool, kMaxFacets> cut;
    cut.fill(false);
    std::array<int, kMaxFacets> o1, o2, o3;
    FacetList cutp;
    int ori = 0, ct1 = 0, ct2 = 0, ct3 = 0;


    for (int i = 0; i < halfspace[cindex].size(); i++) {
        o1[i] = algorithms::orient_3d(
            tri0,
            halfspace[cindex][i][0],
            halfspace[cindex][i][1],
            halfspace[cindex][i][2]);
        o2[i] = algorithms::orient_3d(
            tri1,
            halfspace[cindex][i][0],
            halfspace[cindex][i][1],
            halfspace[cindex][i][2]);
        o3[i] = algorithms::orient_3d(
            tri2,
            halfspace[cindex][i][0],
            halfspace[cindex][i][1],
            halfspace[cindex][i][2]);
        if (o1[i] + o2[i] + o3[i] >= 2) // 1,1,0 case
        {
            return 0;
        }
        if (o1[i] == 1) ct1++;
        if (o2[i] == 1) ct2++;
        if (o3[i] == 1)
            ct3++; // if ct1 or ct2 or ct3 >0, then NOT totally inside, otherwise, totally inside
        if (o1[i] == 0 && o2[i] == 0 && o3[i] == 1) {
            return 0;
        }
        if (o1[i] == 1 && o2[i] == 0 && o3[i] == 0) {
            return 0;
        }
        if (o1[i] == 0 && o2[i] == 1 && o3[i] == 0) {
            return 0;
        }


        if (o1[i] * o2[i] == -1 || o1[i] * o3[i] == -1 || o3[i] * o2[i] == -1)
            cutp.emplace_back(i);
        else if (o1[i] + o2[i] + o3[i] == -1 && o1[i] * o2[i] == 0) { // 0,0,-1 case, we also want
                                                                      // this face,really rare to
                                                                      // happen
            cutp.emplace_back(i);
        }
    }
    for (int i = 0; i < halfspace[cindex].size(); i++) {
        t_last_orientations[0][i] = (int8_t)o1[i];
        t_last_orientations[1][i] = (int8_t)o2[i];
        t_last_orientations[2][i] = (int8_t)o3[i];
    }
    if (cutp.size() == 0) {
        if (ct1 == 0 && ct2 == 0 && ct3 == 0) {
            return 2; // totally inside, or not any edge is on the facet
        }
    }

    if (cutp.size() == 0) {
        return 0;
    }
#ifdef ENVELOPE_WITH_GMP
    cid.assign(cutp.d, cutp.d + cutp.n);
    return 1;
#endif
    LPI_filtered_suppvars sl;
    std::array<Scalar, 2> seg00, seg01, seg02, seg10, seg11, seg12;
    for (int i = 0; i < cutp.size(); i++) {
        int temp = 0;
        if (o1[cutp[i]] * o2[cutp[i]] == -1 || o1[cutp[i]] + o2[cutp[i]] == -1) {
            seg00[temp] = tri0[0];
            seg01[temp] = tri0[1];
            seg02[temp] = tri0[2];
            seg10[temp] = tri1[0];
            seg11[temp] = tri1[1];
            seg12[temp] = tri1[2];
            temp++;
        }
        if (o1[cutp[i]] * o3[cutp[i]] == -1 || o1[cutp[i]] + o3[cutp[i]] == -1) {
            seg00[temp] = tri0[0];
            seg01[temp] = tri0[1];
            seg02[temp] = tri0[2];
            seg10[temp] = tri2[0];
            seg11[temp] = tri2[1];
            seg12[temp] = tri2[2];
            temp++;
        }
        if (o2[cutp[i]] * o3[cutp[i]] == -1 || o2[cutp[i]] + o3[cutp[i]] == -1) {
            seg00[temp] = tri1[0];
            seg01[temp] = tri1[1];
            seg02[temp] = tri1[2];
            seg10[temp] = tri2[0];
            seg11[temp] = tri2[1];
            seg12[temp] = tri2[2];
            temp++;
        }

        for (int k = 0; k < 2; k++) {
            bool precom = orient3D_LPI_prefilter(
                seg00[k],
                seg01[k],
                seg02[k],
                seg10[k],
                seg11[k],
                seg12[k],
                halfspace[cindex][cutp[i]][0][0],
                halfspace[cindex][cutp[i]][0][1],
                halfspace[cindex][cutp[i]][0][2],
                halfspace[cindex][cutp[i]][1][0],
                halfspace[cindex][cutp[i]][1][1],
                halfspace[cindex][cutp[i]][1][2],
                halfspace[cindex][cutp[i]][2][0],
                halfspace[cindex][cutp[i]][2][1],
                halfspace[cindex][cutp[i]][2][2],
                sl);
            if (precom == false) {
                cut[cutp[i]] = true;
                break;
            }
            for (int j = 0; j < cutp.size(); j++) {
                if (i == j) continue;
                ori = LPI_POSTFILTER(sl, seg00[k], seg01[k], seg02[k], cindex, cutp[j]);

                if (ori == 1) break;
            }
            if (ori != 1) {
                cut[cutp[i]] = true;
                break;
            }
        }

        ori = 0; // initialize the orientation to avoid the j loop doesn't happen because
                 // cutp.size()==1
    }

    if (cutp.size() <= 2) {
        for (int i = 0; i < halfspace[cindex].size(); i++) {
            if (cut[i] == true) cid.emplace_back(i);
        }
        return 1;
    }
    // triangle-facet-facet intersection
    TPI_filtered_suppvars st;

    for (int i = 0; i < cutp.size(); i++) {
        for (int j = i + 1; j < cutp.size(); j++) // two facets and the triangle generate a point
        {
            if (cut[cutp[i]] == true && cut[cutp[j]] == true) continue;
            if (USE_ADJACENT_INFORMATION) {
                bool neib = is_two_facets_neighbouring(cindex, cutp[i], cutp[j]);
                if (neib == false) continue;
            }
#ifdef FAST_ENVELOPE_USE_TIMER
            totaltimer.start();
#endif
            // is_3_triangle_cut_float_fast, inlined so that its prefilter -- which the
            // original recomputes from the same 27 doubles right after -- lands in st.
            int inter;
            {
                const bool hoisted = t_query_triangle.valid && t_query_triangle.n_off_plane;
                Vector3 n;
                if (!hoisted) {
                    n = (tri0 - tri1).cross(tri0 - tri2) + tri0;
                    if (algorithms::orient_3d(n, tri0, tri1, tri2) == 0) {
                        n = {{Vector3(rand(), rand(), rand())}};
                    }
                }
                const bool pre_ff = orient3D_TPI_prefilter_rec(
                    t_query_triangle.plane,
                    FACET_RECORDS(cindex)[cutp[i]],
                    FACET_RECORDS_2(cindex)[cutp[j]],
                    t_query_triangle.points,
                    halfspace[cindex][cutp[i]][0].data(),
                    halfspace[cindex][cutp[j]][0].data(),
                    st);
                if (!pre_ff) {
                    inter = 2;
                } else {
                    int fo1, fo2, fo3;
                    if (hoisted) {
                        fo1 = orient3D_TPI_postfilter_rec(
                            st,
                            t_query_triangle.edge_points[0],
                            t_query_triangle.edge_plane[0]);
                        fo2 = orient3D_TPI_postfilter_rec(
                            st,
                            t_query_triangle.edge_points[1],
                            t_query_triangle.edge_plane[1]);
                    } else {
                        fo1 = orient3D_TPI_postfilter(
                            st,
                            n[0],
                            n[1],
                            n[2],
                            tri0[0],
                            tri0[1],
                            tri0[2],
                            tri1[0],
                            tri1[1],
                            tri1[2]);
                        fo2 = orient3D_TPI_postfilter(
                            st,
                            n[0],
                            n[1],
                            n[2],
                            tri1[0],
                            tri1[1],
                            tri1[2],
                            tri2[0],
                            tri2[1],
                            tri2[2]);
                    }
                    if (fo1 * fo2 == -1) {
                        inter = 0;
                    } else {
                        if (hoisted)
                            fo3 = orient3D_TPI_postfilter_rec(
                                st,
                                t_query_triangle.edge_points[2],
                                t_query_triangle.edge_plane[2]);
                        else
                            fo3 = orient3D_TPI_postfilter(
                                st,
                                n[0],
                                n[1],
                                n[2],
                                tri2[0],
                                tri2[1],
                                tri2[2],
                                tri0[0],
                                tri0[1],
                                tri0[2]);
                        if (fo1 * fo3 == -1 || fo2 * fo3 == -1)
                            inter = 0;
                        else if (fo1 * fo2 * fo3 == 0)
                            inter = 2;
                        else
                            inter = 1;
                    }
                }
            }
            if (inter == 2) { // we dont know if point exist or if inside of triangle
                cut[cutp[i]] = true;
                cut[cutp[j]] = true;
                continue;
            }
            if (inter == 0) continue; // sure not inside

            for (int k = 0; k < cutp.size(); k++) {
                if (k == i || k == j) continue;

                ori = TPI_POSTFILTER(st, cindex, cutp[k]);

                if (ori == 1) break;
            }

            if (ori != 1) {
                cut[cutp[i]] = true;
                cut[cutp[j]] = true;
            }
        }
    }

    for (int i = 0; i < halfspace[cindex].size(); i++) {
        if (cut[i] == true) cid.emplace_back(i);
    }
    if (cid.size() == 0)
        return 0; // not cut and facets, and not totally inside, then not intersected
    return 1;
}
bool FastEnvelope::is_tpp_on_polyhedra(
    const std::array<Vector3, 3>& triangle,
    const Vector3& facet10,
    const Vector3& facet11,
    const Vector3& facet12,
    const Vector3& facet20,
    const Vector3& facet21,
    const Vector3& facet22,
    const int& prismid,
    const int& faceid,
    const TPI_filtered_suppvars* shared_st,
    int shared_pre) const
{
    FACET_RECORD_MEMO;
    int ori;
    TPI_filtered_suppvars st_own;
    TPI_exact_suppvars s;
    bool premulti = false;
    bool pre = shared_st ? (shared_pre != 0)
                         : orient3D_TPI_prefilter(
                               triangle[0][0],
                               triangle[0][1],
                               triangle[0][2],
                               triangle[1][0],
                               triangle[1][1],
                               triangle[1][2],
                               triangle[2][0],
                               triangle[2][1],
                               triangle[2][2],
                               facet10[0],
                               facet10[1],
                               facet10[2],
                               facet11[0],
                               facet11[1],
                               facet11[2],
                               facet12[0],
                               facet12[1],
                               facet12[2],
                               facet20[0],
                               facet20[1],
                               facet20[2],
                               facet21[0],
                               facet21[1],
                               facet21[2],
                               facet22[0],
                               facet22[1],
                               facet22[2],
                               st_own);
    const TPI_filtered_suppvars& st = shared_st ? *shared_st : st_own;
    if (pre == true) {
        for (int i = 0; i < halfspace[prismid].size(); i++) {
            /*bool neib = is_two_facets_neighbouring(prismid, i, faceid);// this works only when the polyhedron is convex and no two neighbour facets are coplanar
            if (neib == false) continue;*/
            if (i == faceid) continue;
            ori = TPI_POSTFILTER(st, prismid, i);

            if (ori == 0) {
                if (premulti == false) {
                    premulti = orient3D_TPI_pre_exact(
                        triangle[0][0],
                        triangle[0][1],
                        triangle[0][2],
                        triangle[1][0],
                        triangle[1][1],
                        triangle[1][2],
                        triangle[2][0],
                        triangle[2][1],
                        triangle[2][2],

                        facet10[0],
                        facet10[1],
                        facet10[2],
                        facet11[0],
                        facet11[1],
                        facet11[2],
                        facet12[0],
                        facet12[1],
                        facet12[2],

                        facet20[0],
                        facet20[1],
                        facet20[2],
                        facet21[0],
                        facet21[1],
                        facet21[2],
                        facet22[0],
                        facet22[1],
                        facet22[2],
                        s);
                    if (premulti == false) return false;
                }
                ori = orient3D_TPI_post_exact(
                    s,
                    halfspace[prismid][i][0][0],
                    halfspace[prismid][i][0][1],
                    halfspace[prismid][i][0][2],
                    halfspace[prismid][i][1][0],
                    halfspace[prismid][i][1][1],
                    halfspace[prismid][i][1][2],
                    halfspace[prismid][i][2][0],
                    halfspace[prismid][i][2][1],
                    halfspace[prismid][i][2][2]);
            }

            if (ori == 1) return false;
        }

    } else {
        premulti = orient3D_TPI_pre_exact(
            triangle[0][0],
            triangle[0][1],
            triangle[0][2],
            triangle[1][0],
            triangle[1][1],
            triangle[1][2],
            triangle[2][0],
            triangle[2][1],
            triangle[2][2],

            facet10[0],
            facet10[1],
            facet10[2],
            facet11[0],
            facet11[1],
            facet11[2],
            facet12[0],
            facet12[1],
            facet12[2],

            facet20[0],
            facet20[1],
            facet20[2],
            facet21[0],
            facet21[1],
            facet21[2],
            facet22[0],
            facet22[1],
            facet22[2],
            s);
        if (premulti == false) return false;
        for (int i = 0; i < halfspace[prismid].size(); i++) {
            /*bool neib = is_two_facets_neighbouring(prismid, i, faceid);// this works only when the polyhedron is convex and no two neighbour facets are coplanar
            if (neib == false) continue;*/
            if (i == faceid) continue;
            ori = orient3D_TPI_post_exact(
                s,
                halfspace[prismid][i][0][0],
                halfspace[prismid][i][0][1],
                halfspace[prismid][i][0][2],
                halfspace[prismid][i][1][0],
                halfspace[prismid][i][1][1],
                halfspace[prismid][i][1][2],
                halfspace[prismid][i][2][0],
                halfspace[prismid][i][2][1],
                halfspace[prismid][i][2][2]);

            if (ori == 1) return false;
        }
    }
    return true;
}

#ifdef ENVELOPE_WITH_GMP
bool FastEnvelope::is_tpp_on_polyhedra_Rational(
    const std::array<Vector3, 3>& triangle,
    const Vector3& facet10,
    const Vector3& facet11,
    const Vector3& facet12,
    const Vector3& facet20,
    const Vector3& facet21,
    const Vector3& facet22,
    const int& prismid,
    const int& faceid) const
{
    int ori;
    Rational t00, t01, t02, t10, t11, t12, t20, t21, t22,

        f100, f101, f102, f110, f111, f112, f120, f121, f122,

        f200, f201, f202, f210, f211, f212, f220, f221, f222,

        nr0, nr1, nr2,

        dr, n1r, n2r, n3r;

    t00 = (triangle[0][0]);
    t01 = (triangle[0][1]);
    t02 = (triangle[0][2]);
    t10 = (triangle[1][0]);
    t11 = (triangle[1][1]);
    t12 = (triangle[1][2]);
    t20 = (triangle[2][0]);
    t21 = (triangle[2][1]);
    t22 = (triangle[2][2]);

    f100 = (facet10[0]);
    f101 = (facet10[1]);
    f102 = (facet10[2]);
    f110 = (facet11[0]);
    f111 = (facet11[1]);
    f112 = (facet11[2]);
    f120 = (facet12[0]);
    f121 = (facet12[1]);
    f122 = (facet12[2]);

    f200 = (facet20[0]);
    f201 = (facet20[1]);
    f202 = (facet20[2]);
    f210 = (facet21[0]);
    f211 = (facet21[1]);
    f212 = (facet21[2]);
    f220 = (facet22[0]);
    f221 = (facet22[1]);
    f222 = (facet22[2]);

    bool premulti = orient3D_TPI_prefilter_multiprecision(
        t00,
        t01,
        t02,
        t10,
        t11,
        t12,
        t20,
        t21,
        t22,
        f100,
        f101,
        f102,
        f110,
        f111,
        f112,
        f120,
        f121,
        f122,
        f200,
        f201,
        f202,
        f210,
        f211,
        f212,
        f220,
        f221,
        f222,
        dr,
        n1r,
        n2r,
        n3r,
        check_rational);
    if (premulti == false) return false; // point not exist

    for (int i = 0; i < halfspace[prismid].size(); i++) {
        /*bool neib = is_two_facets_neighbouring(prismid, i, faceid);// this works only when the polyhedron is convex and no two neighbour facets are coplanar
        if (neib == false) continue;*/
        if (i == faceid) continue;
        Rational h00(halfspace[prismid][i][0][0]), h01(halfspace[prismid][i][0][1]),
            h02(halfspace[prismid][i][0][2]), h10(halfspace[prismid][i][1][0]),
            h11(halfspace[prismid][i][1][1]), h12(halfspace[prismid][i][1][2]),
            h20(halfspace[prismid][i][2][0]), h21(halfspace[prismid][i][2][1]),
            h22(halfspace[prismid][i][2][2]);

        ori = orient3D_TPI_postfilter_multiprecision(
            dr,
            n1r,
            n2r,
            n3r,
            h00,
            h01,
            h02,
            h10,
            h11,
            h12,
            h20,
            h21,
            h22,
            check_rational);

        if (ori == 1) return false;
    }

    return true;
}
#endif

bool FastEnvelope::is_seg_cut_polyhedra(
    const int& cindex,
    const Vector3& seg0,
    const Vector3& seg1,
    std::vector<int>& cid) const
{
    FACET_RECORD_MEMO;
    cid.clear();
    assert(halfspace[cindex].size() <= kMaxFacets);
    std::array<bool, kMaxFacets> cut;
    cut.fill(false);
    std::array<int, kMaxFacets> o1, o2;
    int ori = 0, ct1 = 0, ct2 = 0;
    FacetList cutp;

    for (int i = 0; i < halfspace[cindex].size(); i++) {
        o1[i] = algorithms::orient_3d(
            seg0,
            halfspace[cindex][i][0],
            halfspace[cindex][i][1],
            halfspace[cindex][i][2]);
        o2[i] = algorithms::orient_3d(
            seg1,
            halfspace[cindex][i][0],
            halfspace[cindex][i][1],
            halfspace[cindex][i][2]);

        if (o1[i] + o2[i] >= 1) {
            return false;
        }

        if (o1[i] == 0 && o2[i] == 0) {
            return false;
        }

        if (o1[i] * o2[i] == -1) cutp.emplace_back(i);
        if (o1[i] == 1) ct1++;
        if (o2[i] == 1) ct2++; // if ct1 or ct2 >0, then NOT totally inside
    }
    if (cutp.size() == 0 && ct1 == 0 && ct2 == 0)
        return true; // no intersected planes, and each point is either inside of poly,
    // or on one facet, since vertices are checked, then totally inside
    if (cutp.size() == 0) {
        return false;
    }
#ifdef ENVELOPE_WITH_GMP
    cid.assign(cutp.d, cutp.d + cutp.n);
    return true;
#endif
    LPI_filtered_suppvars sf;
    for (int i = 0; i < cutp.size(); i++) {
        bool precom = orient3D_LPI_prefilter( // it is boolean maybe need considering
                seg0[0], seg0[1], seg0[2],
                seg1[0], seg1[1], seg1[2],
                halfspace[cindex][cutp[i]][0][0], halfspace[cindex][cutp[i]][0][1], halfspace[cindex][cutp[i]][0][2],
                halfspace[cindex][cutp[i]][1][0], halfspace[cindex][cutp[i]][1][1], halfspace[cindex][cutp[i]][1][2],
                halfspace[cindex][cutp[i]][2][0], halfspace[cindex][cutp[i]][2][1], halfspace[cindex][cutp[i]][2][2],
                sf);
        if (precom == false) {
            cut[cutp[i]] = true;
            continue;
        }
        for (int j = 0; j < cutp.size(); j++) {
            if (i == j) continue;

            ori = LPI_POSTFILTER(sf, seg0[0], seg0[1], seg0[2], cindex, cutp[j]);

            if (ori == 1) break;
        }
        if (ori != 1) {
            cut[cutp[i]] = true;
        }
    }

    for (int i = 0; i < halfspace[cindex].size(); i++) {
        if (cut[i] == true) cid.emplace_back(i);
    }
    if (cid.size() == 0)
        return false; // if no intersection points, and segment not totally inside, then not
                      // intersected
    return true;
}
bool FastEnvelope::point_out_prism(
    const Vector3& point,
    const std::vector<unsigned int>& prismindex,
    const int& jump) const
{
    int ori;

    for (int i = 0; i < prismindex.size(); i++) {
        if (prismindex[i] == jump) continue;

        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            ori = algorithms::orient_3d(
                halfspace[prismindex[i]][j][0],
                halfspace[prismindex[i]][j][1],
                halfspace[prismindex[i]][j][2],
                point);
            if (ori == -1 || ori == 0) {
                break;
            }
            if (j == halfspace[prismindex[i]].size() - 1) {
                return false;
            }
        }
    }

    return true;
}
bool FastEnvelope::point_out_prism_return_local_id(
    const Vector3& point,
    const std::vector<unsigned int>& prismindex,
    const int& jump,
    int& id) const
{
    Vector3 bmin, bmax;

    int ori;

    for (int i = 0; i < prismindex.size(); i++) {
        bmin = cornerlist[prismindex[i]][0];
        bmax = cornerlist[prismindex[i]][1];
        if (point[0] < bmin[0] || point[1] < bmin[1] || point[2] < bmin[2]) continue;
        if (point[0] > bmax[0] || point[1] > bmax[1] || point[2] > bmax[2]) continue;
        if (prismindex[i] == jump) continue;

        for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            ori = algorithms::orient_3d(
                halfspace[prismindex[i]][j][0],
                halfspace[prismindex[i]][j][1],
                halfspace[prismindex[i]][j][2],
                point);
            if (ori == -1 || ori == 0) {
                break;
            }
        }
        if (ori == 1) {
            id = i;
            return false;
        }
    }

    return true;
}


bool FastEnvelope::is_two_facets_neighbouring(const int& pid, const int& i, const int& j) const
{
    int facesize = halfspace[pid].size();
    if (i == j) return false;
    if (i == 0 && j != 1) return true;
    if (i == 1 && j != 0) return true;
    if (j == 0 && i != 1) return true;
    if (j == 1 && i != 0) return true;
    if (i == 0 && j == 1) return false;
    if (j == 0 && i == 1) return false;
    if (i - j == 1 || j - i == 1) return true;
    if (i == 2 && j == facesize - 1) return true;
    if (j == 2 && i == facesize - 1) return true;
    return false;
}
} // namespace fastEnvelope
