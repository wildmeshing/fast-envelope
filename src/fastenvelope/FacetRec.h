#pragma once
// The sub-expressions of the tier-1 TPI / LPI filters (indirectPredicates/ip_filtered.h) that depend
// on one plane only, computed once per plane by make_facet_rec and reused by the *_rec filters.
struct FacetRec
{
    double tc[3]; // TPI postfilter cofactors
    double lc[3]; // LPI postfilter cofactors (rounded products, then difference)
    double mx, my, mz; // max |q1-q3|, |q2-q3| per axis
    double nw[3]; // TPI prefilter plane normal
    double pdot; // nw . p1
    double mwx, mwy, mwz, mwall; // max |p2-p1|, |p3-p2| per axis, and overall
};
