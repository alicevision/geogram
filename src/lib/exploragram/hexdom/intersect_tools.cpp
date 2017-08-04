/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2015 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact for Graphite: Bruno Levy - Bruno.Levy@inria.fr
 *  Contact for this Plugin: Nicolas Ray - nicolas.ray@inria.fr
 *
 *     Project ALICE
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs.
 *
 * As an exception to the GPL, Graphite can be linked with the following
 * (non-GPL) libraries:
 *     Qt, tetgen, SuperLU, WildMagic and CGAL
 */

#include <exploragram/hexdom/intersect_tools.h>

namespace {
    using namespace GEO;
    
    struct IndexPointCmp {
        IndexPointCmp(vector<vec3>& p_data, index_t p_dim)
   	    : dim(p_dim), data(p_data) {}
        bool operator()(const index_t A, const index_t B) {
            return data[A][dim] > data[B][dim];
        }
        index_t dim;
        vector<vec3>& data;
    };
}

namespace GEO {

    const index_t quad_rand_split[2][3] = { { 0, 1, 2 },{ 0, 2, 3 } };    // random triangulation
    const index_t quad_split[4][3] = { { 0, 1, 2 },{ 0, 2, 3 },{ 1, 2, 3 },{ 1, 3, 0 } }; // convex hull
    const index_t diamon_split[12][3] = {                                 // surface and inside
        { 4, 5, 0 },{ 4, 5, 1 },{ 4, 5, 2 },{ 4, 5, 3 },
        { 0, 1, 4 },{ 1, 2, 4 },{ 2, 3, 4 },{ 3, 0, 4 },
        { 0, 1, 5 },{ 1, 2, 5 },{ 2, 3, 5 },{ 3, 0, 5 } // orientation doesn't matter
    };

    const index_t diamon_quad_split[16][3] = {           // surface and inside
        { 4, 5, 0 },{ 4, 5, 1 },{ 4, 5, 2 },{ 4, 5, 3 },
        { 0, 1, 4 },{ 1, 2, 4 },{ 2, 3, 4 },{ 3, 0, 4 },
        { 0, 1, 5 },{ 1, 2, 5 },{ 2, 3, 5 },{ 3, 0, 5 }, // orientation doesn't matter
        { 0, 1, 2 },{ 0, 2, 3 },{ 1, 2, 3 },{ 1, 3, 0 }
    };
    

    bool BBox::intersect(const BBox& b) const {
	FOR(d, 3) {
	    if (min[d] > b.max[d] || max[d] < b.min[d]) {
		return false;
	    }
	}
	return true;
    }

    bool BBox::contains(const vec3& v) const {
	FOR(d, 3) {
	    if (min[d] > v[d] || max[d] < v[d]) {
		return false;
	    }
	}
	return true;
    }

    bool BBox::is_null() const {
	FOR(d, 3) {
	    if (max[d] - min[d] < 0) {
		return true;
	    }
	}
	return false;
    }

    void BBox::add(const BBox& b) {
	if (b.is_null()) return;
	add(b.min);
	add(b.max);
    }

    void BBox::add(const vec3& P) {
	FOR(d, 3) {
	    min[d] = std::min(min[d], P[d]);
	    max[d] = std::max(max[d], P[d]);
	}
    }

    vec3 BBox::bary() const {
	return 0.5*(min + max);
    }


    /**********************************************************/

    inline unsigned int mylog2( unsigned int x ) {
	unsigned int ans = 0 ;
	while( x>>=1 ) ans++;
	return ans ;
    }

    
    void HBoxes::init(vector<BBox>& inboxes) {
	vector<vec3> G(inboxes.size());
	tree_pos_to_org.resize(inboxes.size());
	FOR(p, G.size())  G[p] = inboxes[p].bary();
	FOR(p, G.size())  tree_pos_to_org[p] = p;
	sort(G, 0, tree_pos_to_org.size());

	offset = index_t(pow(2.0, 1.0 + mylog2(G.size()))) - 1;
	tree.resize(offset + G.size());
	FOR(i, G.size())  tree[offset + i] = inboxes[tree_pos_to_org[i]];
	for (int i = int(offset) - 1; i >= 0; i--) {
	    for (int son = 2 * i + 1; son < 2 * i + 3; son++)
		if (son < int(tree.size())) tree[i].add(tree[son]);
	}

	STAT_nb_visits = 0;
	STAT_nb_leafs = 0;
	STAT_nb_requests = 0;
    }

    void HBoxes::sort(vector<vec3> &G, index_t org, index_t dest) {

	// find the best dim to cut
	index_t dim = 2;
	BBox b;
	for (index_t i = org; i < dest; i++) b.add(G[tree_pos_to_org[i]]);
	FOR(d, 2) if (b.max[d] - b.min[d] > b.max[dim] - b.min[dim]) dim = d;
	// sort
	IndexPointCmp cmp(G, dim);
	std::sort(tree_pos_to_org.begin() + int(org), tree_pos_to_org.begin() + int(dest), cmp);
	if (dest - org <= 2) return;
	index_t m = org + index_t(pow(2.0, int(mylog2(dest - org - 1))));
	sort(G, org, m);
	sort(G, m, dest);
    }

    void HBoxes::intersect(BBox& b, vector<index_t>& primitives, index_t node) {
	if (node == 0) STAT_nb_requests++;
	geo_assert(node < tree.size());
	STAT_nb_visits++;
	if (!tree[node].intersect(b)) return;
	if (node >= offset) {
	    STAT_nb_leafs++;
	    primitives.push_back(tree_pos_to_org[node - offset]);
	} else {
	    for (index_t son = 2 * node + 1; son < 2 * node + 3; son++)
		if (son < tree.size())
		    intersect(b, primitives, son);
	}
    }
    
    /**********************************************************/

    void DynamicHBoxes::init(vector<BBox>& inboxes) {
	hbox.init(inboxes);
	moved.clear();
	movedbbox.clear();
    }

    void DynamicHBoxes::intersect(BBox& b, vector<index_t>& primitives) {
	hbox.intersect(b, primitives);
	FOR(i, moved.size()) {
	    if (movedbbox[i].intersect(b))
		primitives.push_back(moved[i]);
	}
    }

    void DynamicHBoxes::update_bbox(index_t id, BBox b) {
	geo_assert(id<hbox.tree_pos_to_org.size());
	FOR(i, moved.size()) {
	    if (moved[i] == id) {
		movedbbox[i] = b;
		return;
	    }
	}
	moved.push_back(id);
	movedbbox.push_back(b);
    }

    
}
