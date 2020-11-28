#include "../math/math.hpp"
#include "../math/rnd.hpp"
#include "variations.hpp"
#include "extrapolations.hpp"
#include "truncations.hpp"
#include <string>
#include <vector>
#include <map>

#pragma once

namespace aural {
namespace flame {

typedef math::vec(*Variation)(math::vec, double);
typedef math::vec(*Extrapolation)(double, double);
typedef double(*Truncation)(math::vec);
typedef double(*Transformation)(double);

std::map<std::string, Variation> variations;
std::map<std::string, Extrapolation> extrapolations;
std::map<std::string, Truncation> truncations;
std::map<std::string, Transformation> transformations;

std::string randVariation() {
    std::vector<std::string> t;
    for(auto i : variations)
        t.push_back(i.first);
    return t[rnd::random(0, t.size()-1)];
}

std::string randExtrapolation() {
    std::vector<std::string> t;
    for(auto i : extrapolations)
        t.push_back(i.first);
    return t[rnd::random(0, t.size()-1)];
}

std::string randTruncation() {
    std::vector<std::string> t;
    for(auto i : truncations)
        t.push_back(i.first);
    return t[rnd::random(0, t.size()-1)];
}

std::string randTransformation() {
    std::vector<std::string> t;
    for(auto i : transformations)
        t.push_back(i.first);
    return t[rnd::random(0, t.size()-1)];
}

char randOperation() {
    std::vector<char> t = {'+', '-', '*', /* '/' */};
    return t[rnd::random(0, t.size()-1)];
}

math::vec addF(math::vec v1, math::vec v2) { return math::vec(v1.x+v2.x, v1.y+v2.y); }
math::vec subF(math::vec v1, math::vec v2) { return math::vec(v1.x-v2.x, v1.y-v2.y); }
math::vec mulF(math::vec v1, math::vec v2) { return math::vec(v1.x*v2.x, v1.y*v2.y); }
math::vec divF(math::vec v1, math::vec v2) { return math::vec(v2.x == 0 ? 0 : v1.x/v2.x, v2.y == 0 ? 0 : v1.y/v2.y); }

void initFlame() {
    discard_chance = rnd::random(0, 100);
    trochoid_b = rnd::random(1.0, M_PI);
    trochoid_h = rnd::random(1.0, 6.0);
    A = rnd::random(-1.5, 1.5);
    B = rnd::random(-1.5, 1.5);
    C = rnd::gaussian(0.5, 0.2)*(rnd::random(0, 100) < 50 ? 1 : -1);
    D = rnd::random(-1.5, 1.5);
    E = rnd::random(-1.5, 1.5);
    F = rnd::random(-1.5, 1.5);
    fan_t = M_PI*C*C;
    curl_p1 = rnd::random(-0.3, 0.3), curl_p2 = rnd::random(-0.3, 0.3);
    ecycloid_a = rnd::random(1.0, 4.0), ecycloid_b = rnd::random(1.0, 3.0);
    hcycloid_a = rnd::random(sqrt(2), 1.0), hcycloid_b = rnd::random(1.0/3.0, 3.0);

    transformations["sin"] = sin;
    transformations["cos"] = cos;
    transformations["tan"] = tan;
    transformations["exp"] = exp;
    transformations["log"] = log;
    transformations["sqrt"] = sqrtfix;

    truncations["arctangent"] = arctangent;
    truncations["xpowy"] = xpowy;
    truncations["xmody"] = xmody;
    truncations["xdivy"] = xdivy;
    truncations["xmuly"] = xmuly;
    truncations["perlin"] = perlin;

    extrapolations["ranunculoid"] = ranunculoid;
    extrapolations["fishcurve"] = fishcurve;
    extrapolations["cornoid"] = cornoid;
    extrapolations["neoid"] = neoid;
    //extrapolations["maltesecross"] = maltesecross;
    extrapolations["lituus"] = lituus;
    //extrapolations["kilroy"] = kilroy;
    extrapolations["quadratrix"] = quadratrix;
    extrapolations["epitrochoid"] = epitrochoid;
    extrapolations["folium"] = folium;
    extrapolations["cycloid"] = cycloid;
    extrapolations["archispiral"] = archispiral;
    extrapolations["cissoid"] = cissoid;
    extrapolations["kampyle"] = kampyle;
    extrapolations["astroid"] = astroid;
    extrapolations["superformula"] = superformula;
    extrapolations["cartesianXY"] = cartesianXY;
    extrapolations["polarXY"] = polarXY;
    extrapolations["hyperbole"] = hyperbole;
    extrapolations["copy"] = copy;
    extrapolations["sextic"] = sextic;
    extrapolations["eightcurve"] = eightcurve;
    extrapolations["epicycloid"] = epicycloid;
    extrapolations["hypocycloid"] = hypocycloid;
    extrapolations["hypotrochoid"] = hypotrochoid;
    extrapolations["trisectrix"] = trisectrix;
    extrapolations["quadrifolium"] = quadrifolium;
    extrapolations["tschirnhausen"] = tschirnhausen;

    variations["rectangles"] = rectangles;
    variations["curl"] = curl;
    //variations["gaussian"] = gaussian;
    variations["fan"] = fan;
    variations["rings"] = rings;
    variations["cosine"] = cosine;
    variations["julia"] = julia;
    variations["sinusoidal"] = sinusoidal;
    variations["sech"] = sech;
    variations["hyperbolic"] = hyperbolic;
    variations["pdj"] = pdj;
    variations["spherical"] = spherical;
    variations["linear"] = linear;
    variations["swirl"] = swirl;
    variations["polar"] = polar;
    variations["heart"] = heart;
    variations["disc"] = disc;
    variations["spiral"] = spiral;
    variations["ex"] = ex;
    variations["bent"] = bent;
    variations["waves"] = waves;
    variations["fisheye"] = fisheye;
    variations["popcorn"] = popcorn;
    variations["exponential"] = exponential;
    variations["horseshoe"] = horseshoe;
    variations["handkerchief"] = handkerchief;
    variations["diamond"] = diamond;
    variations["power"] = power;
    variations["bubble"] = bubble;
    variations["tangent"] = tangent;
    //variations["square"] = square;
    variations["cross"] = cross;
    variations["blob"] = blob;
    variations["ngon"] = ngon;
    variations["secant"] = secant;

    //chimeras
    variations["hyperpdj"] = hyperpdj;
    variations["hyperjulia"] = hyperjulia;
    variations["sechpdj"] = sechpdj;
    variations["julofhyp"] = julofhyp;
    //variations["hyperten"] = hyperten;
    variations["spherewaves"] = spherewaves;
    variations["juspiral"] = juspiral;
    variations["spherecorn"] = spherecorn;
    variations["nspiral"] = nspiral;
    variations["wavecant"] = wavecant;
}

struct Node {
    int type;
    char operation;
    double weight;
    std::string var;
    Node() {};
    Node(int t, char o, std::string v, double w) : type(t), operation(o), var(v), weight(w) {}
};

struct Formula {
    std::vector<std::vector<int> > tree;
    std::vector<Node> nodes;
    int node_count, max_depth, max_node_count;
    int variation_chance;
    double weight_floor, weight_ceiling;
    Formula() : max_depth(2), max_node_count(5), variation_chance(50), weight_floor(1), weight_ceiling(0) {}
    Formula(int md, int mnc, int vc, double wf, double wc) : max_depth(md), max_node_count(mnc), variation_chance(vc), weight_floor(wf), weight_ceiling(wc) {} 

    void addNode(int u, int depth) { //0 - trunc, 1- var, 2 - ext, 3 - trans, 4 - op
        if(depth > max_depth || (nodes[u].type <= 1 && (rnd::random(0, 100) < math::map(depth, 0, max_depth, 0, 100)))) return;
        for(int i = 1; i <= (nodes[u].type == 4 ? 2 : 1); i++) {
            tree[u].push_back(++node_count);
            int ntype;
            if(depth >= max_depth || node_count >= max_node_count) {
                if(nodes[u].type <= 1 || nodes[u].type == 4) ntype = ((rnd::random(0, 100) <= variation_chance || nodes[u].type == 4) ? 1 : 2);
                else ntype = 0;
            } else {
                if(nodes[u].type <= 1 || nodes[u].type == 4) ntype = (rnd::random(0, 100) < 20 ? ((rnd::random(0, 100) <= variation_chance || nodes[u].type == 4) ? 1 : 2) : 4);
                else ntype = rnd::random(0, 100) < 50 ? 0 : 3;
            }
            nodes[node_count] = Node(
                ntype, 
                randOperation(), 
                ntype == 0 ? (depth == max_depth || node_count >= max_node_count) ? "perlin" : randTruncation() : ntype == 1 ? randVariation() : ntype == 2 ? randExtrapolation() : randTransformation(),
                pow(10, rnd::random(weight_floor, weight_ceiling))
            );
            
            addNode(node_count, depth+1);
        }
    }

    void create() {
        node_count = 1;
        for(int i = 0; i < 128; i++) tree.push_back({}), nodes.push_back(Node());
        int type = (max_depth == -1 ? 1 : rnd::random(0, 100) < 50 ? (rnd::random(0, 100) <= variation_chance ? 1 : 2) : 4);
        Node n(
            type, 
            randOperation(), type == 1 ? randVariation() : randExtrapolation(),
            pow(10, rnd::random(weight_floor, weight_ceiling))
        );
        nodes[node_count] = n;
        addNode(1, 0);
    }

    void create(int maxd, int maxndc, int varc, double wf, double wc) {
        max_depth = maxd;
        max_node_count = maxndc;
        variation_chance = varc;
        weight_floor = wf;
        weight_ceiling = wc;
        create();
    }

    math::vec resolve(math::vec v) {
        return resolveAux(1, v);
    }

    math::vec resolveAux(int u, math::vec v) { //0 - trunc, 1- var, 2 - ext, 3 - trans, 4 - op
        Node U = nodes[u];
        if(!tree[u].size()) {
            if(U.type == 0) return math::vec(truncations[U.var](v), 0);
            else return variations[U.var](v, U.weight);
        } else if(U.type < 4) {
            switch(U.type) {
                case 0:
                    return math::vec(truncations[U.var](resolveAux(tree[u][0], v)), 0);
                    break;
                case 1:
                    return variations[U.var](resolveAux(tree[u][0], v), U.weight);
                    break;
                case 2:
                    return extrapolations[U.var](resolveAux(tree[u][0], v).x, U.weight);
                    break;
                case 3:
                    return math::vec(transformations[U.var](resolveAux(tree[u][0], v).x), 0);
                    break;
                default:
                    return math::vec(0, 0);
                    break;
            }
        } else {
            switch(U.operation) {
                case '+':
                    return (resolveAux(tree[u][0], v) + resolveAux(tree[u][1], v));
                    break;
                case '-':
                    return (resolveAux(tree[u][0], v) - resolveAux(tree[u][1], v));
                    break;
                case '*':
                    return (resolveAux(tree[u][0], v) * resolveAux(tree[u][1], v));
                    break;
                case '/':
                    return (resolveAux(tree[u][0], v) / resolveAux(tree[u][1], v));
                    break;
                default:
                    return math::vec(0, 0);
                    break;
            }
        }
    }

    void stringifyAux(int u, std::string &s) {
        Node U = nodes[u];
        if(!tree[u].size()) s += U.var;
        else if(U.type < 4) {
            s += U.var;
            s += "(";
            stringifyAux(tree[u][0], s);
            s += ")";
        } else {
            switch(U.operation) {
                case '+':
                    s += "(";
                    stringifyAux(tree[u][0], s);
                    s += ") + (";
                    stringifyAux(tree[u][1], s);
                    s += ")";
                    break;
                case '-':
                    s += "(";
                    stringifyAux(tree[u][0], s);
                    s += ") - (";
                    stringifyAux(tree[u][1], s);
                    s += ")";
                    break;
                case '*':
                    s += "(";
                    stringifyAux(tree[u][0], s);
                    s += ") * (";
                    stringifyAux(tree[u][1], s);
                    s += ")";
                    break;
                case '/':
                    s += "(";
                    stringifyAux(tree[u][0], s);
                    s += ") / (";
                    stringifyAux(tree[u][1], s);
                    s += ")";
                    break;
            }
        }
    }

    std::string stringify() {
        std::string res = "";
        stringifyAux(1, res);
        return res;
    }
};

}
}