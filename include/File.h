/*
 *  Copyright (c) 2016 Vijay Ingalalli
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef FILE_H
#define FILE_H

#include <limits>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <deque>
#include <array>
#include <algorithm>
#include <functional>
#include <sstream>
#include <set>
#include <string>
#include <cstdio>
#include <array>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bitset>
#include <math.h>
#include <iterator>

//#define DIM 700
#define MAX_EMB 100
#define EPSILON 0.001

typedef std::vector<int> Vector;
typedef std::set<std::set<int>> Set2D;
typedef std::vector<std::vector<int>> Vector2D;
typedef std::deque<std::deque<int>> Deque2D;
typedef std::set<std::set<int>> Set2D;
typedef std::vector<std::string> VecOfStr;
typedef std::vector<Vector2D> Vector3D;
typedef std::vector<std::set<int>> VecOfSet;
typedef std::vector<std::vector<std::set<int>>> EdgeLabel;
typedef std::map<std::pair<int, int>, std::set<int>> EdgeLabelMap;
typedef std::map<int, std::set<int>> NodeLabelMap;
//typedef std::map< std::pair<int, int>, std::bitset<DIM>> BitLabelMap;
//typedef std::vector<std::vector<std::bitset<DIM>>> EdgeLabelBit;
typedef std::pair<EdgeLabelMap, NodeLabelMap> SubGraph;
typedef std::unordered_map<int, std::set<int>> AttMap;
typedef std::pair<int, int> PairInt;

/*
extern "C" {
  #include "../include/nauty/traces.h"
  #include "../include/nauty/nausparse.h"
}
*/

struct GraphParameter {
    int nNodes;
    std::vector<int> nodes;
    VecOfSet attributes;
    bool nodeLabelsExist;
    GraphParameter() : edges(2) {}
    Vector2D edges;
    EdgeLabel neighbourSign;
    Vector2D adjacencyList;
    EdgeLabelMap eLabelMap;
    std::map<std::set<int>, int> edgeFrequency; // needed for the data graph |
    // Below 2 variables are for query decomposition into core and satellite nodes.
    std::set<int> coreNodes;
    std::map<int,int> old_new, new_old; // original-query-node, new-query-node ids |
    std::set<std::pair<int, int>> src_dst_map;
};

#include <../sumgra/Trie.h>
using namespace std;

class Trie; // Necessary to declare a class within a struct |

struct TimeEval{
    double SgExt = 0;
    double SgFrq = 0;
    double EmbEv = 0;
    double time1 = 0;
    double time2 = 0;
    double time3 = 0;

    int prevISOyes = 0, prevISOno = 0;
    int currISOyes = 0, currISOno = 0;
    int nonFrequent = 0;
    int starPat = 0, normalPat = 0, complexPat = 0;
};


struct Subgraph{
    EdgeLabelMap edgeLabels;
    NodeLabelMap nodeLabels;
    std::vector<int> orderedNodes;
    Vector2D embeddings;
    std::vector<int> canonicalForm;
  //    SG_DECL(canonicalNauty);
    Set2D multiEdges; // set collection of the set of all multiedges of each generated frequent subgraph |
    int frequency;
    GraphParameter parameters;
    bool decomposedPattern = false;
    Set2D nonFreqEdges;
    bool alreadyExtended = false;
    std::map<int, Set2D> nonFreqEdgesN;
    std::map<PairInt, Set2D> nonFreqEdgesE;
    bool nonFreqEdgeExists = false;
    Trie* nonFreqEdge;
    std::map<int, Trie*> nonFreqEdgeN;
    std::map<PairInt, Trie*> nonFreqEdgeE;
    TimeEval timeEval;
    std::map<int, std::vector<int>> invalidCands;
  //    IndexType sgIndex;
};

typedef std::deque<std::deque<Subgraph>> StackFSG;
typedef std::map<int, std::deque<Subgraph>> NFSG;
typedef std::map<int, std::deque<Subgraph>> FSGs;
typedef std::map<int, std::vector<int>> InvalidCands;
typedef std::vector<std::unordered_set<int>>  NodeRepository;

struct FSG{
    EdgeLabelMap structure;
    NodeLabelMap nodeLabels;
};

struct IsoTestMap{ // These maps are maintained for the entire multi-edge set and vertex-labels set for the data graph |
    std::map<std::set<int>, int> multiedgeMap;
    std::map<std::set<int>, int> nodeLabelMap;
};


template <typename T>
std::vector<T> sortIndexDecr(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<T> idx(v.size());
    for (int i = 0; i != idx.size(); ++i) idx[i] = i;

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] > v[i2];});

    return idx;
}

template <typename T>
std::vector<int> sortIndexIncr(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<int> idx(v.size());
    for (int i = 0; i != idx.size(); ++i) idx[i] = i;

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] < v[i2];});

    return idx;
}


template <typename T> /// Remove 'b' from 'a'.
std::vector<T> findDifference(std::vector<T> a, std::vector<T> b)
{
    std::unordered_set<T> bTemp(b.begin(), b.end());
    std::vector<T> solution;
    for(auto it = a.begin(); it != a.end(); ++it)
        if(bTemp.find(*it) == bTemp.end())
            solution.push_back(*it);
    return solution;
}

template <typename T> /// Merge 'b' and 'a' and store in b.
std::vector<T> vectorUnion(std::vector<T> a, std::vector<T> b)
{
    std::unordered_set<T> bTemp(b.begin(), b.end());
    for(auto it = a.begin(); it != a.end(); ++it)
        if(bTemp.find(*it) == bTemp.end())
            b.push_back(*it);
    return b;
}


template <typename T> /// Remove 'b' from 'a'.
std::vector<T> findDifferenceSet(std::vector<T> a, std::unordered_set<T> b)
{
    std::vector<T> solution;
    for(auto it = a.begin(); it != a.end(); ++it)
        if(b.find(*it) == b.end())
            solution.push_back(*it);
    return solution;
}

template <typename T>
std::vector<T> setUnion(std::vector<T> a, std::vector<T> b)
{
    std::vector<T> c(a.size()+b.size());
    typename std::vector<T>::iterator it;

    std::sort (a.begin(), a.end());
    std::sort (b.begin(), b.end());

    it=std::set_union (a.begin(), a.end(), b.begin(), b.end(), c.begin());
    c.resize(it-c.begin());
    return c;
}

template <typename T>
std::vector<T> setIntersection(std::vector<T> a, std::vector<T> b)
{
    std::vector<T> c(a.size()+b.size());
    typename std::vector<T>::iterator it;

    std::sort (a.begin(), a.end());
    std::sort (b.begin(), b.end());

    it=std::set_intersection (a.begin(), a.end(), b.begin(), b.end(), c.begin());
    c.resize(it-c.begin());
    return c;
}

class File
{
    public:
        File();
        virtual ~File();
        void readContents(const std::string& nodeFile, const std::string& edgeFile, GraphParameter& graphInfo);
        void extractInfo(const NodeLabelMap& nodes, const EdgeLabelMap& edges, GraphParameter& graphInfo);
        void printFSG(const std::map<int, std::deque<FSG>> &allFSG, std::string &filePath);
        void printOutput(int &support, int &noOfFsg, double totalTime, std::string &filePath);
    protected:
        void splitString(const std::string& str, char chr, std::vector<std::string>& strs);
};

#endif // FILE_H
