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

#ifndef GRAPHISO_H
#define GRAPHISO_H

#include "File.h"

/// Takes a template argument of any type 'T' and the nodes of type 'int', which are partitioned into 'partitionedNodes'.
template <typename T>
void partitionNodes(const std::map<int, T>& invariant, Vector2D& partitionedNodes)
{
    std::map<T, std::vector<int>> nodeGroup;
    for(auto it = invariant.begin(); it != invariant.end(); ++it){
        auto it_g = nodeGroup.find(it->second);
        if(it_g == nodeGroup.end()){
            std::vector<int> tmp;
            tmp.push_back(it->first);
            nodeGroup.insert(std::make_pair(it->second, tmp));
        }
        else
            it_g->second.push_back(it->first);
    }
    /// Collect the partitions with decreasing vertex degrees.
    for(auto it = nodeGroup.rbegin(); it != nodeGroup.rend(); ++it)
        partitionedNodes.push_back(it->second);
}

class GraphIso
{
    public:
        GraphIso();
        virtual ~GraphIso();
        void generateCanonical(const GraphParameter& subgraph, std::vector<int>& canonicalForm, const IsoTestMap& dataGraph);
    protected:
    private:
        void computePermutationGroups(const GraphParameter& pattern, const IsoTestMap& dataGraph, Vector2D& permutationGroups);
        void computeVertexPermutation(Vector2D& permutationGroups, Vector2D& vertexPermutations);
};

#endif // GRAPHISO_H
