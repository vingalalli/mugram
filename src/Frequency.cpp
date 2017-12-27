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

#include "Frequency.h"

Frequency::Frequency()
{
    //ctor
}

Frequency::~Frequency()
{
    //dtor
}

bool Frequency::isFrequent(const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, Subgraph& newSubgraph)
{
    SubgraphIso subIso;
    if(subIso.findEmbeddingsMNI(supportVal, newSubgraph.parameters, dataGraphInfo, graphIndexes, newSubgraph))
        return true;
    else
        return false;
}

bool containsEdge(Subgraph& A, Subgraph& B) /// Check if all the multiedges of A are contained in at least one multiedge of B
{
    bool contains = false;
    for(auto itA = A.multiEdges.begin(); itA != A.multiEdges.end(); ++itA){
        contains = false;
        for(auto itB = B.multiEdges.begin(); itB != B.multiEdges.end(); ++itB){
            if(std::includes((*itB).begin(), (*itB).end(), (*itA).begin(), (*itA).end())){
                contains = true;
                break;
            }
        }
        if(!contains)
            return false;
    }
    return contains;
}

bool Frequency::isPreviousNonFrequent(NFSG& nonFSG, Subgraph& newSubgraph, IndexType& newSubgraphInd)
{
    SubgraphIso subIso;
    for(auto it = nonFSG.begin(); it != nonFSG.end(); ++it){
    /// OPT-1: Can not be a subgraph |
        if(it->first >= newSubgraph.edgeLabels.size())
            break;
        for(auto itS = it->second.begin(); itS != it->second.end(); ++itS){
    /// OPT-2: A multiedge of NFSG is not contained in new subgraph.
            if(!containsEdge(*itS, newSubgraph))
                continue;
            Vector2D solns;
            if(subIso.findIsoMatch(*itS, newSubgraphInd, solns))
                return true;
        }
    }
    return false;
}

bool Frequency::isAlreadyFSG(std::deque<IndexType>& maximalFSG, Subgraph& newSubgraph)
{
    for(size_t i = 0; i < maximalFSG.size(); ++i){
        if(newSubgraph.parameters.nNodes < maximalFSG[i].neighborTrie.size()){
            Vector2D solns;
            SubgraphIso subIso;
            if(subIso.findIsoMatch(newSubgraph, maximalFSG[i], solns))
                return true;
        }
    }
    return false;
}
