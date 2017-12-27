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

#include "SubgraphIso.h"

SubgraphIso::SubgraphIso()
{
    //ctor
}

SubgraphIso::~SubgraphIso()
{
    //dtor
}


void SubgraphIso::buildIndexes(GraphParameter& dataGraphInfo, IndexType& indexType)
{
    std::vector<int> srt_d_nodes(dataGraphInfo.nNodes);
    Index dataIndex;
    dataIndex.sortSignature(dataGraphInfo.neighbourSign, srt_d_nodes, dataGraphInfo.nNodes); // sort the neighborhood signature with decreasing size of the number of neighbors |

    dataIndex.buildAttHash(dataGraphInfo.attributes, indexType.attributeHash);
    dataIndex.buildSynTrie(dataGraphInfo.neighbourSign, dataGraphInfo.nNodes, indexType.synopsesTrie);

//    EdgeLabelBit dataBitSet(dataGraphInfo.nNodes);
//    dataIndex.BuildBitSign(dataGraphInfo.neighbourSign, dataGraphInfo.nNodes, dataBitSet);

    indexType.neighborTrie.resize(dataGraphInfo.nNodes);
    dataIndex.buildNeighTrie(dataGraphInfo.adjacencyList, dataGraphInfo.eLabelMap, dataGraphInfo.nNodes, indexType.neighborTrie);

}

bool SubgraphIso::findIsoMatch(Subgraph& query, IndexType& indexType, Vector2D& solutions)
{
    const GraphParameter queryGraphInfo = query.parameters;
    Match subGraph;

    /// "orderedNodes"  WILL HAVE BEEN computed, when called from "isPreviousNonFrequent".
    std::vector<int> queryOrdering = query.orderedNodes;

    /// "orderedNodes" MAY/MAY NOT HAVE BEEN computed, when called from "isAlreadyFSG".
    if(queryOrdering.empty())
        subGraph.orderVertices(queryGraphInfo.nNodes, queryGraphInfo.neighbourSign, queryGraphInfo.adjacencyList, queryOrdering);

    Index dataIndex;

    /// Fetch the vertex attribute solutions for all the query vertices to be used when necessary.

    /// With vertex labels.
    VecOfSet nodeMatches(queryGraphInfo.nNodes);

    /// Without vertex labels

    std::vector<int> initialMatches;
    std::vector<int> solution(queryGraphInfo.nNodes);
    int initialVertex = queryOrdering[0];
    if (!queryGraphInfo.neighbourSign[initialVertex].empty())
        dataIndex.querySynTrie(queryGraphInfo.neighbourSign[initialVertex], indexType.synopsesTrie, initialMatches);
        if (!initialMatches.empty()){
            if(subGraph.checkSI(initialMatches, queryGraphInfo, queryOrdering, solution, nodeMatches, indexType.neighborTrie)){
                solutions.push_back(solution);
                return true;
            }
        }
    return false;
}


bool findInitialMatches(const int& freq, IndexType& indexType, std::map<int,int>& new_old, Vector2D& initialMatchesAll, Vector2D& queryOrderPerms, Subgraph& matchedSubgraph, const EdgeLabel& neighbourSign)
{
    Index dataIndex;
    int nNodes = queryOrderPerms[0].size();
    for(size_t i = 0; i < nNodes; ++i){
        int initialVertex;
        if(new_old.empty())
            initialVertex = queryOrderPerms[i][0];
        else
            initialVertex = new_old.find(queryOrderPerms[i][0])->second;
        if (!neighbourSign[initialVertex].empty())
            dataIndex.querySynTrie(neighbourSign[initialVertex], indexType.synopsesTrie, initialMatchesAll[i]);
        if(initialMatchesAll[i].size() < freq)  // Check the size of initial matches BEFORE PRUNING |
            return false;
        else{
            auto it = matchedSubgraph.invalidCands.find(initialVertex); // check if the initial query vertex of i^th permutation  has any invalid candidates |
            if(it != matchedSubgraph.invalidCands.end()){
                    initialMatchesAll[i] = findDifference(initialMatchesAll[i], it->second); // Prune the invalid candidates, learned from the previous levels |
                if(initialMatchesAll[i].size() < freq) // Check the size of initial matches AFTER PRUNING |
                    return false;
            }
        }
    }
    return true;
}


/// Outputs only one embedding per iteration.
bool SubgraphIso::findEmbeddingsMNI(const int& freq, GraphParameter& queryGraphInfo, const GraphParameter& dataGraphInfo, IndexType& indexType, Subgraph& matchedSubgraph)
{
    int nCoreNodes = 0;
    for (size_t m = 0; m < queryGraphInfo.adjacencyList.size(); ++m){
        if (queryGraphInfo.adjacencyList[m].size() > 1)
            ++nCoreNodes ;
    }

    Match subGraph;
    Index dataIndex;

    /// With vertex labels.
    VecOfSet nodeMatches(queryGraphInfo.nNodes);
    std::vector<std::unordered_set<int>> nodeRepository(queryGraphInfo.nNodes); // Maintains distinct node embeddings for each query node.
    int requiredEmbs = 0;

    if(false){
        std::vector<int> orderedQuery;
        subGraph.orderStarNodes(queryGraphInfo.neighbourSign, orderedQuery); // compute query ordering for a star pattern
        matchedSubgraph.orderedNodes = orderedQuery;

        std::vector<int> initialMatches;
        if (!queryGraphInfo.neighbourSign[orderedQuery[0]].empty())
            dataIndex.querySynTrie(queryGraphInfo.neighbourSign[orderedQuery[0]], indexType.synopsesTrie, initialMatches);

        auto it = matchedSubgraph.invalidCands.find(orderedQuery[0]);
        if(it != matchedSubgraph.invalidCands.end())
            initialMatches = findDifference(initialMatches, it->second);

        if(initialMatches.size() < freq)
            return false;

        if(subGraph.findStarMatches(freq, initialMatches, queryGraphInfo, orderedQuery, nodeMatches, indexType.neighborTrie, matchedSubgraph.invalidCands, nodeRepository))
            return true;
        return false;
    }

    /// Compute a set of query sequencing permutations that allows equal distribution of matched nodes for all query nodes. This helps to quickly achieve MNI support.
    Vector2D queryOrderPerms(queryGraphInfo.nNodes);
    for(int i = 0; i < queryGraphInfo.nNodes; ++i)
        subGraph.orderVerticesAllMNI(queryGraphInfo.nNodes, queryGraphInfo.neighbourSign, queryGraphInfo.adjacencyList, queryOrderPerms[i], i);
    matchedSubgraph.orderedNodes = queryOrderPerms[0]; // A pattern needs to a have one permutation of query ordering

    /// Compute the initial candidate solutions for all the query vertices.
    Vector2D initialMatchesAll(queryGraphInfo.nNodes);
    if(!findInitialMatches(freq, indexType, queryGraphInfo.new_old, initialMatchesAll, queryOrderPerms, matchedSubgraph, queryGraphInfo.neighbourSign)){
        return false;
    }

    std::vector<int> imSize(initialMatchesAll.size());
    for(size_t i = 0; i < initialMatchesAll.size(); ++i)
        imSize[i] = initialMatchesAll[i].size();

    std::vector<int> pO = sortIndexIncr(imSize); // order on the choice of permutation - "queryOrderPerms" |

    /// Case 3: If the pattern is a quasi-clique pattern
    bool notFrequent = false;
    for(size_t k = 0; k < queryGraphInfo.nNodes; ++k){
        int i = queryOrderPerms[pO[k]][0];
        if(!nodeRepository[i].empty())
            initialMatchesAll[i] = findDifferenceSet(initialMatchesAll[i], nodeRepository[i]);

        if(subGraph.findMatchesFastMNI(freq, initialMatchesAll[i], queryGraphInfo, queryOrderPerms[i], nodeMatches, indexType.neighborTrie, matchedSubgraph.invalidCands, nodeRepository, notFrequent))
            return true;
        else if(notFrequent)
            return false;
    }
    return false;
}
