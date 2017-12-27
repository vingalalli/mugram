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

#include "GraphIso.h"
#include "SubgraphIso.h"
#include "File.h"

GraphIso::GraphIso()
{
    //ctor
}

GraphIso::~GraphIso()
{
    //dtor
}


/// Find permutations without vertex labels.
void GraphIso::computePermutationGroups(const GraphParameter& pattern, const IsoTestMap& dataGraph, Vector2D& permutationGroups)
{
    std::map<int, int> vertexDegree;
    for(size_t i = 0; i < pattern.adjacencyList.size(); ++i)
        vertexDegree.insert(std::make_pair(i, pattern.adjacencyList[i].size()));

    Vector2D vtxDegreeGroup;
    partitionNodes(vertexDegree, vtxDegreeGroup); // Nodes partitioned according to vertex degrees

    for(size_t i = 0; i < vtxDegreeGroup.size(); ++i){
        if(vtxDegreeGroup[i].size() == 1) // No further partition required
            permutationGroups.push_back(vtxDegreeGroup[i]);
        else{
            /// Fetch the neighborhood information for vertex partitioning.
            std::map<int, std::set<std::vector<int>>> edgeLabels;
            for(size_t s = 0; s < vtxDegreeGroup[i].size(); ++s){ // For each vertex with the same partition |
                std::set<std::vector<int>> adjacentTokens; // An ordered collection of neighborhood information. Order is necessary for correct partition |
                std::vector<int> adjListTemp = pattern.adjacencyList[vtxDegreeGroup[i][s]];
                for(size_t t = 0; t < adjListTemp.size(); ++t){
                    std::vector<int> neighbourInfo(2); // IF vertex labels exist, then this will be of size = 3 |
                    neighbourInfo[0] = dataGraph.multiedgeMap.find(pattern.eLabelMap.find(std::make_pair(vtxDegreeGroup[i][s], adjListTemp[t]))->second)->second; // Get the int-tokenized multiedge between "vtxDegreeGroup[i][s]" and "ADJ(vtxDegreeGroup[i]) = adjListTemp[t]" |
                    neighbourInfo[1] = pattern.adjacencyList[adjListTemp[t]].size(); // Get the degree of 't'th adjacent vertex of vtxDegreeGroup[i][s] |
                    adjacentTokens.insert(neighbourInfo);
                }
                edgeLabels.insert(std::make_pair(vtxDegreeGroup[i][s], adjacentTokens));
            }

            Vector2D edgeLabelGroup;
            partitionNodes(edgeLabels, edgeLabelGroup); // Nodes partitioned according to similarity among the set of multiedges connected to the corresponding vertices
            for(size_t k = 0; k < edgeLabelGroup.size(); ++k)
                permutationGroups.push_back(edgeLabelGroup[k]); // This is the final partition; in case of further ideas, use if/else again
        }
    }
}

void GraphIso::computeVertexPermutation(Vector2D& permutationGroups, Vector2D& vertexPermutations)
{
    Vector3D allPermutationGroups(permutationGroups.size());
    for(size_t i = 0; i < permutationGroups.size(); ++i){
        do{
            allPermutationGroups[i].push_back(permutationGroups[i]);
        }
        while(std::next_permutation(permutationGroups[i].begin(),permutationGroups[i].end()));
    }


    /// Collect all the vertex permutations (this can also be implemented using a recursion instead of iterative approach)
    for(size_t i = 0; i < allPermutationGroups[0].size(); ++i){
        Vector2D singlePermutation;
        singlePermutation.push_back(allPermutationGroups[0][i]); /// Add elements to this vector while Incrementing and remove while Decrementing.
        if(allPermutationGroups.size() == 1){ // If only 1 group exists |
            vertexPermutations.push_back(singlePermutation.back());
            continue;
        }
        // Recursively find all the permutations |
        int depth = 1;
        Vector3D pStack(allPermutationGroups.size()-1);
        bool incr = true; // keeps checking if we are moving forward in search tree space |
        while (depth != 0 ) {
            if (singlePermutation.size() == allPermutationGroups.size()) {
                std::vector<int> nodeSequence; // = singlePermutation[0];
                for(size_t k = 0; k < singlePermutation.size(); ++k)
                    nodeSequence.insert(nodeSequence.end(), singlePermutation[k].begin(), singlePermutation[k].end());

                vertexPermutations.push_back(nodeSequence);
                singlePermutation.pop_back();
                --depth; // this restores the 'depth' value going out of bound |
                pStack[depth-1].pop_back(); // remove the stack element only when entire length is matched |
                incr = false;
            }
            else {
                if (incr)
                    pStack[depth-1] = allPermutationGroups[depth];
                else{
                    singlePermutation.pop_back();
                    pStack[depth-1].pop_back();
                }
            }
            if (!pStack[depth-1].empty()) {
                singlePermutation.push_back(pStack[depth-1].back());
                ++depth;
                incr = true;
            }
            else {
                --depth;
                incr = false;
            }
        }
    }
}

inline bool hasBigPermGroup(const int& permSize, const Vector2D& permutationGroups, int& bigGroupPos)
{
    for(size_t i = 0; i < permutationGroups.size(); ++i){
        if(permutationGroups[i].size() > permSize){
            bigGroupPos = i;
            return true;
        }
    }
    return false;
}

/// Find permutations without vertex labels. FOR induced permutation.
void computePermutationGroupsNew(const std::map<int, std::vector<int>>& adjacencyList, const EdgeLabelMap&  eLabelMap, const IsoTestMap& dataGraph, Vector2D& permutationGroups)
{
    std::map<int, int> vertexDegree;
    for(auto it = adjacencyList.begin(); it != adjacencyList.end(); ++it)
        vertexDegree.insert(std::make_pair(it->first, it->second.size()));
//    for(size_t i = 0; i < pattern.adjacencyList.size(); ++i)
//        vertexDegree.insert(std::make_pair(i, pattern.adjacencyList[i].size()));

    Vector2D vtxDegreeGroup;
    partitionNodes(vertexDegree, vtxDegreeGroup); // Nodes partitioned according to vertex degrees

    for(size_t i = 0; i < vtxDegreeGroup.size(); ++i){
        if(vtxDegreeGroup[i].size() == 1) // No further partition required
            permutationGroups.push_back(vtxDegreeGroup[i]);
        else{
            /// Fetch the neighborhood information for vertex partitioning.
            std::map<int, std::set<std::vector<int>>> edgeLabels;
            for(size_t s = 0; s < vtxDegreeGroup[i].size(); ++s){ // For each vertex with the same partition |
                std::set<std::vector<int>> adjacentTokens; // An ordered collection of neighborhood information. Order is necessary for correct partition |
//                std::vector<int> adjListTemp = pattern.adjacencyList[vtxDegreeGroup[i][s]];
                auto it_s = adjacencyList.find(vtxDegreeGroup[i][s]);
                std::vector<int> adjListTemp = it_s->second;

                for(size_t t = 0; t < adjListTemp.size(); ++t){
                    std::vector<int> neighbourInfo(2); // IF vertex labels exist, then this will be of size = 3 |
                    neighbourInfo[0] = dataGraph.multiedgeMap.find(eLabelMap.find(std::make_pair(vtxDegreeGroup[i][s], adjListTemp[t]))->second)->second; // Get the int-tokenized multiedge between "vtxDegreeGroup[i][s]" and "ADJ(vtxDegreeGroup[i]) = adjListTemp[t]" |
//                    neighbourInfo[1] = pattern.adjacencyList[adjListTemp[t]].size(); // Get the degree of 't'th adjacent vertex of vtxDegreeGroup[i][s] |
                    auto it_t = adjacencyList.find(adjListTemp[t]);
                    neighbourInfo[1] = it_t->second.size();
                    adjacentTokens.insert(neighbourInfo);
                }
                edgeLabels.insert(std::make_pair(vtxDegreeGroup[i][s], adjacentTokens));
            }

            Vector2D edgeLabelGroup;
            partitionNodes(edgeLabels, edgeLabelGroup); // Nodes partitioned according to similarity among the set of multiedges connected to the corresponding vertices
            for(size_t k = 0; k < edgeLabelGroup.size(); ++k)
                permutationGroups.push_back(edgeLabelGroup[k]); // This is the final partition; in case of further ideas, use if/else again
        }
    }
}

void inducePartition(const IsoTestMap& dataGraph, const Vector2D& permutationGroups, int& pivot, const GraphParameter& subgraph, std::deque<Vector2D>& bigPermGroups)
{
    std::vector<int> partitionableGroup = permutationGroups[pivot];
    int groupSize = partitionableGroup.size();
    Vector2D newGroup;
    for(size_t i = 0; i < pivot; ++i) // Fetch the groups BEFORE the pivot |
        newGroup.push_back(permutationGroups[i]);

    for(int i = 0; i < groupSize; ++i){ // Induce partition for each element in the BIG group |
        // Each element forms a singleton group - "pivotGroup" |
        std::vector<int> pivotGroup(1);
        pivotGroup[0] = partitionableGroup[i];
        newGroup.push_back(pivotGroup);
        Vector2D oldAdj = subgraph.adjacencyList;
        /// Fetch the adjacency list after a pivot node 'i' is chosen.
        std::map<int, std::vector<int>> newAdj;
        for(size_t j = 0; j < groupSize; ++j){
            if(partitionableGroup[i] != partitionableGroup[j]){
                std::vector<int> aList = oldAdj[partitionableGroup[j]];
                aList.erase(std::remove(aList.begin(), aList.end(), i), aList.end());
                newAdj.insert(std::make_pair(partitionableGroup[j], aList));
            }
        }
        computePermutationGroupsNew(newAdj, subgraph.eLabelMap, dataGraph, newGroup);
        for(size_t j = pivot+1; j < permutationGroups.size(); ++j) // Fetch the groups AFTER the pivot |
            newGroup.push_back(permutationGroups[j]);
        bigPermGroups.push_back(newGroup);
    }
}

void GraphIso::generateCanonical(const GraphParameter& subgraph, std::vector<int>& canonicalForm, const IsoTestMap& dataGraph)
{
    Vector2D permutationGroups, vertexPermutations;
    computePermutationGroups(subgraph, dataGraph, permutationGroups);

    /// Process the BIG permutation groups into smaller ones.
    std::deque<Vector2D> reducedPermutationGroups;
    int bigGroupPos, permSize = 6;
    if(hasBigPermGroup(permSize, permutationGroups, bigGroupPos) && false ){ // Check for the initial set of permutation groups |
        std::deque<Vector2D> bigPermGroups;
        inducePartition(dataGraph, permutationGroups, bigGroupPos, subgraph, bigPermGroups); // Initial partition |
        while(!bigPermGroups.empty()){ // all groups have been classified
            auto it = bigPermGroups.begin();
            if(hasBigPermGroup(permSize, *it, bigGroupPos)){
                std::deque<Vector2D> inducedPermGroups;
                inducePartition(dataGraph, *it, bigGroupPos, subgraph, inducedPermGroups); // Subsequent partitions |
                bigPermGroups.erase(it); //erase the previous set;
                for(auto itB = inducedPermGroups.begin(); itB != inducedPermGroups.end(); ++itB)
                    bigPermGroups.push_back(*itB);
            }
            else{
                reducedPermutationGroups.push_back(*it);
                bigPermGroups.erase(it);
            }
        }
    }
    else
        reducedPermutationGroups.push_back(permutationGroups);

    for(auto it = reducedPermutationGroups.begin(); it != reducedPermutationGroups.end(); ++it)
        computeVertexPermutation(*it, vertexPermutations);

    /// Method to compute a unique canonical representation.
    std::map<std::pair<int, int>, int> edgeAttMapped; // <edge, mapped distinct integer>
    for(auto it = subgraph.eLabelMap.begin(); it != subgraph.eLabelMap.end(); ++it){
        auto itE = dataGraph.multiedgeMap.find(it->second);
        edgeAttMapped.insert(std::make_pair(it->first, itE->second));
    }
    canonicalForm.push_back(0); // initialize the value to minimum
    for(size_t i = 0; i < vertexPermutations.size(); ++i){
        std::vector<int> canonicalTemp;
        // Perform vertical scan |
        for(size_t n = 0; n < vertexPermutations[i].size(); ++n){ // Get the upper triangular entries of edges in the adjacency matrix
            bool upperTriangular = true;
            int m = 0;
            while(upperTriangular && m < vertexPermutations[i].size()){
                if (n > m){ // This step has a complexity of factorial[vertexPermutations[i].size()] |
                    std::pair<int,int> edge = std::make_pair(vertexPermutations[i][m], vertexPermutations[i][n]);
                    auto it = edgeAttMapped.find(edge);
                    if(it != edgeAttMapped.end())
                        canonicalTemp.push_back(it->second); // Edge exists |
                    else
                        canonicalTemp.push_back(0); // NO edge exists |
                    ++m;
                }
                else
                    upperTriangular = false;
            }
        }
        if (canonicalTemp > canonicalForm) // Collect the maximal canonical form
            canonicalForm = canonicalTemp;
    }
}
