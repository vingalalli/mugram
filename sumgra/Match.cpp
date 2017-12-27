/*
 *  Copyright (c) 2015 Vijay Ingalalli
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

#include "Match.h"

Match::Match()
{
    //ctor
}

Match::~Match()
{
    //dtor
}

void collectInvalidCands(InvalidCands& invalidCands, const int& patternNode, const int& matchedNode)
{
    auto it = invalidCands.find(patternNode);
    if(it == invalidCands.end()){
        std::vector<int> tmp {matchedNode};
        invalidCands.insert(std::make_pair(patternNode, tmp));
    }
    else{
        it->second.push_back(matchedNode);
    }
}

bool fetchIsoEmbRules(const int& freq, const int& iVsolution, const Vector2D& satSolutions, NodeRepository&  nodeRepository, Flags& flags)
{
    /// The following rules help us to tap the Cartesian product of the satellite solutions that respect the condition of isomorphic embedding. When the rules fail, the default DFS search is conducted.

    int nNodes = satSolutions.size()+1;
    std::vector<int> satSolutionSize(nNodes-1);
//    VecOfSet satSolnSet(nNodes-1); // this is a little bit expensive |
    std::set<int> s1Matches;
    for(size_t i = 0; i < satSolutions.size(); ++i){
        satSolutionSize[i] = satSolutions[i].size();
        /// Rule 0. If 2 nodes of a pattern have only 1 solution, and they are the same, then we can never have an isomorphic match. Hence its an invalid match.
        if(satSolutionSize[i] == 1){
            int prev = 0;
            if(!s1Matches.empty())
                prev = s1Matches.size();
            s1Matches.insert(satSolutions[i][0]);
            if(prev == s1Matches.size()){
                flags.invalidMatch = true;
                return false;
            }
        }
    }

    /// Rule 1. For every satellite node, the solutions size should be as big as the #satNodes.
    bool minSizeExists = true;
    for(size_t s = 0; s < satSolutions.size(); ++s){
        if(satSolutionSize[s] < nNodes){
            minSizeExists = false;
            break;
        }
    }

    /// Rule 2.
    /// 1, BUILD A TREE INDEX.
    /// 1. If a solution group of size n exits more than n times, then the whole pattern is invalid.
//        bool requiredSizeExists = true;
//        std::vector<int> noOfSatSolSize(nNodes-1, 0); // keep its 0'th element as dummy
//        for(size_t s = 0; s < nNodes; ++s){
//            if(satSolutionSize[s] < nNodes){
//                noOfSatSolSize[satSolutionSize[s]]+=1;
//            }
//        }
//        for(size_t s = 1; s < noOfSatSolSize.size(); ++s){
//            if(noOfSatSolSize[s] <= s && noOfSatSolSize[s] != 0){
//                for(size_t t = 0; t < )
//            }
////            else{
////                // check if they have distinct subsets in them |
////            }
//        }

    if(minSizeExists){
        nodeRepository[0].insert(iVsolution); // core node can be added on its own, as there is at least one isomorphic embedding |
        if(!flags.satBinFilled){
            bool binsFilled = true;
            for(size_t s = 0; s < satSolutions.size(); ++s){
                for(size_t t = 0; t < satSolutions[s].size(); ++t)
                    nodeRepository[s+1].insert(satSolutions[s][t]);
                if(nodeRepository[s+1].size() < freq)
                    binsFilled = false;
            }
            if(binsFilled)
                flags.satBinFilled = true;
        }
        if(flags.satBinFilled && nodeRepository[0].size() >= freq)
            return true;
    }

    if(minSizeExists) // at least one rule holds true;
        flags.validRule = true;
    return false;
}

/// Maintain stack indexes instead of stack variables. (This will be very time efficient.)
bool fetchIsoEmbDFS(const int& freq, const int& iVsolution, const Vector2D& satSolutions, NodeRepository&  nodeRepository, Flags& flags)
{
    int nNodes = satSolutions.size()+1;
    std::vector<int> solution; // the solution is not necessarily isomorphic
    solution.push_back(iVsolution);
    std::unordered_set<int> solutionIso; // tests for isomorphism of the solution |
    solutionIso.insert(iVsolution);

    int depth = 1;
    Vector2D pStack(nNodes-1);
    bool rootMatched = true;
    bool incr = true; // keeps checking if we are moving forward in search tree space |
    while (depth != 0){
        if (solution.size() == nNodes){
            if(solutionIso.size() == nNodes){ // Its an isomorphic embedding
                flags.embFound = true;
                if(rootMatched){ // make sure that "iVsolution" is inserted only once |
                    nodeRepository[0].insert(iVsolution);
                    rootMatched = false;
                }
                if(!flags.satBinFilled){
                    int satSupReached = 0;
                    for(size_t k = 1; k < nNodes; ++k){
                        if(nodeRepository[k].size() < freq)
                            nodeRepository[k].insert(solution[k]);
                        else
                            ++satSupReached;
                    }
                    if(satSupReached == nNodes-1) // MNI of each satellite node bin is filled |
                        flags.satBinFilled = true;
                }
                else{ // all satellite nodes have been matched |
                    if(nodeRepository[0].size() >= freq)
                        return true;
                }
            }

            solutionIso.erase(solution.back());
            solution.pop_back();
            --depth; // this restores the 'depth' value going out of bound |
            pStack[depth-1].pop_back(); // remove the stack element only when entire length is matched |
            incr = false;
        }
        else {
            if (incr){
                pStack[depth-1] = satSolutions[depth-1];
            }
            else{
                solutionIso.erase(solution.back());
                solution.pop_back();
                pStack[depth-1].pop_back();
            }
        }
        if (!pStack[depth-1].empty()) {
            solutionIso.insert(pStack[depth-1].back());
            solution.push_back(pStack[depth-1].back());

            /// Fetch the values until isomorphic value is accessed |
            while(solution.size() != solutionIso.size() && !pStack[depth-1].empty()){
                solution.pop_back(); // remove repeated solution from a possible embedding |
                pStack[depth-1].pop_back(); // remove the repeated solution from the stack |
                if (!pStack[depth-1].empty()){
                    solutionIso.insert(pStack[depth-1].back()); // add the new solution |
                    solution.push_back(pStack[depth-1].back());
                }
            }
            /// A stack level can get empty at any level, because of the while loop before.
            if(pStack[depth-1].empty()){
                --depth;
                incr = false;
                continue;
            }

            if(solution.size() == solutionIso.size()){ // an isomorphic match is found |
                ++depth;
                incr = true;
            }
            else{ // no iso match is found => "pStack[depth-1]" has gone empty |
                --depth;
                incr = false;
            }
        }
        else {
            --depth;
            incr = false;
        }
    }
    return false;
}

bool fetchIsoEmbDFSOrd(const int& freq, const int& initVertex, const Vector2D& satSolutions, const Vector& orderIndex, NodeRepository&  nodeRepository, Flags& flags)
{
    int nNodes = satSolutions.size()+1;
    std::vector<int> solution; // the solution is not necessarily isomorphic
    solution.push_back(initVertex);
    std::unordered_set<int> solutionIso; // tests for isomorphism of the solution |
    solutionIso.insert(initVertex);

    int depth = 1;
    Vector2D pStack(nNodes-1);
    bool rootMatched = true;
    bool incr = true; // keeps checking if we are moving forward in search tree space |
    while (depth != 0){
        if (solution.size() == nNodes){
            if(solutionIso.size() == nNodes){ // Its an isomorphic embedding
                flags.embFound = true;
/// instead of just checking only "rootMatched", make sure to insert all the size-1 solution vertices, only once. This is possinle since our "satSolutions" are ordered in ascending size.
                if(rootMatched){ // make sure that "satSolutions[0][0]" is inserted only once |
                    nodeRepository[0].insert(initVertex);
                    rootMatched = false;
                }
                if(!flags.satBinFilled){
                    int satSupReached = 0;
                    for(size_t k = 1; k < nNodes; ++k){
                        if(nodeRepository[orderIndex[k-1] + 1].size() < freq)
                            nodeRepository[orderIndex[k-1] + 1].insert(solution[k]);
                        else
                            ++satSupReached;
                    }
                    if(satSupReached == nNodes-1) // MNI of each satellite node bin is filled |
                        flags.satBinFilled = true;
                }
                else{ // all satellite nodes have been matched |
                    if(nodeRepository[0].size() >= freq)
                        return true;
                }
            }

            solutionIso.erase(solution.back());
            solution.pop_back();
            --depth; // this restores the 'depth' value going out of bound |
            pStack[depth-1].pop_back(); // remove the stack element only when entire length is matched |
            incr = false;
        }
        else {
            if (incr){
                pStack[depth-1] = satSolutions[depth-1];
            }
            else{
                solutionIso.erase(solution.back());
                solution.pop_back();
                pStack[depth-1].pop_back();
            }
        }
        if (!pStack[depth-1].empty()) {
            solutionIso.insert(pStack[depth-1].back());
            solution.push_back(pStack[depth-1].back());
            while(solution.size() != solutionIso.size() && !pStack[depth-1].empty()){
                solution.pop_back(); // remove repeated solution from a possible embedding |
                pStack[depth-1].pop_back(); // remove the repeated solution from the stack |
                if (!pStack[depth-1].empty()){
                    solutionIso.insert(pStack[depth-1].back()); // add the new solution |
                    solution.push_back(pStack[depth-1].back());
                }
            }
            /// A stack level can get empty at any level, because of the while loop before.
            if(pStack[depth-1].empty()){
                --depth;
                incr = false;
                continue;
            }
            if(solution.size() == solutionIso.size()){ // an isomorphic match is found |
                ++depth;
                incr = true;
            }
            else{ // no iso match is found => "pStack[depth-1]" has gone empty |
                --depth;
                incr = false;
            }
        }
        else {
            --depth;
            incr = false;
        }
    }
    return false;
}

void subgraphSearch(const std::vector<int>& querySequence, const int& p_m, const VecOfSet& nodeMatches, const Vector2D& matchedQueryNeighbours, const std::vector<int>& matchedDataVertices, const std::vector<Trie*>& dataNeighTrie, const EdgeLabelMap& queryEdges, Vector2D& exact_stack)
{
    Index index;
    int nextVertex = querySequence[p_m];
    std::vector<int> matchedDataNeighbours(matchedQueryNeighbours[p_m-1].size());

    for(size_t i = 0; i < matchedDataNeighbours.size(); ++i)
        matchedDataNeighbours[i] = matchedDataVertices[find(querySequence.begin(), querySequence.end(), matchedQueryNeighbours[p_m-1][i]) - querySequence.begin()];

    bool allNbrsMatched = true;
    std::vector<int> edgeMatches;
    for(size_t i = 0; i < matchedQueryNeighbours[p_m-1].size(); ++i) {
        auto query_it = queryEdges.find(std::make_pair(matchedQueryNeighbours[p_m-1][i], nextVertex));
        std::vector<int> nbrLblMatches;
        index.queryNeighTrie(dataNeighTrie[matchedDataNeighbours[i]], query_it->second, nbrLblMatches);
        if (!nbrLblMatches.empty()){
            if (i == 0) {
                edgeMatches = nbrLblMatches;
            }
            else {
                std::unordered_set<int> s(edgeMatches.begin(), edgeMatches.end()); // perform intersection
                edgeMatches.clear();
                for (auto it = nbrLblMatches.begin(); it != nbrLblMatches.end(); ++it) {
                    if (s.find((*it)) != s.end())
                        edgeMatches.emplace_back((*it));
                }
            }
        }
        else{
            allNbrsMatched = false;
            break;
        }
    }

    if (allNbrsMatched) {
        std::vector<int> exactMatches;
        if (!nodeMatches[p_m].empty()) { // Find the common elements between edge labels and node labels |
            for (size_t i = 0; i < edgeMatches.size(); ++i){
                if (nodeMatches[p_m].find(edgeMatches[i]) != nodeMatches[p_m].end())
                    exactMatches.push_back(edgeMatches[i]);
            }
        }
        else
            exactMatches=edgeMatches;

        for(auto it_e = exactMatches.begin(); it_e != exactMatches.end(); ++it_e) {
            auto dataEnd_it = matchedDataVertices.begin()+p_m;
            if (find(matchedDataVertices.begin(), dataEnd_it, (*it_e)) == dataEnd_it)
                    exact_stack[p_m-1].emplace_back((*it_e));
        }
    }
}

void subgraphSearchMNI(const std::vector<int>& querySequence, const int& p_m, const VecOfSet& nodeMatches, const Vector2D& matchedQueryNeighbours, const std::vector<int>& matchedDataVertices, const std::vector<Trie*>& dataNeighTrie, const EdgeLabelMap& queryEdges, Vector2D& exact_stack, std::map<int, std::vector<int>>& invalidCands)
{
    Index index;
    int nextVertex = querySequence[p_m];
    std::vector<int> matchedDataNeighbours(matchedQueryNeighbours[p_m-1].size());

    for(size_t i = 0; i < matchedDataNeighbours.size(); ++i)
        matchedDataNeighbours[i] = matchedDataVertices[find(querySequence.begin(), querySequence.end(), matchedQueryNeighbours[p_m-1][i]) - querySequence.begin()];

    bool allNbrsMatched = true;
    std::vector<int> edgeMatches;
    for(size_t i = 0; i < matchedQueryNeighbours[p_m-1].size(); ++i) {
        auto query_it = queryEdges.find(std::make_pair(matchedQueryNeighbours[p_m-1][i], nextVertex));
        std::vector<int> nbrLblMatches;
        index.queryNeighTrie(dataNeighTrie[matchedDataNeighbours[i]], query_it->second, nbrLblMatches);
        if (!nbrLblMatches.empty()){
            if (i == 0) {
                edgeMatches = nbrLblMatches;
            }
            else {
                std::unordered_set<int> s(edgeMatches.begin(), edgeMatches.end()); // perform intersection
                edgeMatches.clear();
                for (auto it = nbrLblMatches.begin(); it != nbrLblMatches.end(); ++it) {
                    if (s.find((*it)) != s.end())
                        edgeMatches.emplace_back((*it));
                }
            }
        }
        else{
            allNbrsMatched = false;
            break;
        }
    }

    if (allNbrsMatched) {
        std::vector<int> exactMatches;
        if (!nodeMatches[p_m].empty()) { // Find the common elements between edge labels and node labels |
            for (size_t i = 0; i < edgeMatches.size(); ++i){
                if (nodeMatches[p_m].find(edgeMatches[i]) != nodeMatches[p_m].end())
                    exactMatches.push_back(edgeMatches[i]);
            }
        }
        else
            exactMatches=edgeMatches;

/// The below commented code integrates pruning on the fly, rather than computing separately after the 'exactMatches' for loop. However the tests show that this is inefficient owing to the two if conditions in the for loop.

//        auto it = invalidCands.find(querySequence[p_m]);
//        std::unordered_set<int> invalid;
//        if(it != invalidCands.end())
//            std::copy(it->second.begin(), it->second.end(), std::inserter(invalid, invalid.end()));

        for(auto it_e = exactMatches.begin(); it_e != exactMatches.end(); ++it_e) {
            auto dataEnd_it = matchedDataVertices.begin()+p_m;
            if (find(matchedDataVertices.begin(), dataEnd_it, (*it_e)) == dataEnd_it){
//                if(!invalid.empty()){
//                    if(invalid.find(*it_e) == invalid.end())
//                        exact_stack[p_m-1].emplace_back((*it_e));
//                }
//                else
                    exact_stack[p_m-1].push_back((*it_e));

            }
        }
//if(it != invalidCands.end())
//cout << "After:  " << exact_stack[p_m-1].size() << endl;

        /// prune using the propagated invalid cands |

//        auto it = invalidCands.find(querySequence[p_m]);
//        if(it != invalidCands.end()){
//            exact_stack[p_m-1] = findDifference(exact_stack[p_m-1], it->second);
//        }

    }
}

void Match::chooseFrontier(const std::vector<int>& alreadyM, const Vector2D& queryAdjacencyList, std::vector<int>& frontier)
{
    for (size_t i = 0; i < alreadyM.size(); ++i) {
        if (i == 0)
            frontier = queryAdjacencyList[alreadyM[i]];
        else
            frontier = setUnion(frontier, queryAdjacencyList[alreadyM[i]]);
    }
    frontier = findDifference(frontier, alreadyM);
}

void Match::orderStarNodes(const EdgeLabel& neighbourSign, Vector& orderedQuery)
{
    std::vector<int> n_of_edges;
    for (size_t i = 0; i < neighbourSign.size(); ++i){
        int allDim = 0;
        for (int j = 0; j < neighbourSign[i].size(); ++j)
            allDim += neighbourSign[i][j].size();
        n_of_edges.push_back(allDim);
    }
    orderedQuery = sortIndexDecr(n_of_edges); //; Decreasing order
}

void Match::orderVertices(const int& queryNodes, const EdgeLabel& queryNeighbourSign, const Vector2D& queryAdjacencyList, std::vector<int>& querySequence)
{
    std::vector<int> signRank;
    for (size_t i = 0; i < queryNodes; ++i) {
        int allDim = 0;
        for (int j = 0; j < queryNeighbourSign[i].size(); ++j)
            allDim += queryNeighbourSign[i][j].size();
        signRank.push_back(allDim);
    }


/// Using only neighbourhood signature size |
    querySequence.push_back(sortIndexDecr(signRank)[0]); //; Decreasing order of signRank
    while (querySequence.size() != queryNodes) {
        std::vector<int> frontier;
        chooseFrontier(querySequence, queryAdjacencyList, frontier);
        if (frontier.size() == 1)
            querySequence.push_back(frontier[0]);
        else {
            std::vector<int> rank_1(frontier.size());

            for (size_t i = 0; i < frontier.size(); ++i)
                rank_1[i] = 0;
            std::vector<int> srt_rank_1 = sortIndexDecr(rank_1);

/// Check for 2 constraints;
//
            if (rank_1[srt_rank_1[0]] != rank_1[srt_rank_1[1]])
                querySequence.push_back(frontier[srt_rank_1[0]]);
            else {
                std::vector<int> same_rank_1;
                same_rank_1.push_back(frontier[srt_rank_1[0]]);
                int i = 0;
                while (i < frontier.size()-1){
                    if (rank_1[srt_rank_1[i]] == rank_1[srt_rank_1[i+1]])
                        same_rank_1.push_back(frontier[srt_rank_1[i+1]]);
                    else
                        break;
                    ++i;
                }
                std::vector<int> rank_2(same_rank_1.size());
                for (size_t i = 0; i < same_rank_1.size(); ++i)
                    rank_2[i] = signRank[same_rank_1[i]];
                std::vector<int> srt_rank_2 = sortIndexDecr(rank_2);
                querySequence.push_back(same_rank_1[srt_rank_2[0]]);
            }
//
        }
    }
    int initialVertex = querySequence[0];
}


void Match::orderVerticesAllMNI(const int& queryNodes, const EdgeLabel& queryNeighbourSign, const Vector2D& queryAdjacencyList, std::vector<int>& orderedNodes, int& initialVertex)
{
//    std::vector<int> querySequence;

    std::vector<int> signRank;
    for (size_t i = 0; i < queryNodes; ++i) {
        int allDim = 0;
        for (int j = 0; j < queryNeighbourSign[i].size(); ++j)
            allDim += queryNeighbourSign[i][j].size();
        signRank.push_back(allDim);
    }

/// Using only neighbourhood signature size |
//    orderedNodes.push_back(sortIndexDecr(signRank)[0]); //; Decreasing order of
    orderedNodes.push_back(initialVertex); //; Decreasing order of


    while (orderedNodes.size() != queryNodes) {
        std::vector<int> frontier;
        chooseFrontier(orderedNodes, queryAdjacencyList, frontier);
        if (frontier.size() == 1)
            orderedNodes.push_back(frontier[0]);
        else {
            std::vector<int> rank_1(frontier.size());

            for (size_t i = 0; i < frontier.size(); ++i)
                rank_1[i] = 0;
            std::vector<int> srt_rank_1 = sortIndexDecr(rank_1);
/// Check for 2 constraints;
//
            if (rank_1[srt_rank_1[0]] != rank_1[srt_rank_1[1]])
                orderedNodes.push_back(frontier[srt_rank_1[0]]);
            else {
                std::vector<int> same_rank_1;
                same_rank_1.push_back(frontier[srt_rank_1[0]]);
                int i = 0;
                while (i < frontier.size()-1){
                    if (rank_1[srt_rank_1[i]] == rank_1[srt_rank_1[i+1]])
                        same_rank_1.push_back(frontier[srt_rank_1[i+1]]);
                    else
                        break;
                    ++i;
                }
                std::vector<int> rank_2(same_rank_1.size());
                for (size_t i = 0; i < same_rank_1.size(); ++i)
                    rank_2[i] = signRank[same_rank_1[i]];
                std::vector<int> srt_rank_2 = sortIndexDecr(rank_2);
                orderedNodes.push_back(same_rank_1[srt_rank_2[0]]);
            }
//
        }
    }
}


/// Check if a subgraph is found
bool Match::checkSI(std::vector<int>& initialMatches, const GraphParameter& queryGraphInfo, const std::vector<int>& querySequence, std::vector<int>& matchedDataVertices, const VecOfSet& nodeMatches, const std::vector<Trie*>& dataNeighTrie)
{

    bool found = false;
    Vector2D matchedQueryNeighbours(queryGraphInfo.nNodes-1);
    std::vector<int> q_v;
    q_v.push_back(querySequence[0]);
    for(size_t i = 0; i < queryGraphInfo.nNodes-1; ++i) {
        int nextVertex = querySequence[i+1];
        matchedQueryNeighbours[i] = setIntersection(queryGraphInfo.adjacencyList[nextVertex], q_v);
        q_v.push_back(nextVertex);
    }


    for (auto it = initialMatches.begin(); it != initialMatches.end(); ++it) {
        int p_m = 0;
        matchedDataVertices[p_m] = *it;
        ++p_m;
        Vector2D exactStack(queryGraphInfo.nNodes-1);
        bool incr = true; // keeps checking if we are moving forward in search tree space |
        while (p_m != 0 ) {
            if (p_m == querySequence.size()) {
                found = true;
                return true;

                --p_m;
                incr = false;
            }
            else {
                if (incr)
                    subgraphSearch(querySequence, p_m, nodeMatches, matchedQueryNeighbours, matchedDataVertices, dataNeighTrie, queryGraphInfo.eLabelMap, exactStack);
            }

            if (!exactStack[p_m-1].empty()) {
                matchedDataVertices[p_m] = exactStack[p_m-1].back();
                exactStack[p_m-1].pop_back();
                ++p_m;
                incr = true;
            }
            else {
                --p_m;
                incr = false;
            }
        }
        if(found)
            return true;
        /// Collect the invalid candidates for the initial vertex.
        return found;
    }
}


bool Match::findStarMatches(const int& freq, const Vector& initialMatches, GraphParameter& queryGraphInfo, const Vector& orderedQuery, const VecOfSet& nodeMatches, const std::vector<Trie*>& dataNeighTrie, InvalidCands& invalidCands, NodeRepository&  nodeRepository)
{
    Index index;
    int dynamicOffset = initialMatches.size(); // decrements for every invalid root node |
    int nNodes = orderedQuery.size();
    Flags flags;

    for (size_t v = 0; v < initialMatches.size(); ++v){
        if(dynamicOffset < freq)
            return false;

        /// Collect the solutions for satellite nodes.
        Vector2D satSolutions(nNodes-1);
        int nM = 1; // This counts the initialVertex |
        std::map<std::set<int>, int> repeatedMultidges;
        for(; nM < nNodes; ++nM){
            std::pair<int, int> edge = std::make_pair(orderedQuery[0], orderedQuery[nM]);
            std::set<int> multiedge = queryGraphInfo.eLabelMap.find(edge)->second;
            auto it = repeatedMultidges.find(multiedge);
            if(it != repeatedMultidges.end()){
                satSolutions[nM-1] = satSolutions[it->second];
                continue;
            }
            else
                repeatedMultidges.insert(std::make_pair(multiedge, nM-1));
            index.queryNeighTrie(dataNeighTrie[initialMatches[v]], multiedge, satSolutions[nM-1]);
            if (satSolutions[nM-1].empty()) // Star pattern structure not found
                break;
        }

        /// No embeddings found; root solution is invalid.
        if (nM < nNodes){
            --dynamicOffset;
            collectInvalidCands(invalidCands, orderedQuery[0], initialMatches[v]);
            continue;
        }

        flags.invalidMatch = false;
        flags.embFound = false;
        flags.validRule = false;

        /// Priority method: Use a set of rules to find isomorphic solutions from the satellite matches.

        if(fetchIsoEmbRules(freq, initialMatches[v], satSolutions, nodeRepository, flags))
            return true;
        else{
            if(flags.validRule){ // a valid rule exists, and hence skip DFS approach |
                continue;
            }
            else if(flags.invalidMatch){
                --dynamicOffset;
                collectInvalidCands(invalidCands, orderedQuery[0], initialMatches[v]);
                continue;
            }
        }

        /// Default method: DFS approach to find isomorphic solutions from the satellite matches.
        /// Order the solutions in the increasing order |
//        Vector solSize(satSolutions.size());
//        for(size_t i = 0; i < satSolutions.size(); ++i)
//            solSize[i] = satSolutions[i].size();
//        Vector orderIndex = SortIndexIncr(solSize);
//        Vector2D satSolutionsOrd(satSolutions.size());
//        for(size_t i = 0; i < satSolutions.size(); ++i)
//            satSolutionsOrd[i] = satSolutions[orderIndex[i]];
//
//
//        if(fetchIsoEmbDFSOrd(freq, initialMatches[v], satSolutionsOrd, orderIndex, nodeRepository, flags)){
        if(fetchIsoEmbDFS(freq, initialMatches[v], satSolutions, nodeRepository, flags)){
            return true;
        }


        /// No embeddings found; root solution is invalid.
        if(!flags.embFound){
            --dynamicOffset;
            collectInvalidCands(invalidCands, orderedQuery[0], initialMatches[v]);
        }

    }
    return false;
}

bool Match::findMatchesFastMNI(const int& freq, const std::vector<int>& initialMatches, const GraphParameter& queryGraphInfo, const std::vector<int>& querySequence, const VecOfSet& nodeMatches, const std::vector<Trie*>& dataNeighTrie, std::map<int, std::vector<int>>& invalidCands, std::vector<std::unordered_set<int>>&  nodeRepository, bool& notFrequent)
{
//    bool supportReached = false;
    Vector2D matchedQueryNeighbours(queryGraphInfo.nNodes-1);
    std::vector<int> q_v;
    q_v.push_back(querySequence[0]);
    for(size_t i = 0; i < queryGraphInfo.nNodes-1; ++i) {
        int nextVertex = querySequence[i+1];
        matchedQueryNeighbours[i] = setIntersection(queryGraphInfo.adjacencyList[nextVertex], q_v);
        q_v.push_back(nextVertex);
    }

    std::vector<int> matchedDataVertices(queryGraphInfo.nNodes);
    clock_t stop_time = clock();

    /// Iterative Approach
    int updatedMatchSize = initialMatches.size() + nodeRepository[querySequence[0]].size(); // potential candidates + already matched vertices |

    for (auto it = initialMatches.begin(); it != initialMatches.end(); ++it) {
        if(updatedMatchSize < freq){ // quit if "initialMatches" shrinks smaller than support |
            notFrequent = true;
            return false;
        }
        bool initNotMatched = true;
        int p_m = 0;
        matchedDataVertices[p_m] = (*it);
        ++p_m;
        Vector2D exactStack(queryGraphInfo.nNodes-1);
        bool incr = true; // keeps checking if we are moving forward in search tree space |
        while (p_m != 0 ) {
            if (p_m == querySequence.size()) {
                initNotMatched = false;
/// MNI tested during subgraph search.
// /*
                bool MniReached = true;
                for(size_t r = 0; r < matchedDataVertices.size(); ++r){
                    if (nodeRepository[querySequence[r]].size() < freq)
                        nodeRepository[querySequence[r]].insert(matchedDataVertices[r]);
                    if(nodeRepository[querySequence[r]].size() < freq)
                        MniReached = false;
                }
// */
                if(MniReached)
                    return true;
                break; // one match is found. escape while loop |
            }
            else {
                if (incr)
                    subgraphSearchMNI(querySequence, p_m, nodeMatches, matchedQueryNeighbours, matchedDataVertices, dataNeighTrie, queryGraphInfo.eLabelMap, exactStack, invalidCands);
            }

            if (!exactStack[p_m-1].empty()) {
                matchedDataVertices[p_m] = exactStack[p_m-1].back();
                exactStack[p_m-1].pop_back();
                ++p_m;
                incr = true;
            }
            else {
                --p_m;
                incr = false;
            }
        }
//        if(supportReached)
//            return true;
        /// Collect the invalid candidates for the initial vertex.
        if(initNotMatched){
            --updatedMatchSize;
            auto it1 = invalidCands.find(querySequence[0]);
                if(it1 == invalidCands.end()){
                    std::vector<int> tmp {*it};
                    invalidCands.insert(std::make_pair(querySequence[0], tmp));
                }
                else
                    it1->second.push_back(*it);
        }
    }
    return false;
}
