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

#include "Index.h"

Index::Index()
{
    //ctor
}

Index::~Index()
{
    //dtor
}


void Index::buildAttHash(const VecOfSet& attSign, AttMap& attributeHash)  // A vertex with no attribute should have '-1' assigned to it in the text file |
{
    for (size_t i = 0; i < attSign.size(); ++i) {
        if ( !(attSign[i].size() == 1 && (*attSign[i].begin()) == -1) ) { // map only vertexes that have labels |
            for (auto it_a = attSign[i].begin(); it_a != attSign[i].end(); ++it_a) {
                auto it_m = attributeHash.find((*it_a));
                if (it_m == attributeHash.end()) {
                    std::set<int> attTmp;
                    attTmp.insert(i);
                    attributeHash.insert(make_pair((*it_a), attTmp));
                }
                else
                    it_m->second.insert(i);
            }
        }
    }
}

void Index::queryAttHash(const VecOfSet& queryAtt, const AttMap& attributeHash, VecOfSet&  attMatches)
{
    for (size_t i = 0; i < queryAtt.size(); ++i) {
        if ( !(queryAtt[i].size() == 1 && (*queryAtt[i].begin()) == -1) ) { // only if the vertex has labels |
            size_t k = 0;
            for (auto it_q = queryAtt[i].begin(); it_q != queryAtt[i].end(); ++it_q) {
                if (k == 0)
                    attMatches[i] = attributeHash.find((*it_q))->second;
                else {
//                    attMatches[i] = setIntersection(attMatches[i], attributeHash.find((*it_q))->second);
                    std::set<int> intersect;
                    set_intersection(attMatches[i].begin(),attMatches[i].end(),attributeHash.find((*it_q))->second.begin(),attributeHash.find((*it_q))->second.end(), std::inserter(intersect,intersect.begin()));
                    attMatches[i] = intersect;
                }
                ++k;
            }
        }
    }
}

void Index::createSynopses(const std::vector<std::set<int>>& signature, std::vector<short>& synopses)
{
    std::set<int> uniqueDim;
    int allDim = 0;
    std::vector<int> itemSize;
    for (int j=0; j<signature.size(); ++j) {
        for(auto it = signature[j].begin(); it != signature[j].end(); ++it)
            uniqueDim.insert((*it));
        allDim += signature[j].size();
        itemSize.push_back(signature[j].size());
    }
//    Add "1" to each synopses field since, the RTree does not consider "zero" values; it does not even work with negative values

    synopses[0] = signature.size()+1;
    synopses[1] = uniqueDim.size()+1;
    synopses[2] = allDim+1;
    synopses[3] = *uniqueDim.begin()+1;
    synopses[4] = *uniqueDim.rbegin()+1;
    synopses[5] = *std::max_element(std::begin(itemSize), std::end(itemSize))+1;
}

BoundingBox Index::bounds(std::vector<short> synopses)
{
  	BoundingBox bb;
        for (size_t i = 0; i < synopses.size(); ++i) {
            bb.edges[i].first  = 0;
            bb.edges[i].second = synopses[i];
        }
  	return bb;
}

void Index::sortSignature(EdgeLabel& neighbourSign, std::vector<int>& sortedNodes, const int& elements)
{
    std::vector<int> nEdges(elements);
    for(size_t i = 0; i < elements; ++i) {
        nEdges[i] = 0;
        std::vector<std::set<int>> sortedSign;
        std::vector<int> nNeighbours;
        for(auto it = neighbourSign[i].begin(); it!=neighbourSign[i].end(); ++it) {
            nEdges[i] += (*it).size();
            nNeighbours.push_back((*it).size());
        }
        std::vector<int> sortedNeighIndex = sortIndexDecr(nNeighbours); // sort each signature within itself wrt subsignature size |
        for(size_t j = 0; j < nNeighbours.size(); ++j)
            sortedSign.push_back(neighbourSign[i][sortedNeighIndex[j]]);
//        neighbourSign[i] = sortedSign;
    }
    sortedNodes = sortIndexDecr(nEdges); // sort all the data vertices wrt data adjaceny list of each vertex |
}

void Index::buildSynTrie(const EdgeLabel& dataSignature, const int& dataNodes, RTree& synopsesTrie)
{
    std::vector<short> zero(SYN_SIZE);
    int i =0;
    for(auto it = dataSignature.begin(); it != dataSignature.end(); ++it) {
        if((*it).empty())
            synopsesTrie.Insert(i, bounds(zero)); // if the node has no edges
        else {
            std::vector<short> synopses(SYN_SIZE);
            createSynopses((*it), synopses);
            synopsesTrie.Insert(i, bounds(synopses));
        }
        ++i;
    }
}


void Index::buildNeighTrie(const Vector2D& adjacencyList, const EdgeLabelMap& eLabelMap, const int& dataNodes, std::vector<Trie*>& nbrIndex)
{
    for(int m = 0; m < dataNodes; ++m) {
        Trie *trieSignature = new Trie();
        for(int n = 0; n < adjacencyList[m].size(); ++n) {
            auto it = eLabelMap.find(std::make_pair(m, adjacencyList[m][n]));
            std::vector<int> multi_edge;
            std::copy(it->second.begin(), it->second.end(), std::back_inserter(multi_edge));
            trieSignature->addSignatureDim(multi_edge, adjacencyList[m][n]); //Add the multi-edge and its corresponding vertex id to the trie;
        }
        trieSignature->updateHashTable(trieSignature->Root(), trieSignature->LabelMap);
        std::map<int, Node*> myMapInt;
        trieSignature->assignParent(trieSignature->Root());
        nbrIndex[m] = trieSignature;
    }
}


void Index::querySynTrie(const std::vector<std::set<int>>& initSignature, RTree& synopsesTrie, std::vector<int>& initialMatches)
{
    std::vector<short> qSynopses(SYN_SIZE);
    createSynopses(initSignature, qSynopses);
    BoundingBox bound = bounds(qSynopses);
    Visitor x = synopsesTrie.Query(RTree::AcceptEnclosing(bound), Visitor());
    if (!x.edgeIndices.empty())
        initialMatches = x.edgeIndices;
}


void Index::queryNeighTrie(Trie* t, const std::set<int>& multi_e, std::vector<int>& MatchedIds)
{
    /// Check if all the elements of multiedges are found in this tree
    for(auto it = multi_e.begin(); it != multi_e.end(); ++it)
        if (t->LabelMap.find(*it+1) == t->LabelMap.end())
            return; // if any edge is not present in the tree, quit |

    auto it = multi_e.rbegin();
    std::vector<Node *> matches = t->LabelMap.find(*it+1)->second; // get all the pointers for the end character of multiedge |

    for(size_t i = 0; i < matches.size(); ++i){
        Node* currentNode = matches[i]->nodeParent();
        if(multi_e.size() > 1) {
            int mtch = 0;
            auto rit = multi_e.rbegin();
            std::advance (rit,1); // skip the matching of the last element in 'multi_e' as it is already done |
            for (; rit != multi_e.rend(); rit++) {
                bool found = false;
                while(currentNode->contentInt() != 0) { // Trace back until root is reached |
                    if(currentNode->contentInt() == (*rit)+1){
                        currentNode = currentNode->nodeParent();
                        ++mtch;
                        found = true;
                        break;
                    }
                    else
                        currentNode = currentNode->nodeParent();
                }
                if(!found)
                    break;
            }
            if(mtch == multi_e.size()-1) {
                std::vector<int> tm(matches[i]->vertexId());
                MatchedIds.insert(MatchedIds.end(), tm.begin(), tm.end());
            }
        }
        else {
            std::vector<int> tm(matches[i]->vertexId());
            MatchedIds.insert(MatchedIds.end(), tm.begin(), tm.end());
        }
    }
}
