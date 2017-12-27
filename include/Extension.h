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

#ifndef EXTENSION_H
#define EXTENSION_H

#include "File.h"
#include "Frequency.h"
#include "SubgraphIso.h"
#include "GraphIso.h"
#include <../sumgra/Trie.h>

typedef std::map<int, std::set<std::vector<int>>> IsoMeasure;

class Extension
{
    public:
        Extension();
        virtual ~Extension();
        void subgraphExtension(const IsoTestMap& dataGraph, Subgraph& oldFSG, const std::set<int>& newEdge, IsoMeasure& canonicalAllLevel, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::deque<Subgraph>& newFSGs, std::deque<Subgraph> oldFSGs, std::deque<IndexType>& maximalFSG, FSGs& fsgs);
        void extendWithNode(const IsoTestMap& dataGraph, Subgraph& oldFSG, const set<int>& newEdge, IsoMeasure& canonicalAllLevel, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::deque<Subgraph>& newFSGs, std::deque<Subgraph> oldFSGs, std::deque<IndexType>& maximalFSG, FSGs& fsgs);
        void extendWithoutNode(const IsoTestMap& dataGraph, Subgraph& oldFSG, const set<int>& newEdge, IsoMeasure& canonicalAllLevel, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::deque<Subgraph>& newFSGs, std::deque<Subgraph> oldFSGs, std::deque<IndexType>& maximalFSG, FSGs& fsgs);
        void extendWithNodeNew(const IsoTestMap& dataGraph, Subgraph& oldFSG, const set<int>& newEdge, IsoMeasure& canonicalAllLevel, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::deque<Subgraph>& newFSGs, std::deque<Subgraph> oldFSGs, std::deque<IndexType>& maximalFSG, FSGs& fsgs);
        void extendWithoutNodeNew(const IsoTestMap& dataGraph, Subgraph& oldFSG, const set<int>& newEdge, IsoMeasure& canonicalAllLevel, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::deque<Subgraph>& newFSGs, std::deque<Subgraph> oldFSGs, std::deque<IndexType>& maximalFSG, FSGs& fsgs);
        void getInvalidCands(std::deque<Subgraph>& oldFSGs, FSGs& fsgs, Subgraph& newSubgraph, IndexType& newSubgraphInd);
        bool isIsomorphicCurrent(Subgraph& newSubgraph, std::set<std::vector<int>>& extraNodeCanonical, const IsoTestMap& dataGraph);
        bool isIsomorphicPrevious(const int& subgraphSize, const std::vector<int>& canonical, IsoMeasure& canonicalAllLevel);
        //bool isFrequent(const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, Subgraph& newSubgraph);
        void extractInfo(const Subgraph& subgraph, GraphParameter& subgraphInfo);
    protected:
    private:
};

#endif // EXTENSION_H
