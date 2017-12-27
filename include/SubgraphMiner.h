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

#ifndef SUBGRAPHMINER_H
#define SUBGRAPHMINER_H

#include "File.h"
#include "SubgraphIso.h"
#include "Extension.h"


class SubgraphMiner
{
    public:
        SubgraphMiner();
        virtual ~SubgraphMiner();
        void collectSubgraphSeeds(const GraphParameter& dataGraphInfo, std::deque<Subgraph>& subgraphSeeds);
        void findFrequentSeeds(const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, std::deque<Subgraph>& newSubgraphs, std::deque<Subgraph>& newFSGs);
        void mineFrequentSubgraphs(const GraphParameter& dataGraphInfo, IndexType& graphIndexes, const int& supportVal, std::map<int, std::deque<FSG>>& allFSG, int& noOfFsg);
        void dfsRecursiveMining(StackFSG& fsgClassStack, IsoMeasure& canonicalAllLevel, std::deque<Subgraph>& extendableEdges, int nE, const IsoTestMap& dataGraph, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::map<int, std::deque<FSG>>& allFSG, int& noOfFsg, std::deque<IndexType>& maximalFSG, FSGs& fsgs, int& tg);
    protected:
    private:

};

#endif // SUBGRAPHMINER_H
