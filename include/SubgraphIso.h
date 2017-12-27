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

#ifndef SUBGRAPHISO_H
#define SUBGRAPHISO_H

#include "File.h"
#include "../sumgra/Index.h"
#include "../sumgra/Match.h"

#define MAX_ALLOWED_TIME 0.1
#define TOTAL_MAX_ALLOWED_TIME 60 // in seconds

struct IndexType{
    AttMap attributeHash;
    RTree synopsesTrie;
    std::vector<Trie*> neighborTrie;
    // IndexType() : dataBitSet(dataGraphInfo.nNodes) {}
    // EdgeLabelBit dataBitSet;
};

class SubgraphIso
{
    public:
        SubgraphIso();
        virtual ~SubgraphIso();
        void buildIndexes(GraphParameter& dataGraphInfo, IndexType& graphIndexes);
        bool findEmbeddingsMNI(const int& freq, GraphParameter& queryGraphInfo, const GraphParameter& dataGraphInfo, IndexType& indexType, Subgraph& matchedSubgraph);
        bool findIsoMatch(Subgraph& query, IndexType& indexType, Vector2D& solutions);
    protected:
    private:
};

#endif // SUBGRAPHISO_H
