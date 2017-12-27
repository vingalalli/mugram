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

#ifndef FREQUENCY_H
#define FREQUENCY_H

#include "File.h"
#include "SubgraphIso.h"

class Frequency
{
    public:
        Frequency();
        virtual ~Frequency();
            bool isFrequent(const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, Subgraph& newSubgraph);
            bool isPreviousNonFrequent(NFSG& nonFSG, Subgraph& newSubgraph, IndexType& newSubgraphInd);
            bool isAlreadyFSG(std::deque<IndexType>& maximalFSG, Subgraph& newSubgraph);
    protected:
    private:
};

#endif // FREQUENCY_H
