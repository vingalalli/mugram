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

#ifndef MATCH_H
#define MATCH_H

#include "../include/File.h"
#include "Index.h"

struct Flags
{
    bool embFound;
    bool satBinFilled = false;
    bool validRule;
    bool invalidMatch;
    bool temp = false; // needed for some tests.
};

struct StarSolutions
{
    Vector coreSolutions;
    Vector3D satSolutions;
};

struct InitialVertex
{
    int pattern;
    int data;
};

struct Solutions
{
    Vector2D coreSolutions; // each row is a core vertex; if coreSolutions.size() > 1 => there are more than one core vertices with satellite solutoins |
    Vector3D sattSolutions; //
};

class Match
{
    public:
        Match();
        virtual ~Match();
        void orderVertices(const int& queryNodes, const EdgeLabel& queryNeighbourSign, const Vector2D& queryAdjacencyList, std::vector<int>& querySequence);
        void orderVerticesAllMNI(const int& queryNodes, const EdgeLabel& queryNeighbourSign, const Vector2D& queryAdjacencyList, std::vector<int>& orderedNodes, int& initialVertex);
        void orderStarNodes(const EdgeLabel& neighbourSign, std::vector<int>& orderedQuery);
        bool findMatchesFastMNI(const int& freq, const std::vector<int>& initialMatches, const GraphParameter& queryGraphInfo, const std::vector<int>& querySequence, const VecOfSet& nodeMatches, const std::vector<Trie*>& dataNeighTrie, std::map<int, std::vector<int>>& invalidCands, std::vector<std::unordered_set<int>>&  nodeRepository, bool& notFrequent);
        bool findStarMatches(const int& freq, const Vector& initialMatches, GraphParameter& queryGraphInfo, const Vector& orderedQuery, const VecOfSet& nodeMatches, const std::vector<Trie*>& dataNeighTrie, InvalidCands& invalidCands, NodeRepository&  nodeRepository);
        bool checkSI(std::vector<int>& initialMatches, const GraphParameter& queryGraphInfo, const std::vector<int>& querySequence, std::vector<int>& matchedDataVertices, const VecOfSet& nodeMatches, const std::vector<Trie*>& dataNeighTrie);
    protected:
    private:
        void chooseFrontier(const std::vector<int>& already_m, const Vector2D& queryAdjacencyList, std::vector<int>& frontier);
};

#endif // MATCH_H
