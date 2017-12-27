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

#include "include/File.h"
#include "include/Extension.h"
#include "include/Frequency.h"
#include "include/SubgraphIso.h"
#include "include/GraphIso.h"
#include "include/SubgraphMiner.h"

int main(int argc, char* argv[])
{
    std::string filePath = argv[1];
    int support = std::stoi(argv[2]);

    /// Read data graph
    std::string nodeFile = filePath + "nodes.txt";
    std::string edgeFile = filePath + "edges.txt";
    GraphParameter dataGraphInfo;
    File dataGraph;
    dataGraph.readContents(nodeFile, edgeFile, dataGraphInfo);

    /// Build data graph indexes
    SubgraphIso subIsomorphism;
    IndexType graphIndexes;
    subIsomorphism.buildIndexes(dataGraphInfo, graphIndexes);

    /// Mine the multigraph
    std::map<int, std::deque<FSG>> allFSG;
    int noOfFsg;
    SubgraphMiner subgraphMiner;
    clock_t Time = clock();
    subgraphMiner.mineFrequentSubgraphs(dataGraphInfo, graphIndexes, support, allFSG, noOfFsg);
    double totalTime = double(clock()-Time)/CLOCKS_PER_SEC;
    cout << "Total time: " << totalTime << endl;
    // cout << "No. of levels: " << allFSG.size() << endl;

    /// Output results
    File results;
    results.printOutput(support, noOfFsg, totalTime, filePath);
    if(strcmp(argv[2],"-yes")) // output FSGs to a file
        results.printFSG(allFSG, filePath);
    cout << "FSG summary (Level: #FSG)" << endl;
    for(auto it = allFSG.begin(); it != allFSG.end(); ++it)
        cout << it->first << ": " << it->second.size() << endl;
}
