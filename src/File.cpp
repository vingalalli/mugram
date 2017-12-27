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

#include "File.h"

File::File()
{
    //ctor
}

File::~File()
{
    //dtor
}

void File::splitString(const std::string& str, char chr, std::vector<std::string>& strs)
{
    std::string::const_iterator first = str.cbegin();
    std::string::const_iterator second = std::find(first+1, str.cend(), chr);
    while(second != str.cend()){
        strs.emplace_back(first, second);
        first = second+1;
        second = std::find(second+1, str.cend(), chr);
    }
    strs.emplace_back(first, str.cend());
}

void File::readContents(const std::string& nodeFile, const std::string& edgeFile, GraphParameter& graphInfo)
{
    /// Read node file
    graphInfo.nodeLabelsExist = false; // Do not consider node labels
    const char * nF = nodeFile.c_str();
    std::ifstream nFile (nF);
    if (nFile.is_open()) {
        std::string nodeAtt;
        while (getline(nFile, nodeAtt)) {
            std::vector<std::string> nodeAttAll;
            splitString(nodeAtt, ',', nodeAttAll);
            std::set<int> setAtt;
            for (size_t i = 0; i < nodeAttAll.size(); ++i)
                setAtt.insert(stoi(nodeAttAll[i]));
            graphInfo.attributes.push_back(setAtt);
        }
        nFile.close();
        graphInfo.nNodes = graphInfo.attributes.size();
      //        graphInfo.nodeLabelsExist = true;
    }
    else{
        graphInfo.nodeLabelsExist = false;
        std::cout << "Unable to open node file!" << std::endl;
    }
    /// Read edge file
    const char * eF = edgeFile.c_str();
    std::ifstream eFile (eF);
    std::map<std::string, int> node_m; // < original_id, mapped_id >

    if (eFile.is_open()) {
        std::string multiEdge;
        int n = 0;
        std::vector<std::vector<int>> adjListTemp(graphInfo.nNodes);
        EdgeLabel edgeLabels(graphInfo.nNodes);
        while (getline(eFile, multiEdge)) {
            std::vector<std::string> edgeContent;
            splitString(multiEdge, ' ', edgeContent);
            int node1; int node2;
            auto it = node_m.find(edgeContent.at(0));
            if (it == node_m.end()) {
                node_m.insert(make_pair(edgeContent.at(0), n));
                node1 = n;
                ++n;
            }
            else
                node1 = it->second;
            it = node_m.find(edgeContent.at(1));
            if (it == node_m.end()) {
                node_m.insert(make_pair(edgeContent.at(1), n));
                node2 = n;
                ++n;
            }
            else
                node2 = it->second;

            graphInfo.edges[0].push_back(node1);
            graphInfo.edges[1].push_back(node2);

            std::vector<std::string> labels;
            splitString(edgeContent.at(2), ',', labels);
            std::set<int> temp;
            for(auto it = labels.begin(); it != labels.end(); ++it)
                temp.insert(stoi((*it)));

            edgeLabels[node1].push_back(temp);
            edgeLabels[node2].push_back(temp);

            // Frequency of multiedges |
            auto it_e = graphInfo.edgeFrequency.find(temp);
            if (it_e == graphInfo.edgeFrequency.end())
                graphInfo.edgeFrequency.insert(make_pair(temp, 1));
            else
                it_e->second++;

            std::pair<int, int> p_1 = std::make_pair(node1,node2);
            auto it_1 = graphInfo.eLabelMap.find(p_1);
            if (it_1 == graphInfo.eLabelMap.end())
                graphInfo.eLabelMap.insert(std::make_pair(p_1, temp));
            else
                for (auto it = temp.begin(); it != temp.end(); ++it)
                    it_1->second.insert(*it);
            std::pair<int, int> p_2 = std::make_pair(node2,node1);
            auto it_2 = graphInfo.eLabelMap.find(p_2);
            if (it_2 == graphInfo.eLabelMap.end())
                graphInfo.eLabelMap.insert(std::make_pair(p_2, temp));
            else
                for (auto it = temp.begin(); it != temp.end(); ++it)
                    it_2->second.insert(*it);

            if (find(adjListTemp[node1].begin(), adjListTemp[node1].end(), node2) == adjListTemp[node1].end())
                adjListTemp[node1].push_back(node2);
            if (find(adjListTemp[node2].begin(), adjListTemp[node2].end(), node1) == adjListTemp[node2].end())
                adjListTemp[node2].push_back(node1);
      //            adjListTemp[node1].push_back(node2);
      //            adjListTemp[node2].push_back(node1);
        }
        graphInfo.adjacencyList = adjListTemp;
        graphInfo.neighbourSign = edgeLabels;
    }
    else
      std::cout << "Unable to open edge file!" << std::endl;
}

void File::extractInfo(const NodeLabelMap& nodes, const EdgeLabelMap& edges, GraphParameter& graphInfo)
{
    for (auto it = nodes.begin(); it != nodes.end(); ++it) {
        graphInfo.nodes.push_back(it->first);
        graphInfo.attributes.push_back(it->second);
    }
    graphInfo.nNodes = graphInfo.nodes.size();

    graphInfo.eLabelMap = edges;
    Vector2D adjListTemp(graphInfo.nNodes);
    EdgeLabel edgeLabels(graphInfo.nNodes);
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        edgeLabels[it->first.first].push_back(it->second);
        edgeLabels[it->first.second].push_back(it->second);
        graphInfo.edges[0].push_back(it->first.first);
        graphInfo.edges[1].push_back(it->first.second);
        if (find(adjListTemp[it->first.first].begin(), adjListTemp[it->first.first].end(), it->first.second) == adjListTemp[it->first.first].end())
            adjListTemp[it->first.first].push_back(it->first.second);
        if (find(adjListTemp[it->first.second].begin(), adjListTemp[it->first.second].end(), it->first.first) == adjListTemp[it->first.second].end())
            adjListTemp[it->first.second].push_back(it->first.first);
    }
    graphInfo.adjacencyList = adjListTemp;
    graphInfo.neighbourSign = edgeLabels;
}

void File::printFSG(const std::map<int, std::deque<FSG>> &allFSG, std::string &filePath)
{
    ofstream outFsg;
    std::string patternFile = filePath + "fsg.txt";
    outFsg.open(patternFile);
    for(auto it = allFSG.begin(); it != allFSG.end(); ++it){
        for(auto it1 = it->second.begin(); it1 != it->second.end(); ++it1){
            for(auto it2 = (*it1).structure.begin(); it2 != (*it1).structure.end(); ++it2){
                outFsg << it2->first.first << "\t" << it2->first.second << "\t";
                int mSize = it2->second.size();
                for(auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3){
                    if(mSize <= 1)
                        outFsg << *it3 << endl;
                    else{
                        outFsg << *it3 << ",";
                        --mSize;
                    }
                }
            }
            outFsg << endl;            
        }
    }
}

void File::printOutput(int &support, int &noOfFsg, double totalTime, std::string &filePath)
{
    ofstream outFile;
    std::string oFile = filePath + "output.txt";
    outFile.open(oFile);
    // outFile.open(oFile, std::ios_base::app); // append the results to the existing file |
    outFile << "support\t" << "# FSG\t" << "time (secs)" << endl;
    outFile << support << "\t" << noOfFsg << "\t" << totalTime << endl;
}
