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

#include "SubgraphMiner.h"

SubgraphMiner::SubgraphMiner()
{
    //ctor
}

SubgraphMiner::~SubgraphMiner()
{
    //dtor
}


/// Decompose the set into all possible (2^n-1) subsets and update to the collection of all multiedges.
void retrieveAllItems(std::set<int> inputSet, VecOfSet& allSubsets)
{
    /// Initilize the solution |
    std::set<int> init;
    auto it = inputSet.begin();
    init.insert(*it);
    allSubsets.push_back(init);

    if(inputSet.size() > 1){
        std::advance(it, 1);
        for(; it != inputSet.end(); ++it){
            /// Add the new current element as a new subset|
            std::set<int> single;
            single.insert(*it);
            allSubsets.push_back(single);
            /// Grow the already existing subsets with the new element |
            VecOfSet alreadyCollected = allSubsets;
            for(size_t k = 0; k < alreadyCollected.size(); ++k){
                alreadyCollected[k].insert(*it);
                allSubsets.push_back(alreadyCollected[k]);
            }
        }
    }
}

void collectAllDataEdges(const std::map<std::set<int>, int>& distinctDataEdges, std::set<std::set<int>>& distinctEdges)
{
    for(auto it = distinctDataEdges.begin(); it != distinctDataEdges.end(); ++it){
        VecOfSet allSubsets;
        retrieveAllItems(it->first, allSubsets);
        for(auto itS = allSubsets.begin(); itS != allSubsets.end(); ++itS)
            distinctEdges.insert(*itS);
    }
}


void SubgraphMiner::collectSubgraphSeeds(const GraphParameter& dataGraphInfo, std::deque<Subgraph>& subgraphSeeds)
{
    if(dataGraphInfo.nodeLabelsExist){ // vertex labels are accounted for |
        cout << "This version does not support vertex labels!" << endl;
        return;
    }

    /// Collect all distinct multiedges in the datagraph and count their number of occurrences.
    std::map<std::set<int>, int> distinctDataEdges; // <multiedge, its frequency>
    for(auto it = dataGraphInfo.eLabelMap.begin(); it != dataGraphInfo.eLabelMap.end(); ++it){
        auto itU = distinctDataEdges.find(it->second);
        if(itU == distinctDataEdges.end())
            distinctDataEdges.insert(std::make_pair(it->second, 1));
        else
            ++itU->second;
    }

    /// Collecting the sub-multiedges.
    std::set<std::set<int>> distinctEdges; // These edges are lexicographically ordered |
    collectAllDataEdges(distinctDataEdges, distinctEdges);

    /// Create structure for the edge subgraphSeeds
    for(auto it = distinctEdges.begin(); it != distinctEdges.end(); ++it){
        Subgraph subgraph;
        subgraph.edgeLabels.insert(std::make_pair(std::make_pair(0,1), *it)); // add the pattern with nodes 0 and 1 |
        subgraph.multiEdges.insert(*it);
        subgraphSeeds.push_back(subgraph);
    }
}


void SubgraphMiner::findFrequentSeeds(const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, std::deque<Subgraph>& newSubgraphs, std::deque<Subgraph>& newFSGs)
{
    Extension extender;
    for(size_t i = 0; i < newSubgraphs.size(); ++i){
        GraphParameter subgraphInfo;
        extender.extractInfo(newSubgraphs[i], subgraphInfo);
        SubgraphIso subIsomorphism;
        /// Find support using MIS on the fly by modifying SubgraphIso.
        if(subIsomorphism.findEmbeddingsMNI(supportVal, subgraphInfo, dataGraphInfo, graphIndexes, newSubgraphs[i]))
            newFSGs.push_back(newSubgraphs[i]);
    }
}

void collectFSG(const int& sSize, const std::deque<Subgraph> newFSGs, std::map<int, std::deque<FSG>>& allFSG)
{
    std::deque<FSG> newFSG;
    for(size_t i = 0; i < newFSGs.size(); ++i){
        FSG iNewFSG;
        iNewFSG.structure = newFSGs[i].edgeLabels;
        newFSG.push_back(iNewFSG);
    }
    auto it = allFSG.find(sSize);
    if(it == allFSG.end())
        allFSG.insert(std::make_pair(sSize, newFSG));
    else
        it->second.insert(it->second.end(), newFSG.begin(), newFSG.end());
}


void SubgraphMiner::dfsRecursiveMining(StackFSG& fsgClassStack, IsoMeasure& canonicalAllLevel, std::deque<Subgraph>& extendableEdges, int nE, const IsoTestMap& dataGraph, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::map<int, std::deque<FSG>>& allFSG, int& noOfFsg, std::deque<IndexType>& maximalFSG, FSGs& fsgs, int& tg)
{
    if (nE == extendableEdges.size())
        return; // go to next extension |

    while(!fsgClassStack.empty()){ // iterate until the entire DFS tree is spanned |
        std::deque<Subgraph> newFSGs;
        std::set<int> newEdge = extendableEdges[nE].edgeLabels.begin()->second;
		Extension extendor;
        extendor.subgraphExtension(dataGraph, fsgClassStack.front().front(), newEdge, canonicalAllLevel, supportVal, dataGraphInfo, graphIndexes, nonFSG, newFSGs, fsgClassStack.front(), maximalFSG, fsgs);

        if (!newFSGs.empty()){ // Grow the DFS tree until we find frequent subgraphs ; else go back to the next stack element
            fsgClassStack.push_front(newFSGs); // Add the new FSG to the stack
            auto itN = fsgs.find(newFSGs.front().edgeLabels.size());
            if(itN == fsgs.end()){
                fsgs.insert(std::make_pair(newFSGs.front().edgeLabels.size(), newFSGs));
            }
            else{
                itN->second.insert(itN->second.end(), newFSGs.begin(), newFSGs.end());
            }
            collectFSG(newFSGs.front().edgeLabels.size(), newFSGs, allFSG);
            noOfFsg += newFSGs.size(); // count all the subgraphs
        }
        else{ // Adding "newEdge" to "oldFSG" did not produce newFSGs
            dfsRecursiveMining(fsgClassStack, canonicalAllLevel, extendableEdges, nE+1, dataGraph, supportVal, dataGraphInfo, graphIndexes, nonFSG, allFSG, noOfFsg, maximalFSG, fsgs, tg);

            if(!fsgClassStack.empty()){
                /// Collect maximal patterns [Efficient without it].

                if(fsgClassStack.size() >= tg){ // check if it is maximal pattern
                    tg = fsgClassStack.size();
                    if(fsgClassStack.front().front().edgeLabels.size() > 2){ // check if the pattern has at least 3 edges
                        // store only the index structure of maximal patterns
                        SubgraphIso subIso;
                        IndexType FSGindex;
                        subIso.buildIndexes(fsgClassStack.front().front().parameters, FSGindex);
                        maximalFSG.push_back(FSGindex);

                    }
                }

                if(fsgClassStack.front().size() == 1){
                    fsgClassStack.pop_front(); // remove the entire "deque"
                }
                else{
                    fsgClassStack.front().pop_front(); // remove only one deque-element
                }
            }

        }
        nE = 0; // reset the pointer for next Extension - nE
    }
}

void findFsgSeeds(const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, std::deque<Subgraph>& newSubgraphs, std::deque<Subgraph>& newFSGs)
{
  if(dataGraphInfo.nodeLabelsExist){ // vertex labels are accounted for |
      cout << "This version does not support vertex labels!" << endl;
      exit;
  }

  /// Collect all distinct multiedges in the datagraph and count their number of occurrences.
  std::map<std::set<int>, int> distinctDataEdges; // <multiedge, its frequency>
  for(auto it = dataGraphInfo.eLabelMap.begin(); it != dataGraphInfo.eLabelMap.end(); ++it){
      auto itU = distinctDataEdges.find(it->second);
      if(itU == distinctDataEdges.end())
          distinctDataEdges.insert(std::make_pair(it->second, 1));
      else
          ++itU->second;
  }

  /// Collecting the sub-multiedges.
  std::set<std::set<int>> distinctEdges; // These edges are lexicographically ordered |
  collectAllDataEdges(distinctDataEdges, distinctEdges);

  /// Create structure for the edge subgraphSeeds
  for(auto it = distinctEdges.begin(); it != distinctEdges.end(); ++it){
      Subgraph subgraph;
      subgraph.edgeLabels.insert(std::make_pair(std::make_pair(0,1), *it)); // add the pattern with nodes 0 and 1 |
      subgraph.multiEdges.insert(*it);
      newSubgraphs.push_back(subgraph);
  }
}

void SubgraphMiner::mineFrequentSubgraphs(const GraphParameter& dataGraphInfo,  IndexType& graphIndexes, const int& supportVal, std::map<int, std::deque<FSG>>& allFSG, int& noOfFsg)
{
    /// Create initial seeds for mining.
    std::deque<Subgraph> subgraphSeeds, fsgSeeds;
//findFsgSeeds(supportVal, dataGraphInfo, graphIndexes, subgraphSeeds, fsgSeeds);
    clock_t Time = clock();
    collectSubgraphSeeds(dataGraphInfo, subgraphSeeds);
    double seedTime = double(clock()-Time)/CLOCKS_PER_SEC;
    // cout << "Subgraph seeds: " << subgraphSeeds.size() << endl;
    // cout << "Time: " << seedTime << endl;
    Time = clock();
    findFrequentSeeds(supportVal, dataGraphInfo, graphIndexes, subgraphSeeds, fsgSeeds);
    double totalTime = double(clock()-Time)/CLOCKS_PER_SEC;
    cout << "Frequent seeds: " << fsgSeeds.size() << endl;
    // cout << "Time: " << totalTime << endl;

    /// Collect the discovered FSG seed for outputting.
    collectFSG(1, fsgSeeds, allFSG);
    std::deque<FSG> frequentSubgraphs; // for better ordering of patterns use above data structure |
    noOfFsg = fsgSeeds.size(); // initialize the frequent subgraphs to the number of size-1 frequent subgraphs |

    /// Maintain a multiedge and vertex attribute mapping required for canonical representation of frequent subgraphs.
    IsoTestMap dataGraph;
    std::vector<int> edgeSize(fsgSeeds.size());
    for(size_t i = 0; i < fsgSeeds.size(); ++i){
        std::set<int> freqEdge = fsgSeeds[i].edgeLabels.begin()->second;
        dataGraph.multiedgeMap.insert(std::make_pair(freqEdge, i+1));
        edgeSize[i] = freqEdge.size();
    }

    /// Order the frequent seeds with their increasing size.
    std::vector<int> ind = sortIndexIncr(edgeSize);
    std::deque<Subgraph> fsgSeedsSorted(fsgSeeds.size());
    for(size_t i = 0; i < fsgSeeds.size(); ++i)
        fsgSeedsSorted[i] = fsgSeeds[ind[i]];

    std::deque<IndexType> maximalFSG;
    FSGs fsgs; // For invalid candidates.

    std::cout.flush() << "Iteration No: " ;
    for (size_t iS = 0; iS < fsgSeedsSorted.size(); ++iS) { // iS -> subgraph iterator |
		std::cout.flush() << iS << "...";
        /// Add only the lexicographically GE (>=) elements.
        std::deque<Subgraph> extendableEdges(fsgSeedsSorted.begin()+iS, fsgSeedsSorted.end());
        /// Add all elements.
        IsoMeasure canonicalAllLevel; // For the DFS subtree that grows from "fsgSeedsSorted[iS]" |
        StackFSG fsgClassStack; // FSG-class stack for DFS tree |
        fsgClassStack.push_front(std::deque<Subgraph>(fsgSeedsSorted.begin()+iS, fsgSeedsSorted.begin()+iS+1));
        NFSG nonFSG;
        int tg = 0;
        int nE = 0;
        dfsRecursiveMining(fsgClassStack, canonicalAllLevel, extendableEdges, nE, dataGraph, supportVal, dataGraphInfo, graphIndexes, nonFSG, allFSG, noOfFsg, maximalFSG, fsgs, tg);
    }
    cout << endl << "Done." << endl << "# FSG: " << noOfFsg << endl;
}
