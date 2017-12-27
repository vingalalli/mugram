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

#include "Extension.h"

Extension::Extension()
{
		//ctor
}

Extension::~Extension()
{
  	//dtor
}


void Extension::extractInfo(const Subgraph& subgraph, GraphParameter& subgraphInfo)
{
/// Node labels exist
	/*
    Update: subgraphInfo.attributes
    */
///
	  std::set<int> nodes;
	  for(auto it = subgraph.edgeLabels.begin(); it != subgraph.edgeLabels.end(); ++it){
		    nodes.insert(it->first.first);
		    nodes.insert(it->first.second);
	  }
	  subgraphInfo.nNodes = nodes.size();

	  subgraphInfo.eLabelMap = subgraph.edgeLabels; // add the original edges |
	  for(auto it = subgraph.edgeLabels.begin(); it != subgraph.edgeLabels.end(); ++it){ // add the reversed edges, for the undirected graphs |
		    std::pair<int, int> edgeRev = std::make_pair(it->first.second, it->first.first);
		    subgraphInfo.eLabelMap.insert(std::make_pair(edgeRev, it->second));
	  }

	  Vector2D adjListTemp(subgraphInfo.nNodes);
	  EdgeLabel nbrSignTemp(subgraphInfo.nNodes);
	  for(auto it = subgraph.edgeLabels.begin(); it != subgraph.edgeLabels.end(); ++it){
		    nbrSignTemp[it->first.first].push_back(it->second);
		    nbrSignTemp[it->first.second].push_back(it->second);
		    if (find(adjListTemp[it->first.first].begin(), adjListTemp[it->first.first].end(), it->first.second) == adjListTemp[it->first.first].end())
		      	adjListTemp[it->first.first].push_back(it->first.second);
		    if (find(adjListTemp[it->first.second].begin(), adjListTemp[it->first.second].end(), it->first.first) == adjListTemp[it->first.second].end())
		      	adjListTemp[it->first.second].push_back(it->first.first);
	  }
	  subgraphInfo.adjacencyList = adjListTemp;
	  subgraphInfo.neighbourSign = nbrSignTemp;
}

bool Extension::isIsomorphicCurrent(Subgraph& newSubgraph, std::set<std::vector<int>>& canonicalSet, const IsoTestMap& dataGraph)
{
	  // Call TestIsoInvariants to avoid generating canonical OR compute automorphic group
	  extractInfo(newSubgraph, newSubgraph.parameters);
	  GraphIso subgraph;
	  subgraph.generateCanonical(newSubgraph.parameters, newSubgraph.canonicalForm, dataGraph);
	  if(canonicalSet.find(newSubgraph.canonicalForm) == canonicalSet.end()){
		    canonicalSet.insert(newSubgraph.canonicalForm);
		    return false;
	  }
	  else
	    	return true;
}

bool Extension::isIsomorphicPrevious(const int& subgraphSize, const std::vector<int>& canonical, IsoMeasure& canonicalAllLevel)
{
	  auto it = canonicalAllLevel.find(subgraphSize);
	  if(it == canonicalAllLevel.end()){ // Create a new level and add canonical elements to it
		    std::set<std::vector<int>> tmp;
		    tmp.insert(canonical);
		    canonicalAllLevel.insert(std::make_pair(subgraphSize, tmp));
		    return false;
	  }
	  else{
		    if(it->second.find(canonical) == it->second.end()){ // The subgraph in the same level does not exist
			      it->second.insert(canonical);
			      return false;
		    }
		    else
		      	return true; // Subgraph already found
	  }
}

inline bool notAvalidExtension(const Set2D& nonFreqEdges, const std::set<int>& newEdge)
{
	  for(auto it = nonFreqEdges.begin(); it != nonFreqEdges.end(); ++it)
	    	if(std::includes(newEdge.begin(), newEdge.end(), (*it).begin(), (*it).end()))
	        	return true;
	  return false;
}


void Extension::subgraphExtension(const IsoTestMap& dataGraph, Subgraph& oldFSG, const set<int>& newEdge, IsoMeasure& canonicalAllLevel, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::deque<Subgraph>& newFSGs, std::deque<Subgraph> oldFSGs, std::deque<IndexType>& maximalFSG, FSGs& fsgs)
{
    if(oldFSG.alreadyExtended){ // If the subgraph has already been extended
        // Ext 1. Adding an edge WITH an extra node.
        extendWithNode(dataGraph, oldFSG, newEdge, canonicalAllLevel, supportVal, dataGraphInfo, graphIndexes, nonFSG, newFSGs, oldFSGs, maximalFSG, fsgs);
        // Ext 2. Adding an edge WITHOUT an extra node.
        extendWithoutNode(dataGraph, oldFSG, newEdge, canonicalAllLevel, supportVal, dataGraphInfo, graphIndexes, nonFSG, newFSGs, oldFSGs, maximalFSG, fsgs);
    }
    else{
        // Ext 1. Adding an edge WITH an extra node.
        extendWithNodeNew(dataGraph, oldFSG, newEdge, canonicalAllLevel, supportVal, dataGraphInfo, graphIndexes, nonFSG, newFSGs, oldFSGs, maximalFSG, fsgs);
        // Ext 2. Adding an edge WITHOUT an extra node.
        extendWithoutNodeNew(dataGraph, oldFSG, newEdge, canonicalAllLevel, supportVal, dataGraphInfo, graphIndexes, nonFSG, newFSGs, oldFSGs, maximalFSG, fsgs);
        // Comment the below line to UNPRUNE the unfrequent-edge-type-extension approach.
        oldFSG.alreadyExtended = true; // "oldFSG" has just been extended.
    }
}


void Extension::extendWithNode(const IsoTestMap& dataGraph, Subgraph& oldFSG, const set<int>& newEdge, IsoMeasure& canonicalAllLevel, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::deque<Subgraph>& newFSGs, std::deque<Subgraph> oldFSGs, std::deque<IndexType>& maximalFSG, FSGs& fsgs)
{
    for(auto it = oldFSG.nonFreqEdgesN.begin(); it != oldFSG.nonFreqEdgesN.end(); ++it){
        Set2D nonFreqEdges = it->second; // "nonFreqEdges" specific to the possible extension for node *it
        if(!nonFreqEdges.empty()) // Non-frequent edge extensions are available for this extension |
            if(notAvalidExtension(nonFreqEdges, newEdge))
                continue; // try the next extension |

        std::pair<int, int> tmpEdgeID = std::make_pair(it->first, oldFSG.orderedNodes.size());
        EdgeLabelMap edgeLabels = oldFSG.edgeLabels;
        if(edgeLabels.find(tmpEdgeID) == edgeLabels.end())
            edgeLabels.insert(std::make_pair(tmpEdgeID, newEdge));
        Set2D multiEdges = oldFSG.multiEdges;
        multiEdges.insert(newEdge); // maintain the set of multiedges |

        Subgraph newSubgraph;
        newSubgraph.edgeLabels = edgeLabels;
        newSubgraph.multiEdges = multiEdges;
        extractInfo(newSubgraph, newSubgraph.parameters);
        GraphIso subgraph;
				Frequency frequency;
        subgraph.generateCanonical(newSubgraph.parameters, newSubgraph.canonicalForm, dataGraph);
        if(frequency.isAlreadyFSG(maximalFSG, newSubgraph))
            continue;
        if(isIsomorphicPrevious(newSubgraph.edgeLabels.size(), newSubgraph.canonicalForm, canonicalAllLevel)){
                continue;
        }

        SubgraphIso subIso;
        IndexType newSubgraphInd;
        subIso.buildIndexes(newSubgraph.parameters, newSubgraphInd);
        if(frequency.isPreviousNonFrequent(nonFSG, newSubgraph, newSubgraphInd)){
            continue;
        }

        getInvalidCands(oldFSGs, fsgs, newSubgraph, newSubgraphInd);

        if(frequency.isFrequent(supportVal, dataGraphInfo, graphIndexes, newSubgraph)){
            newFSGs.push_back(newSubgraph);
        }
        else{
            it->second.insert(newEdge);
            /// Collect non-frequent SGs
            auto itN = nonFSG.find(newSubgraph.edgeLabels.size());
            if(itN == nonFSG.end()){
                std::deque<Subgraph> tmp;
                tmp.push_back(newSubgraph);
                nonFSG.insert(std::make_pair(newSubgraph.edgeLabels.size(), tmp));
            }
            else{
                // Always ONLY a NEW non-frequent subgraph is added. No need to check again |
                itN->second.push_back(newSubgraph);
            }
        }
    }
}

void Extension::extendWithoutNode(const IsoTestMap& dataGraph, Subgraph& oldFSG, const set<int>& newEdge, IsoMeasure& canonicalAllLevel, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::deque<Subgraph>& newFSGs, std::deque<Subgraph> oldFSGs, std::deque<IndexType>& maximalFSG, FSGs& fsgs)
{
    for(auto it = oldFSG.nonFreqEdgesE.begin(); it != oldFSG.nonFreqEdgesE.end(); ++it){
        Set2D nonFreqEdges = it->second; // "nonFreqEdges" specific to the possible extension for edge *it
        if(!nonFreqEdges.empty()){ // Non-frequent edge extensions are available for this extension |
            if(notAvalidExtension(nonFreqEdges, newEdge))
                continue; // try the next extension |
        }
        std::pair<int,int> newEdgeIDRev = std::make_pair(it->first.second, it->first.first);
        if(oldFSG.edgeLabels.find(it->first) != oldFSG.edgeLabels.end() || oldFSG.edgeLabels.find(newEdgeIDRev) != oldFSG.edgeLabels.end())
            continue;

        EdgeLabelMap edgeLabels = oldFSG.edgeLabels;
        edgeLabels.insert(std::make_pair(it->first, newEdge)); // add a new edge: // add the new edge only in one direction "newEdgeID" |
        Set2D multiEdges = oldFSG.multiEdges;
        multiEdges.insert(newEdge); // maintain the set of multiedges |

        Subgraph newSubgraph;
        newSubgraph.edgeLabels = edgeLabels;
        newSubgraph.multiEdges = multiEdges;
        extractInfo(newSubgraph, newSubgraph.parameters);
        GraphIso subgraph;
				Frequency frequency;
        subgraph.generateCanonical(newSubgraph.parameters, newSubgraph.canonicalForm, dataGraph);
        if(frequency.isAlreadyFSG(maximalFSG, newSubgraph))
            continue;
        if(isIsomorphicPrevious(newSubgraph.edgeLabels.size(), newSubgraph.canonicalForm, canonicalAllLevel)){
            continue;
        }

        SubgraphIso subIso;
        IndexType newSubgraphInd;
        subIso.buildIndexes(newSubgraph.parameters, newSubgraphInd);

        if(frequency.isPreviousNonFrequent(nonFSG, newSubgraph, newSubgraphInd)){
            continue;
        }

        getInvalidCands(oldFSGs, fsgs, newSubgraph, newSubgraphInd);


        if(frequency.isFrequent(supportVal, dataGraphInfo, graphIndexes, newSubgraph)){
            newFSGs.push_back(newSubgraph);
        }
        else{
            it->second.insert(newEdge);
            /// Collect non-frequent SGs
            auto itN = nonFSG.find(newSubgraph.edgeLabels.size());
            if(itN == nonFSG.end()){
                std::deque<Subgraph> tmp;
                tmp.push_back(newSubgraph);
                nonFSG.insert(std::make_pair(newSubgraph.edgeLabels.size(), tmp));
            }
            else{
                // Always ONLY a NEW non-frequent subgraph is added. No need to check again |
                itN->second.push_back(newSubgraph);
            }
        }
    }
}

void Extension::extendWithNodeNew(const IsoTestMap& dataGraph, Subgraph& oldFSG, const set<int>& newEdge, IsoMeasure& canonicalAllLevel, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::deque<Subgraph>& newFSGs, std::deque<Subgraph> oldFSGs, std::deque<IndexType>& maximalFSG, FSGs& fsgs)
{
    Set2D temp;
    int nNodes = oldFSG.orderedNodes.size(); // no of nodes in the input subgraph |
    std::set<std::vector<int>> extraNodeCanonical;
    for(auto it = oldFSG.orderedNodes.begin(); it != oldFSG.orderedNodes.end(); ++it){
        EdgeLabelMap edgeLabels = oldFSG.edgeLabels;
        std::pair<int, int> tmpEdgeID = std::make_pair(*it, nNodes);
        if(edgeLabels.find(tmpEdgeID) == edgeLabels.end())
            edgeLabels.insert(std::make_pair(tmpEdgeID, newEdge)); // add a new edge: a new node "nNodes" forming an edge with an already existing node "*it" |
        Set2D multiEdges = oldFSG.multiEdges;
        multiEdges.insert(newEdge); // maintain the set of multiedges |

        Subgraph newSubgraph;
        newSubgraph.edgeLabels = edgeLabels;
        newSubgraph.multiEdges = multiEdges;
        if(isIsomorphicCurrent(newSubgraph, extraNodeCanonical, dataGraph)){ // Check if this new set of edges "edgeLabels" are isomorphic to the already existing "newSubgraphs"
            continue;
        }
				Frequency frequency;
        if(frequency.isAlreadyFSG(maximalFSG, newSubgraph))
            continue;

        oldFSG.nonFreqEdgesN.insert(std::make_pair((*it), temp)); // This step assigns a set of unique nodes that are to be extended. Each node is selected on a first-come-first-serve basis, which is a representative node from its permutation group.

        if(isIsomorphicPrevious(newSubgraph.edgeLabels.size(), newSubgraph.canonicalForm, canonicalAllLevel)){

                continue;
        }

        SubgraphIso subIso;
        IndexType newSubgraphInd;
        subIso.buildIndexes(newSubgraph.parameters, newSubgraphInd);

        if(frequency.isPreviousNonFrequent(nonFSG, newSubgraph, newSubgraphInd)){
            continue;
        }
        getInvalidCands(oldFSGs, fsgs, newSubgraph, newSubgraphInd);

        if(frequency.isFrequent(supportVal, dataGraphInfo, graphIndexes, newSubgraph)){
            newFSGs.push_back(newSubgraph);
        }
        else{
            oldFSG.nonFreqEdgesN.find((*it))->second.insert(newEdge);
            /// Collect non-frequent SGs
            auto itN = nonFSG.find(newSubgraph.edgeLabels.size());
            if(itN == nonFSG.end()){
                std::deque<Subgraph> tmp;
                tmp.push_back(newSubgraph);
                nonFSG.insert(std::make_pair(newSubgraph.edgeLabels.size(), tmp));
            }
            else{
                // Always ONLY a NEW non-frequent subgraph is added. No need to check again |
                itN->second.push_back(newSubgraph);
            }
        }
    }
}


void Extension::extendWithoutNodeNew(const IsoTestMap& dataGraph, Subgraph& oldFSG, const set<int>& newEdge, IsoMeasure& canonicalAllLevel, const int& supportVal, const GraphParameter& dataGraphInfo, IndexType& graphIndexes, NFSG& nonFSG, std::deque<Subgraph>& newFSGs, std::deque<Subgraph> oldFSGs, std::deque<IndexType>& maximalFSG, FSGs& fsgs)
{
    Set2D temp;
    std::set<std::pair<int,int>> extendableEdges;
    std::set<std::vector<int>> noNodeCanonical;
    for(size_t i = 0; i < oldFSG.orderedNodes.size()-1; ++i){
        for(size_t j = i+1; j < oldFSG.orderedNodes.size(); ++j){
            std::pair<int,int> newEdgeID = std::make_pair(oldFSG.orderedNodes[i], oldFSG.orderedNodes[j]); // a possible edge in the subgraph clique |
            std::pair<int,int> newEdgeIDRev = std::make_pair(oldFSG.orderedNodes[j], oldFSG.orderedNodes[i]); // a possible edge in the subgraph clique |
            if(oldFSG.edgeLabels.find(newEdgeID) != oldFSG.edgeLabels.end() || oldFSG.edgeLabels.find(newEdgeIDRev) != oldFSG.edgeLabels.end())
                continue;
            EdgeLabelMap edgeLabels = oldFSG.edgeLabels;
            edgeLabels.insert(std::make_pair(newEdgeID, newEdge)); // add a new edge: // add the new edge only in one direction "newEdgeID" |
            Set2D multiEdges = oldFSG.multiEdges;
            multiEdges.insert(newEdge); // maintain the set of multiedges |

            Subgraph newSubgraph;
            newSubgraph.edgeLabels = edgeLabels;
            newSubgraph.multiEdges = multiEdges;
            if(isIsomorphicCurrent(newSubgraph, noNodeCanonical, dataGraph)){ // Check if this new set of edges "edgeLabels" are isomorphic to the already existing "newSubgraphs"
                continue;
            }
						Frequency frequency;
            if(frequency.isAlreadyFSG(maximalFSG, newSubgraph))
                continue;

            oldFSG.nonFreqEdgesE.insert(std::make_pair(newEdgeID, temp)); // This step assigns a set of unique node-pairs to which an edge can be added. Each node-pair is selected on a first-come-first-serve basis, which is a representative node-pair from its permutation group.
            if(isIsomorphicPrevious(newSubgraph.edgeLabels.size(), newSubgraph.canonicalForm, canonicalAllLevel)){
                continue;
            }

        SubgraphIso subIso;
        IndexType newSubgraphInd;
        subIso.buildIndexes(newSubgraph.parameters, newSubgraphInd);

            if(frequency.isPreviousNonFrequent(nonFSG, newSubgraph, newSubgraphInd)){
                continue;
            }

            getInvalidCands(oldFSGs, fsgs, newSubgraph, newSubgraphInd);

            if(frequency.isFrequent(supportVal, dataGraphInfo, graphIndexes, newSubgraph)){
                newFSGs.push_back(newSubgraph);
            }
            else{
                oldFSG.nonFreqEdgesE.find(newEdgeID)->second.insert(newEdge);
                /// Collect non-frequent SGs
                auto itN = nonFSG.find(newSubgraph.edgeLabels.size());
                if(itN == nonFSG.end()){
                    std::deque<Subgraph> tmp;
                    tmp.push_back(newSubgraph);
                    nonFSG.insert(std::make_pair(newSubgraph.edgeLabels.size(), tmp));
                }
                else{
                    // Always ONLY a NEW non-frequent subgraph is added. No need to check again |
                    itN->second.push_back(newSubgraph);
                }
            }
        }
    }
}

void fetchInvalids(Subgraph& newSubgraph, Subgraph& oldSubgraph, IndexType& newSubgraphInd)
{
	  SubgraphIso subIso;
	  Vector2D solutions; // Holds a set of isomorphic solutions for "oldSubgraph.orderedNodes" |
	  if(!subIso.findIsoMatch(oldSubgraph, newSubgraphInd, solutions))
				return;
	  for(size_t i = 0; i < solutions.size(); ++i){
		    for(auto itO = oldSubgraph.orderedNodes.begin(), itN = solutions[i].begin(); itO != oldSubgraph.orderedNodes.end(); ++itO, ++itN){
			      auto itIO = oldSubgraph.invalidCands.find(*itO);
			      auto itIN = newSubgraph.invalidCands.find(*itN);
			      if(itIN != newSubgraph.invalidCands.end() && itIO != oldSubgraph.invalidCands.end())
			        	itIN->second = vectorUnion(itIO->second , itIN->second);
		    }
	  }
}

void Extension::getInvalidCands(std::deque<Subgraph>& oldFSGs, FSGs& fsgs, Subgraph& newSubgraph, IndexType& newSubgraphInd)
{
		/// Propagate the invalid candidates directly form a parent to its child.
	  newSubgraph.invalidCands = oldFSGs.front().invalidCands;

	  if(newSubgraph.edgeLabels.size() <= 2) // oldFSG has to be at least size 2 |
				return;
	  /// Propagate the invalid candidates from already discovered FSG that are subgraphs of "newSubgraph".
	  for(auto it = fsgs.begin(); it != fsgs.end(); ++it){
	    	if(it->first >= newSubgraph.edgeLabels.size())
	      		continue;
		    for(size_t i = 0; i < it->second.size(); ++i)
		    		fetchInvalids(newSubgraph, it->second[i], newSubgraphInd);
	  }
	  /// Propagate invalid candidates from the class of patterns.
	  for(size_t i = 1; i < oldFSGs.size(); ++i)
	    	fetchInvalids(newSubgraph, oldFSGs[i], newSubgraphInd);
}
