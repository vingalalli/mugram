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

#include "Trie.h"

Trie::Trie()
{
    root = new Node();
}

Trie::~Trie()
{
    destroyTrie(Root());
}

void Trie::assignParent(Node* trieNode)
{
    if (trieNode->children().empty())
        return;
    for (size_t i =0; i < trieNode->children().size(); ++i){
        trieNode->children().at(i)->setNodeParent(trieNode);
        assignParent(trieNode->children().at(i));
    }
    return;
}


void Trie::destroyTrie(Node* n) {
    Node* current = n;
    int i = 0;
    while (current->children().size() != i) {
       destroyTrie(current->children().at(i));
       i = i + 1;
    }
    delete current;
}

void Trie::preOrder(Node* n) {
    Node* current = n;
    cout << current->contentInt() << endl; // postOrder
    if(current->wordMarker())
        cout << ": word" << endl;
    else
        cout << endl;
    int i = 0;
    while (current->children().size() != i) {
       preOrder(current->children().at(i));
       i = i + 1;
    }
}

void Trie::postOrder(Node* n) {
    Node* current = n;
    int i = 0;
    while (current->children().size() != i) {
       preOrder(current->children().at(i));
       i = i + 1;
    }
    cout << current->contentInt() << endl; // postOrder
//    if(current->wordMarker())
//        cout << ": word" << endl;
//    else
//        cout << endl;
}

Node* Node::findChildInt(int c)
{
    for(int i = 0; i < mChildren.size(); i++){
        Node* tmp = mChildren.at(i);
        if (tmp->contentInt() == c)
            return tmp;
    }
    return NULL;
}

void Trie::updateHashTable(Node* trieNode, LabelHash& LabelMap)
{
    Node* current = trieNode;
    auto it = LabelMap.find(current->contentInt());
    if (it == LabelMap.end()){ // new label to be added
        std::vector<Node *> tmp;
        tmp.push_back(current);
        LabelMap.insert(std::make_pair(current->contentInt(), tmp));
    }
    else { // label already exists
        it->second.push_back(current);
    }
    int i = 0;
    while(current->children().size() != i){
       updateHashTable(current->children().at(i), LabelMap);
       i = i + 1;
    }
}

void Trie::addSignatureDim(const std::vector<int>& s, const int& id)
{
    Node* current = root;
    if (s.empty()){
        current->setWordMarker(); // an empty word
        return;
    }
    int i = 0;
    for(auto it = s.begin(); it!=s.end(); ++it){
        Node* child = current->findChildInt((*it)+1); // label = id +1;
        if(child != NULL){
            child->setVertexId(id);//add the vertex id for the already existing node label
            current = child;
        }
        else{
            Node* tmp = new Node();
            tmp->setContentInt((*it)+1); // label = id +1;
            tmp->setVertexId(id);// add anew node label and its vertex id
            current->appendChild(tmp);
            current = tmp;
        }
        if(i == s.size() - 1)
            current->setWordMarker();
        ++i;
    }
}


void Trie::addMultiedges(const std::set<int>& multiedge)
{
    Node* current = root;
//    current->setContentInt(0); // to check if the Trie is empty or not |
    if (multiedge.empty()){
        current->setWordMarker(); // an empty word
        return;
    }
    int i = 0;
    for(auto it = multiedge.begin(); it!=multiedge.end(); ++it) {
        Node* child = current->findChildInt((*it)+1); // label = id +1;
        if ( child != NULL ){
            current->countOccurrences();
            current = child;
        }
        else{
            Node* tmp = new Node();
            tmp->setContentInt((*it)+1); // label = id +1;
            tmp->assignMultiedge(multiedge);
            tmp->countOccurrences();
            current->appendChild(tmp);
            current = tmp;
        }
        if ( i == multiedge.size() - 1 ){
            current->setWordMarker();
//            cout << "set" << endl;
        }
        ++i;
    }
}

void Trie::addNonFreqEdges(const std::set<int>& multiedge)
{
    Node* current = root;
    if (multiedge.empty()){
        current->setWordMarker(); // an empty word
        return;
    }
    int i = 0;
    for(auto it = multiedge.begin(); it!=multiedge.end(); ++it) {
        Node* child = current->findChildInt((*it)+1); // label = id +1;
        if (child != NULL){
            current = child;
        }
        else{
            Node* tmp = new Node();
            tmp->setContentInt((*it)+1); // label = id +1;
            current->appendChild(tmp);
            current = tmp;
        }
        if (i == multiedge.size() - 1)
            current->setWordMarker();
        ++i;
    }
}

/// This is a subset containment approach.
/*
void Trie::IfEdgeExists(Node* n, std::set<int> s, bool& m_found){
	Node* current = n;
    if(current->contentInt() == (*s.begin())+1)
        s.erase(s.begin());
    if (s.empty()) {
        m_found = true;
    }
    else if (!m_found){
        int j = 0;
        while (current->children().size() != j) {
            IfEdgeExists(current->children().at(j), s, m_found);
            if (m_found)
                break;
            j = j + 1;
        }
    }
}
*/


/// This is a superset containment approach
bool Trie::edgeExists(Node* n, std::set<int> s, int pos, bool& found){
//    if(pos >= s.size() && !found) //
//        return false;
    auto it = s.begin();
    std::advance(it, pos);
    int i = 0;

    while (i < n->children().size() && !found) {
        Node* child = n->children().at(i);
        if(child->contentInt() == (*it)+1){
            if(child->wordMarker()){
                found = true;
//                return true;
            }
            else if(pos < s.size()-1){
                ++pos;
                edgeExists(child, s, pos, found);
            }
        }
        ++i;
    }
//    return false;
}
