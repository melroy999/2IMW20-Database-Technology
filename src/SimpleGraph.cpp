//
// Created by Nikolay Yakovets on 2018-01-31.
//

#include <cmath>
#include "SimpleGraph.h"

SimpleGraph::SimpleGraph(uint32_t n)   {
    setNoVertices(n);
}

uint32_t SimpleGraph::getNoVertices() const {
    return V;
}

void SimpleGraph::setDataStructureSizes(bool isJoin) {
    E = 0;
    numEdges.resize(getNoLabels());
    sources.resize(getNoLabels(), std::vector<uint64_t>(static_cast<unsigned long>(std::ceil((double) V / 64))));
    targets.resize(getNoLabels(), std::vector<uint64_t>(static_cast<unsigned long>(std::ceil((double) V / 64))));

    adj.resize(getNoLabels());
    for (auto &l : adj)
        l.resize(V);

    if(!isJoin) {
        reverse_adj.resize(getNoLabels());
        for (auto &l : reverse_adj)
            l.resize(V);
    }
}

void SimpleGraph::setNoVertices(uint32_t n) {
    V = n;
}

uint32_t SimpleGraph::getNoEdges() const {
    return E;
}

uint32_t SimpleGraph::getNoDistinctEdges() const {
    return E;
}

uint32_t SimpleGraph::getNoLabels() const {
    return L;
}

void SimpleGraph::setNoLabels(uint32_t noLabels) {
    L = noLabels;
}

void SimpleGraph::addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) {
    if(from >= V || to >= V || edgeLabel >= L)
        throw std::runtime_error(std::string("Edge data out of bounds: ") +
                                         "(" + std::to_string(from) + "," + std::to_string(to) + "," +
                                         std::to_string(edgeLabel) + ")");

    sources[edgeLabel][from/64] |= SET_BIT(from % 64);
    targets[edgeLabel][to/64] |= SET_BIT(to % 64);

    E += 1;
    numEdges[edgeLabel] += 1;
    adj[edgeLabel][from].emplace_back(to);

    if(!reverse_adj.empty()) {
        reverse_adj[edgeLabel][to].emplace_back(from);
    }
}

void SimpleGraph::addEdges(std::shared_ptr<SimpleGraph> &in, uint32_t projectLabel, bool isInverse) {

    E = in->numEdges[projectLabel];

    if(!isInverse) {
        adj_ptr = &in->adj[projectLabel];
        sources_ptr = &in->sources[projectLabel];
        targets_ptr = &in->targets[projectLabel];
    } else {
        adj_ptr = &in->reverse_adj[projectLabel];
        sources_ptr = &in->targets[projectLabel];
        targets_ptr = &in->sources[projectLabel];
    }
}

void SimpleGraph::readFromContiguousFile(const std::string &fileName) {

    std::fstream file(fileName, std::ios_base::in);

    std::string graphStructure;
    file >> graphStructure;

    std::istringstream iss(graphStructure);
    std::string line;
    uint32_t noNodes = 0;
    uint32_t noEdges = 0;
    uint32_t noLabels = 0;
    int j = 0;
    while(std::getline(iss, line, ',')) {
        if(j == 0) {
            noNodes = (uint32_t) std::stoul(line);
        }  else if(j == 1) {
            noEdges = (uint32_t) std::stoul(line);
        } else if(j == 2) {
            noLabels = (uint32_t) std::stoul(line);
        }
        j++;
    }

    setNoLabels(noLabels);
    setNoVertices(noNodes);
    setDataStructureSizes(false);

    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edges;
    edges.reserve(noEdges);

    for(int i = 0; i < noEdges; i++) {
        uint32_t subject, predicate, object;

        file >> subject >> predicate >> object;
        file.ignore(INT32_MAX, '\n');

        // We place them in the order of importance. We want the order predicate > subject > object.
        edges.emplace_back(predicate, subject, object);
    }

    // Sort the vector.
    std::sort(edges.begin(), edges.end());

    // Add the edges in sorted order, do not add duplicates.
    std::tuple<uint32_t, uint32_t, uint32_t> lastTuple;

    bool first = true;
    uint32_t lastSource = 0;
    uint32_t lastTarget = 0;
    uint32_t lastLabel = 0;

    for(auto const &tuple : edges) {
        uint32_t source = std::get<1>(tuple);
        uint32_t target = std::get<2>(tuple);
        uint32_t label = std::get<0>(tuple);

        if(first || !(lastSource == source && lastTarget == target && lastLabel == label)) {
            first = false;
            lastSource = source;
            lastTarget = target;
            lastLabel = label;

            addEdge(source, target, label);
        }
    }

    // For the leaf level, we want the reverse adjacency to be sorted as well.
    for(auto &adj_list : reverse_adj) {
        for(auto &vertices : adj_list) {
            std::sort(vertices.begin(), vertices.end());
        }
    }

    file.close();
}