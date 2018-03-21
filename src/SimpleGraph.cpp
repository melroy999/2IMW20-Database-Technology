//
// Created by Nikolay Yakovets on 2018-01-31.
//

#include "SimpleGraph.h"

SimpleGraph::SimpleGraph(uint32_t n)   {
    setNoVertices(n);
}

uint32_t SimpleGraph::getNoVertices() const {
    return V;
}

void SimpleGraph::setDataStructureSizes() {
    adj.resize(getNoLabels());
    reverse_adj.resize(getNoLabels());

    for (auto &l : adj)
        l.resize(V);

    for (auto &l : reverse_adj)
        l.resize(V);
}

void SimpleGraph::setNoVertices(uint32_t n) {
    V = n;
}

uint32_t SimpleGraph::getNoEdges() const {
    uint32_t sum = 0;
    for (const auto &l : adj)
        for(const auto &v : l)
            sum += v.size();
    return sum;
}

uint32_t SimpleGraph::getNoDistinctEdges() const {

    uint32_t sum = 0;

    if(adj_ptr) {
        for (auto sourceVec : *adj_ptr) {

            std::sort(sourceVec.begin(), sourceVec.end());

            uint32_t prevTarget = 0;
            bool first = true;

            for (const auto &target : sourceVec) {
                if (first || prevTarget != target) {
                    first = false;
                    sum++;
                    prevTarget = target;
                }
            }
        }
    }

    for (const auto &labelVec : adj) {
        for (auto sourceVec : labelVec) {

            std::sort(sourceVec.begin(), sourceVec.end());

            uint32_t prevTarget = 0;
            bool first = true;

            for (const auto &target : sourceVec) {
                if (first || prevTarget != target) {
                    first = false;
                    sum++;
                    prevTarget = target;
                }
            }
        }
    }

    return sum;
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

    adj[edgeLabel][from].emplace_back(to);
    reverse_adj[edgeLabel][to].emplace_back(from);
}

void SimpleGraph::addEdges(std::shared_ptr<SimpleGraph> &in, uint32_t projectLabel, bool isInverse) {

    auto _adj = &in->adj[projectLabel];
    auto _reverse_adj = &in->reverse_adj[projectLabel];

    if(!isInverse) {
        adj_ptr = _adj;
        reverse_adj_ptr = _reverse_adj;
    } else {
        adj_ptr = _reverse_adj;
        reverse_adj_ptr = _adj;
    }
}

void SimpleGraph::readFromContiguousFile(const std::string &fileName) {

    std::ifstream graphFile { fileName };

    // Split on spaces.
    std::istream_iterator<std::string> beg(graphFile), end;
    std::vector<std::string> tokens(beg, end); // done!

    // Parse the data.
    std::regex headerPat (R"((\d+),(\d+),(\d+))"); // noNodes,noEdges,noLabels
    std::smatch matches;

    if(tokens.empty()) {
        throw std::runtime_error(std::string("Could not find graph file!"));
    }

    if(std::regex_search(tokens.front(), matches, headerPat)) {
        uint32_t noNodes = (uint32_t) std::stoul(matches[1]);
        uint32_t noLabels = (uint32_t) std::stoul(matches[3]);

        setNoLabels(noLabels);
        setNoVertices(noNodes);
        setDataStructureSizes();
    } else {
        throw std::runtime_error(std::string("Invalid graph header!"));
    }

    for(uint32_t i = 1; i < tokens.size(); i += 4) {
        auto subject = (uint32_t) std::stoul(tokens[i]);
        auto predicate = (uint32_t) std::stoul(tokens[i + 1]);
        auto object = (uint32_t) std::stoul(tokens[i + 2]);

        addEdge(subject, object, predicate);
    }

    graphFile.close();
}