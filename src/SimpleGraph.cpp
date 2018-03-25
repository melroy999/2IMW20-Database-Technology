//
// Created by Nikolay Yakovets on 2018-01-31.
//

#include <sstream>
#include "SimpleGraph.h"

SimpleGraph::SimpleGraph(uint32_t n, uint32_t l)   {
    setNoVertices(n);
    setNoLabels(l);
    setDataStructureSizes();
    sources = std::vector<uint64_t>(static_cast<unsigned long>(std::ceil((double) n / 64)));
    targets = std::vector<uint64_t>(static_cast<unsigned long>(std::ceil((double) n / 64)));
}

SimpleGraph::SimpleGraph(const std::shared_ptr<SimpleGraph> &g, uint32_t l) {
    tree = &g->trees[l];
    setNoVertices(g->getNoVertices());
}

uint32_t SimpleGraph::getNoVertices() const {
    return V;
}

void SimpleGraph::setNoVertices(uint32_t n) {
    V = n;
}

uint32_t SimpleGraph::getNoEdges() const {
    return tree ? tree->getNoEdges() : E;
}

uint32_t SimpleGraph::getNoDistinctEdges() const {
    return tree ? tree->getNoEdges() : E;
}

uint32_t SimpleGraph::getNoLabels() const {
    return L;
}

void SimpleGraph::setNoLabels(uint32_t noLabels) {
    L = noLabels;
}

void SimpleGraph::setDataStructureSizes() {
    trees.resize(L, K2Tree(V));
}

void SimpleGraph::addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) {
    if(from >= V || to >= V || edgeLabel >= L)
        throw std::runtime_error(std::string("Edge data out of bounds: ") +
                                         "(" + std::to_string(from) + "," + std::to_string(to) + "," +
                                         std::to_string(edgeLabel) + ")");

    trees[edgeLabel].addEdge(from, to);
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
    setDataStructureSizes();

    for(int i = 0; i < noEdges; i++) {
        uint32_t subject, predicate, object;

        file >> subject >> predicate >> object;
        file.ignore(INT32_MAX, '\n');

        addEdge(subject, object, predicate);
    }

    sources = std::vector<uint64_t>(static_cast<unsigned long>(std::ceil((double) noNodes / 64)));
    targets = std::vector<uint64_t>(static_cast<unsigned long>(std::ceil((double) noNodes / 64)));
    finalizeTree();

    file.close();
}

void SimpleGraph::finalizeTree() {
    // Calculate the number of edges, sources and target vertices.
    E = 0;
    for(auto &tree : trees) {
        tree.finalizeTree();
        E += tree.getNoEdges();

        doOrIP(&sources, &tree.getSources());
        doOrIP(&targets, &tree.getTargets());
    }

    S = popcountVector(&sources);
    T = popcountVector(&targets);
}

uint32_t SimpleGraph::getNoSources() const {
    if(inverse) {
        return tree ? tree->getNoTargets() : T;
    } else {
        return tree ? tree->getNoSources() : S;
    }
}

uint32_t SimpleGraph::getNoTargets() const {
    if(inverse) {
        return tree ? tree->getNoSources() : S;
    } else {
        return tree ? tree->getNoTargets() : T;
    }
}

const std::vector<uint64_t> &SimpleGraph::getSources() const {
    return tree ? tree->getSources() : sources;
}

const std::vector<uint64_t> &SimpleGraph::getTargets() const {
    return tree ? tree->getTargets() : targets;
}


