//
// Created by Nikolay Yakovets on 2018-01-31.
//

#include <cmath>
#include "SimpleGraph.h"

void SimpleJoinStorage::addEdges(uint32_t from, std::vector<uint32_t> &data, std::vector<uint32_t>::iterator &end) {
    if(!adj[from]) {
        adj[from] = new std::vector<uint32_t>(static_cast<unsigned long>(std::distance(data.begin(), end)));
    }

    std::copy(data.begin(), end, adj[from]->begin());
    noEdges += adj[from]->size();

    sources[from >> 6] |= SET_BIT(from & 63);
    for(const auto to : *adj[from]) {
        targets[to >> 6] |= SET_BIT(to & 63);
    }
}

void SimpleJoinStorage::setNoVertices(uint32_t n) {
    adj.resize(n);
}

void SimpleJoinStorage::setNoBuckets(uint32_t N) {
    sources.resize(N);
    targets.resize(N);
}

SimpleJoinStorage::SimpleJoinStorage(uint32_t n, uint32_t N) {
    setNoVertices(n);
    setNoBuckets(N);
    this->n = n;
    this->N = N;
}

SimpleJoinStorage::SimpleJoinStorage(SimpleGraphStorage &storage, bool isInverse) {
    noEdges = storage.noEdges;
    n = storage.n;
    N = storage.N;

    if(isInverse) {
        adj_ptr = &storage.reverse_adj;
        noSources = storage.noTargets;
        noTargets = storage.noSources;
        minSource = storage.minTarget;
        maxSource = storage.maxTarget;
    } else {
        adj_ptr = &storage.adj;
        noSources = storage.noSources;
        noTargets = storage.noTargets;
        minSource = storage.minSource;
        maxSource = storage.maxSource;
    }
}

void SimpleJoinStorage::finalize() {
    noSources = countBitsSet(sources);
    noTargets = countBitsSet(targets);

    // Find the minimum and maximum source values.
    minSource = 0;
    for(const uint64_t &s : sources) {
        if(s) {
            minSource += __builtin_ctzll(s);
            break;
        }

        minSource += 64;
    }

    maxSource = static_cast<uint32_t>(sources.size() * 64) - 1;
    for(auto it = sources.rbegin(); it != sources.rend(); ++it) {
        if(*it) {
            maxSource -= __builtin_clzll(*it);
            break;
        }

        maxSource -= 64;
    }

    // Clear unnecessary memory containers.
    std::vector<uint64_t>().swap(targets);
    std::vector<uint64_t>().swap(sources);
}


SimpleJoinStorage::~SimpleJoinStorage() {
    for(const auto &entry : adj) {
        delete(entry);
    }
}

void SimpleGraphStorage::setNoVertices(uint32_t n) {
    adj.resize(n);
    reverse_adj.resize(n);
}

void SimpleGraphStorage::setNoBuckets(uint32_t N) {
    sources.resize(N);
    targets.resize(N);
}

SimpleGraphStorage::SimpleGraphStorage(uint32_t n, uint32_t N) {
    setNoVertices(n);
    setNoBuckets(N);
    this->n = n;
    this->N = N;
}

void SimpleGraphStorage::finalize() {
    noSources = countBitsSet(sources);
    noTargets = countBitsSet(targets);

    // Find the minimum and maximum source values.
    minSource = 0;
    for(const uint64_t &s : sources) {
        if(s) {
            minSource += __builtin_ctzll(s);
            break;
        }

        minSource += 64;
    }

    maxSource = static_cast<uint32_t>(sources.size() * 64) - 1;
    for(auto it = sources.rbegin(); it != sources.rend(); ++it) {
        if(*it) {
            maxSource -= __builtin_clzll(*it);
            break;
        }

        maxSource -= 64;
    }


    // Find the minimum and maximum target values.
    minTarget = 0;
    for(const uint64_t &s : targets) {
        if(s) {
            minTarget += __builtin_ctzll(s);
            break;
        }

        minTarget += 64;
    }

    maxTarget = static_cast<uint32_t>(targets.size() * 64) - 1;
    for(auto it = targets.rbegin(); it != targets.rend(); ++it) {
        if(*it) {
            maxTarget -= __builtin_clzll(*it);
            break;
        }

        maxTarget -= 64;
    }
}

void SimpleGraphStorage::addEdge(uint32_t from, uint32_t to) {
    sources[from >> 6] |= SET_BIT(from & 63);
    targets[to >> 6] |= SET_BIT(to & 63);

    ++noEdges;

    if(!adj[from]) {
        adj[from] = new std::vector<uint32_t>();
    }

    adj[from]->emplace_back(to);

    if(!reverse_adj.empty()) {
        if(!reverse_adj[to]) {
            reverse_adj[to] = new std::vector<uint32_t>();
        }

        reverse_adj[to]->emplace_back(from);
    }
}

SimpleGraphStorage::~SimpleGraphStorage() {
    for(const auto &entry : adj) {
        delete(entry);
    }

    for(const auto &entry : reverse_adj) {
        delete(entry);
    }
}


SimpleGraph::SimpleGraph(uint32_t n)   {
    setNoVertices(n);
}

uint32_t SimpleGraph::getNoVertices() const {
    return V;
}

void SimpleGraph::setDataStructureSizes() {
    E = 0;
    auto N = static_cast<uint32_t>(std::ceil((double) V / 64));
    graphs.resize(L, SimpleGraphStorage{V, N});
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
    SimpleGraphStorage *storage = nullptr;

    for(auto const &tuple : edges) {
        uint32_t source = std::get<1>(tuple);
        uint32_t target = std::get<2>(tuple);
        uint32_t label = std::get<0>(tuple);

        if(first || lastLabel != label) {
            // get the current storage block we are looking at.
            storage = &graphs[label];
        }

        if(first || !(lastSource == source && lastTarget == target && lastLabel == label)) {
            first = false;
            lastSource = source;
            lastTarget = target;
            lastLabel = label;

            storage->addEdge(source, target);
        }
    }

    // For the leaf level, we want the reverse adjacency to be sorted as well.
    for(auto &graph : graphs) {
        for(auto &vertices : graph.reverse_adj) {
            if(vertices) {
                std::sort(vertices->begin(), vertices->end());
            }
        }
        graph.finalize();
    }

    file.close();
}

void SimpleGraph::addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) {
    // We don't use this.
    throw std::runtime_error(std::string("This function is not used."));

}

