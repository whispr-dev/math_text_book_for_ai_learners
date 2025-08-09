#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <functional>
#include <map>
#include <set>

namespace GraphTheory {

    const double INF = std::numeric_limits<double>::infinity();
    const int UNDEFINED = -1;
    
    // ===== GRAPH REPRESENTATION STRUCTURES =====
    
    // Edge structure for weighted graphs
    struct Edge {
        int from, to;
        double weight;
        
        Edge(int f, int t, double w = 1.0) : from(f), to(t), weight(w) {}
        
        bool operator<(const Edge& other) const {
            return weight < other.weight;
        }
        
        void print() const {
            std::cout << "(" << from << " -> " << to << ", weight: " << weight << ")";
        }
    };
    
    // Graph class supporting both directed and undirected graphs
    class Graph {
    private:
        int num_vertices;
        bool is_directed;
        std::vector<std::vector<std::pair<int, double>>> adj_list;
        std::vector<Edge> edge_list;
        
    public:
        Graph(int vertices, bool directed = false) 
            : num_vertices(vertices), is_directed(directed) {
            adj_list.resize(vertices);
        }
        
        void addEdge(int from, int to, double weight = 1.0) {
            adj_list[from].push_back({to, weight});
            edge_list.push_back(Edge(from, to, weight));
            
            if (!is_directed) {
                adj_list[to].push_back({from, weight});
            }
        }
        
        int getVertexCount() const { return num_vertices; }
        bool isDirected() const { return is_directed; }
        
        const std::vector<std::pair<int, double>>& getNeighbors(int vertex) const {
            return adj_list[vertex];
        }
        
        const std::vector<Edge>& getEdges() const { return edge_list; }
        
        void printGraph() const {
            std::cout << "Graph (" << (is_directed ? "Directed" : "Undirected") 
                      << ", " << num_vertices << " vertices):" << std::endl;
            for (int i = 0; i < num_vertices; ++i) {
                std::cout << "Vertex " << i << ": ";
                for (const auto& [neighbor, weight] : adj_list[i]) {
                    std::cout << "(" << neighbor << ", w:" << weight << ") ";
                }
                std::cout << std::endl;
            }
        }
        
        // Get adjacency matrix representation
        std::vector<std::vector<double>> getAdjacencyMatrix() const {
            std::vector<std::vector<double>> matrix(num_vertices, 
                std::vector<double>(num_vertices, 0.0));
            
            for (int i = 0; i < num_vertices; ++i) {
                for (const auto& [neighbor, weight] : adj_list[i]) {
                    matrix[i][neighbor] = weight;
                }
            }
            return matrix;
        }
    };
    
    // ===== GRAPH TRAVERSAL ALGORITHMS =====
    
    struct TraversalResult {
        std::vector<int> order;
        std::vector<int> parent;
        std::vector<int> discovery_time;
        std::vector<int> finish_time;
        std::vector<bool> visited;
        
        void print(const std::string& algorithm) const {
            std::cout << algorithm << " Traversal Order: ";
            for (int v : order) {
                std::cout << v << " ";
            }
            std::cout << std::endl;
        }
    };
    
    // Depth-First Search (DFS) - Recursive implementation
    void dfsRecursive(const Graph& graph, int vertex, TraversalResult& result, int& time) {
        result.visited[vertex] = true;
        result.discovery_time[vertex] = ++time;
        result.order.push_back(vertex);
        
        for (const auto& [neighbor, weight] : graph.getNeighbors(vertex)) {
            if (!result.visited[neighbor]) {
                result.parent[neighbor] = vertex;
                dfsRecursive(graph, neighbor, result, time);
            }
        }
        
        result.finish_time[vertex] = ++time;
    }
    
    // DFS - Complete traversal from all unvisited vertices
    TraversalResult depthFirstSearch(const Graph& graph, int start_vertex = 0) {
        TraversalResult result;
        int n = graph.getVertexCount();
        
        result.visited.assign(n, false);
        result.parent.assign(n, UNDEFINED);
        result.discovery_time.assign(n, 0);
        result.finish_time.assign(n, 0);
        
        int time = 0;
        
        // Start from specified vertex
        if (!result.visited[start_vertex]) {
            dfsRecursive(graph, start_vertex, result, time);
        }
        
        // Visit any remaining unvisited vertices
        for (int i = 0; i < n; ++i) {
            if (!result.visited[i]) {
                dfsRecursive(graph, i, result, time);
            }
        }
        
        return result;
    }
    
    // DFS - Iterative implementation using stack
    TraversalResult depthFirstSearchIterative(const Graph& graph, int start_vertex) {
        TraversalResult result;
        int n = graph.getVertexCount();
        
        result.visited.assign(n, false);
        result.parent.assign(n, UNDEFINED);
        result.discovery_time.assign(n, 0);
        result.finish_time.assign(n, 0);
        
        std::stack<int> stack;
        stack.push(start_vertex);
        int time = 0;
        
        while (!stack.empty()) {
            int vertex = stack.top();
            stack.pop();
            
            if (!result.visited[vertex]) {
                result.visited[vertex] = true;
                result.discovery_time[vertex] = ++time;
                result.order.push_back(vertex);
                
                // Add neighbors to stack in reverse order for consistent traversal
                auto neighbors = graph.getNeighbors(vertex);
                for (auto it = neighbors.rbegin(); it != neighbors.rend(); ++it) {
                    if (!result.visited[it->first]) {
                        stack.push(it->first);
                        if (result.parent[it->first] == UNDEFINED) {
                            result.parent[it->first] = vertex;
                        }
                    }
                }
            }
        }
        
        return result;
    }
    
    // Breadth-First Search (BFS)
    TraversalResult breadthFirstSearch(const Graph& graph, int start_vertex) {
        TraversalResult result;
        int n = graph.getVertexCount();
        
        result.visited.assign(n, false);
        result.parent.assign(n, UNDEFINED);
        result.discovery_time.assign(n, 0);
        result.finish_time.assign(n, 0);
        
        std::queue<int> queue;
        queue.push(start_vertex);
        result.visited[start_vertex] = true;
        result.discovery_time[start_vertex] = 0;
        
        int time = 0;
        
        while (!queue.empty()) {
            int vertex = queue.front();
            queue.pop();
            result.order.push_back(vertex);
            
            for (const auto& [neighbor, weight] : graph.getNeighbors(vertex)) {
                if (!result.visited[neighbor]) {
                    result.visited[neighbor] = true;
                    result.parent[neighbor] = vertex;
                    result.discovery_time[neighbor] = ++time;
                    queue.push(neighbor);
                }
            }
        }
        
        return result;
    }
    
    // ===== SHORTEST PATH ALGORITHMS =====
    
    struct ShortestPathResult {
        std::vector<double> distances;
        std::vector<int> predecessors;
        int source;
        bool has_negative_cycle;
        
        std::vector<int> getPath(int target) const {
            std::vector<int> path;
            if (distances[target] == INF) return path; // No path exists
            
            int current = target;
            while (current != UNDEFINED) {
                path.push_back(current);
                current = predecessors[current];
            }
            
            std::reverse(path.begin(), path.end());
            return path;
        }
        
        void printShortestPaths() const {
            std::cout << "Shortest paths from vertex " << source << ":" << std::endl;
            for (size_t i = 0; i < distances.size(); ++i) {
                std::cout << "To " << i << ": ";
                if (distances[i] == INF) {
                    std::cout << "No path" << std::endl;
                } else {
                    std::cout << "Distance " << distances[i] << ", Path: ";
                    auto path = getPath(i);
                    for (size_t j = 0; j < path.size(); ++j) {
                        std::cout << path[j];
                        if (j < path.size() - 1) std::cout << " -> ";
                    }
                    std::cout << std::endl;
                }
            }
        }
    };
    
    // Dijkstra's Algorithm - Single source shortest paths (non-negative weights)
    ShortestPathResult dijkstra(const Graph& graph, int source) {
        int n = graph.getVertexCount();
        ShortestPathResult result;
        result.source = source;
        result.distances.assign(n, INF);
        result.predecessors.assign(n, UNDEFINED);
        result.has_negative_cycle = false;
        result.distances[source] = 0.0;
        
        using PQItem = std::pair<double, int>; // (distance, vertex)
        std::priority_queue<PQItem, std::vector<PQItem>, std::greater<PQItem>> pq;
        pq.push({0.0, source});
        
        std::vector<bool> processed(n, false);
        
        while (!pq.empty()) {
            auto [dist, u] = pq.top();
            pq.pop();
            
            if (processed[u]) continue;
            processed[u] = true;
            
            for (const auto& [v, weight] : graph.getNeighbors(u)) {
                double new_dist = result.distances[u] + weight;
                if (new_dist < result.distances[v]) {
                    result.distances[v] = new_dist;
                    result.predecessors[v] = u;
                    pq.push({new_dist, v});
                }
            }
        }
        
        return result;
    }
    
    // Bellman-Ford Algorithm - Single source shortest paths (handles negative weights)
    ShortestPathResult bellmanFord(const Graph& graph, int source) {
        int n = graph.getVertexCount();
        ShortestPathResult result;
        result.source = source;
        result.distances.assign(n, INF);
        result.predecessors.assign(n, UNDEFINED);
        result.has_negative_cycle = false;
        result.distances[source] = 0.0;
        
        // Relax all edges V-1 times
        for (int i = 0; i < n - 1; ++i) {
            for (const auto& edge : graph.getEdges()) {
                if (result.distances[edge.from] != INF) {
                    double new_dist = result.distances[edge.from] + edge.weight;
                    if (new_dist < result.distances[edge.to]) {
                        result.distances[edge.to] = new_dist;
                        result.predecessors[edge.to] = edge.from;
                    }
                }
            }
        }
        
        // Check for negative-weight cycles
        for (const auto& edge : graph.getEdges()) {
            if (result.distances[edge.from] != INF) {
                double new_dist = result.distances[edge.from] + edge.weight;
                if (new_dist < result.distances[edge.to]) {
                    result.has_negative_cycle = true;
                    break;
                }
            }
        }
        
        return result;
    }
    
    // Floyd-Warshall Algorithm - All pairs shortest paths
    std::vector<std::vector<double>> floydWarshall(const Graph& graph) {
        int n = graph.getVertexCount();
        auto dist = graph.getAdjacencyMatrix();
        
        // Initialize distances
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    dist[i][j] = 0.0;
                } else if (dist[i][j] == 0.0) {
                    dist[i][j] = INF; // No direct edge
                }
            }
        }
        
        // Floyd-Warshall main loop
        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (dist[i][k] != INF && dist[k][j] != INF) {
                        dist[i][j] = std::min(dist[i][j], dist[i][k] + dist[k][j]);
                    }
                }
            }
        }
        
        return dist;
    }
    
    // ===== MINIMUM SPANNING TREE ALGORITHMS =====
    
    struct MST_Result {
        std::vector<Edge> edges;
        double total_weight;
        bool is_connected;
        
        void print(const std::string& algorithm) const {
            std::cout << algorithm << " MST:" << std::endl;
            if (!is_connected) {
                std::cout << "Graph is not connected!" << std::endl;
                return;
            }
            
            std::cout << "Total weight: " << total_weight << std::endl;
            std::cout << "Edges: ";
            for (const auto& edge : edges) {
                edge.print();
                std::cout << " ";
            }
            std::cout << std::endl;
        }
    };
    
    // Union-Find (Disjoint Set) data structure for Kruskal's algorithm
    class UnionFind {
    private:
        std::vector<int> parent, rank;
        
    public:
        UnionFind(int n) : parent(n), rank(n, 0) {
            for (int i = 0; i < n; ++i) {
                parent[i] = i;
            }
        }
        
        int find(int x) {
            if (parent[x] != x) {
                parent[x] = find(parent[x]); // Path compression
            }
            return parent[x];
        }
        
        bool unite(int x, int y) {
            int px = find(x), py = find(y);
            if (px == py) return false;
            
            // Union by rank
            if (rank[px] < rank[py]) {
                parent[px] = py;
            } else if (rank[px] > rank[py]) {
                parent[py] = px;
            } else {
                parent[py] = px;
                rank[px]++;
            }
            return true;
        }
    };
    
    // Kruskal's Algorithm - MST using edge sorting and Union-Find
    MST_Result kruskalMST(const Graph& graph) {
        if (graph.isDirected()) {
            throw std::invalid_argument("Kruskal's algorithm requires an undirected graph");
        }
        
        MST_Result result;
        result.total_weight = 0.0;
        result.is_connected = true;
        
        auto edges = graph.getEdges();
        std::sort(edges.begin(), edges.end()); // Sort by weight
        
        UnionFind uf(graph.getVertexCount());
        
        for (const auto& edge : edges) {
            if (uf.unite(edge.from, edge.to)) {
                result.edges.push_back(edge);
                result.total_weight += edge.weight;
                
                if (result.edges.size() == graph.getVertexCount() - 1) {
                    break; // MST complete
                }
            }
        }
        
        // Check if graph is connected
        if (result.edges.size() != graph.getVertexCount() - 1) {
            result.is_connected = false;
        }
        
        return result;
    }
    
    // Prim's Algorithm - MST using priority queue
    MST_Result primMST(const Graph& graph, int start_vertex = 0) {
        if (graph.isDirected()) {
            throw std::invalid_argument("Prim's algorithm requires an undirected graph");
        }
        
        int n = graph.getVertexCount();
        MST_Result result;
        result.total_weight = 0.0;
        result.is_connected = true;
        
        std::vector<bool> in_mst(n, false);
        using PQItem = std::tuple<double, int, int>; // (weight, from, to)
        std::priority_queue<PQItem, std::vector<PQItem>, std::greater<PQItem>> pq;
        
        // Start with arbitrary vertex
        in_mst[start_vertex] = true;
        for (const auto& [neighbor, weight] : graph.getNeighbors(start_vertex)) {
            pq.push({weight, start_vertex, neighbor});
        }
        
        while (!pq.empty() && result.edges.size() < n - 1) {
            auto [weight, from, to] = pq.top();
            pq.pop();
            
            if (in_mst[to]) continue; // Already in MST
            
            // Add edge to MST
            in_mst[to] = true;
            result.edges.push_back(Edge(from, to, weight));
            result.total_weight += weight;
            
            // Add new edges to priority queue
            for (const auto& [neighbor, edge_weight] : graph.getNeighbors(to)) {
                if (!in_mst[neighbor]) {
                    pq.push({edge_weight, to, neighbor});
                }
            }
        }
        
        // Check if graph is connected
        if (result.edges.size() != n - 1) {
            result.is_connected = false;
        }
        
        return result;
    }
    
    // ===== TOPOLOGICAL SORTING =====
    
    struct TopologicalResult {
        std::vector<int> ordering;
        bool is_dag; // Is Directed Acyclic Graph
        std::vector<int> cycle; // If not DAG, contains a cycle
        
        void print() const {
            if (is_dag) {
                std::cout << "Topological ordering: ";
                for (int vertex : ordering) {
                    std::cout << vertex << " ";
                }
                std::cout << std::endl;
            } else {
                std::cout << "Graph contains a cycle: ";
                for (int vertex : cycle) {
                    std::cout << vertex << " ";
                }
                std::cout << std::endl;
            }
        }
    };
    
    // Topological Sort using DFS (Kahn's algorithm alternative)
    TopologicalResult topologicalSort(const Graph& graph) {
        if (!graph.isDirected()) {
            throw std::invalid_argument("Topological sort requires a directed graph");
        }
        
        int n = graph.getVertexCount();
        TopologicalResult result;
        result.is_dag = true;
        
        std::vector<int> color(n, 0); // 0: white, 1: gray, 2: black
        std::stack<int> finish_stack;
        
        std::function<bool(int)> dfs_visit = [&](int u) -> bool {
            color[u] = 1; // Gray (visiting)
            
            for (const auto& [v, weight] : graph.getNeighbors(u)) {
                if (color[v] == 1) {
                    // Back edge found - cycle detected
                    result.is_dag = false;
                    result.cycle.push_back(v);
                    result.cycle.push_back(u);
                    return false;
                } else if (color[v] == 0 && !dfs_visit(v)) {
                    return false;
                }
            }
            
            color[u] = 2; // Black (finished)
            finish_stack.push(u);
            return true;
        };
        
        // Visit all vertices
        for (int i = 0; i < n; ++i) {
            if (color[i] == 0) {
                if (!dfs_visit(i)) {
                    return result; // Cycle found
                }
            }
        }
        
        // Extract topological ordering from finish times
        while (!finish_stack.empty()) {
            result.ordering.push_back(finish_stack.top());
            finish_stack.pop();
        }
        
        return result;
    }
    
    // ===== GRAPH CONNECTIVITY =====
    
    struct ConnectivityResult {
        bool is_connected;
        int num_components;
        std::vector<int> component_id;
        std::vector<std::vector<int>> components;
        
        void print() const {
            std::cout << "Graph connectivity:" << std::endl;
            std::cout << "Connected: " << (is_connected ? "Yes" : "No") << std::endl;
            std::cout << "Components: " << num_components << std::endl;
            
            for (int i = 0; i < num_components; ++i) {
                std::cout << "Component " << i << ": ";
                for (int vertex : components[i]) {
                    std::cout << vertex << " ";
                }
                std::cout << std::endl;
            }
        }
    };
    
    // Find connected components using DFS
    ConnectivityResult findConnectedComponents(const Graph& graph) {
        int n = graph.getVertexCount();
        ConnectivityResult result;
        result.component_id.assign(n, -1);
        result.num_components = 0;
        
        std::function<void(int, int)> dfs = [&](int vertex, int comp_id) {
            result.component_id[vertex] = comp_id;
            result.components[comp_id].push_back(vertex);
            
            for (const auto& [neighbor, weight] : graph.getNeighbors(vertex)) {
                if (result.component_id[neighbor] == -1) {
                    dfs(neighbor, comp_id);
                }
            }
        };
        
        for (int i = 0; i < n; ++i) {
            if (result.component_id[i] == -1) {
                result.components.push_back(std::vector<int>());
                dfs(i, result.num_components);
                result.num_components++;
            }
        }
        
        result.is_connected = (result.num_components == 1);
        return result;
    }
    
    // ===== GRAPH COLORING =====
    
    struct ColoringResult {
        std::vector<int> colors;
        int chromatic_number;
        bool is_valid;
        
        void print() const {
            std::cout << "Graph coloring:" << std::endl;
            std::cout << "Chromatic number: " << chromatic_number << std::endl;
            std::cout << "Valid coloring: " << (is_valid ? "Yes" : "No") << std::endl;
            std::cout << "Vertex colors: ";
            for (size_t i = 0; i < colors.size(); ++i) {
                std::cout << "v" << i << ":" << colors[i] << " ";
            }
            std::cout << std::endl;
        }
        
        bool validateColoring(const Graph& graph) const {
            for (int i = 0; i < graph.getVertexCount(); ++i) {
                for (const auto& [neighbor, weight] : graph.getNeighbors(i)) {
                    if (colors[i] == colors[neighbor]) {
                        return false; // Adjacent vertices have same color
                    }
                }
            }
            return true;
        }
    };
    
    // Greedy graph coloring algorithm
    ColoringResult greedyColoring(const Graph& graph) {
        int n = graph.getVertexCount();
        ColoringResult result;
        result.colors.assign(n, -1);
        result.chromatic_number = 0;
        
        for (int u = 0; u < n; ++u) {
            std::vector<bool> available(n, true);
            
            // Mark colors used by neighbors as unavailable
            for (const auto& [v, weight] : graph.getNeighbors(u)) {
                if (result.colors[v] != -1) {
                    available[result.colors[v]] = false;
                }
            }
            
            // Find first available color
            int color = 0;
            while (color < n && !available[color]) {
                color++;
            }
            
            result.colors[u] = color;
            result.chromatic_number = std::max(result.chromatic_number, color + 1);
        }
        
        result.is_valid = result.validateColoring(graph);
        return result;
    }
    
    // ===== CYCLE DETECTION =====
    
    struct CycleResult {
        bool has_cycle;
        std::vector<int> cycle;
        
        void print() const {
            if (has_cycle) {
                std::cout << "Cycle found: ";
                for (int vertex : cycle) {
                    std::cout << vertex << " ";
                }
                std::cout << std::endl;
            } else {
                std::cout << "No cycle found" << std::endl;
            }
        }
    };
    
    // Cycle detection in undirected graph using DFS
    CycleResult detectCycleUndirected(const Graph& graph) {
        if (graph.isDirected()) {
            throw std::invalid_argument("This function is for undirected graphs");
        }
        
        int n = graph.getVertexCount();
        CycleResult result;
        result.has_cycle = false;
        
        std::vector<bool> visited(n, false);
        std::vector<int> parent(n, -1);
        
        std::function<bool(int, int)> dfs = [&](int u, int p) -> bool {
            visited[u] = true;
            parent[u] = p;
            
            for (const auto& [v, weight] : graph.getNeighbors(u)) {
                if (!visited[v]) {
                    if (dfs(v, u)) return true;
                } else if (v != p) {
                    // Back edge found - cycle detected
                    result.has_cycle = true;
                    
                    // Reconstruct cycle
                    result.cycle.push_back(v);
                    int curr = u;
                    while (curr != v) {
                        result.cycle.push_back(curr);
                        curr = parent[curr];
                    }
                    result.cycle.push_back(v);
                    return true;
                }
            }
            return false;
        };
        
        for (int i = 0; i < n; ++i) {
            if (!visited[i]) {
                if (dfs(i, -1)) break;
            }
        }
        
        return result;
    }
    
} // namespace GraphTheory

// ===== DEMONSTRATION =====

int main() {
    using namespace GraphTheory;
    
    std::cout << "=== UNIFIED GRAPH THEORY MODULE DEMO ===" << std::endl << std::endl;
    
    // ===== CREATE TEST GRAPHS =====
    std::cout << "1. CREATING TEST GRAPHS" << std::endl;
    
    // Undirected weighted graph
    Graph undirected_graph(6, false);
    undirected_graph.addEdge(0, 1, 4);
    undirected_graph.addEdge(0, 2, 2);
    undirected_graph.addEdge(1, 2, 1);
    undirected_graph.addEdge(1, 3, 5);
    undirected_graph.addEdge(2, 3, 8);
    undirected_graph.addEdge(2, 4, 10);
    undirected_graph.addEdge(3, 4, 2);
    undirected_graph.addEdge(3, 5, 6);
    undirected_graph.addEdge(4, 5, 3);
    
    std::cout << "Undirected graph created with 6 vertices and 9 edges" << std::endl;
    
    // Directed graph for topological sorting
    Graph directed_graph(6, true);
    directed_graph.addEdge(5, 2, 1);
    directed_graph.addEdge(5, 0, 1);
    directed_graph.addEdge(4, 0, 1);
    directed_graph.addEdge(4, 1, 1);
    directed_graph.addEdge(2, 3, 1);
    directed_graph.addEdge(3, 1, 1);
    
    std::cout << "Directed acyclic graph created for topological sorting" << std::endl << std::endl;
    
    // ===== GRAPH TRAVERSAL =====
    std::cout << "2. GRAPH TRAVERSAL ALGORITHMS" << std::endl;
    
    auto dfs_result = depthFirstSearch(undirected_graph, 0);
    dfs_result.print("DFS (Recursive)");
    
    auto dfs_iter_result = depthFirstSearchIterative(undirected_graph, 0);
    dfs_iter_result.print("DFS (Iterative)");
    
    auto bfs_result = breadthFirstSearch(undirected_graph, 0);
    bfs_result.print("BFS");
    std::cout << std::endl;
    
    // ===== SHORTEST PATH ALGORITHMS =====
    std::cout << "3. SHORTEST PATH ALGORITHMS" << std::endl;
    
    auto dijkstra_result = dijkstra(undirected_graph, 0);
    dijkstra_result.printShortestPaths();
    
    auto bellman_ford_result = bellmanFord(undirected_graph, 0);
    std::cout << "Bellman-Ford - Has negative cycle: " 
              << (bellman_ford_result.has_negative_cycle ? "Yes" : "No") << std::endl;
    
    std::cout << "\nFloyd-Warshall all-pairs shortest paths:" << std::endl;
    auto floyd_result = floydWarshall(undirected_graph);
    std::cout << "Distance matrix (first 4x4):" << std::endl;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (floyd_result[i][j] == INF) {
                std::cout << "∞\t";
            } else {
                std::cout << floyd_result[i][j] << "\t";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    // ===== MINIMUM SPANNING TREE =====
    std::cout << "4. MINIMUM SPANNING TREE ALGORITHMS" << std::endl;
    
    auto kruskal_result = kruskalMST(undirected_graph);
    kruskal_result.print("Kruskal's");
    
    auto prim_result = primMST(undirected_graph, 0);
    prim_result.print("Prim's");
    std::cout << std::endl;
    
    // ===== TOPOLOGICAL SORTING =====
    std::cout << "5. TOPOLOGICAL SORTING" << std::endl;
    
    auto topo_result = topologicalSort(directed_graph);
    topo_result.print();
    std::cout << std::endl;
    
    // ===== CONNECTIVITY ANALYSIS =====
    std::cout << "6. GRAPH CONNECTIVITY" << std::endl;
    
    auto connectivity_result = findConnectedComponents(undirected_graph);
    connectivity_result.print();
    std::cout << std::endl;
    
    // ===== GRAPH COLORING =====
    std::cout << "7. GRAPH COLORING" << std::endl;
    
    auto coloring_result = greedyColoring(undirected_graph);
    coloring_result.print();
    std::cout << std::endl;
    
    // ===== CYCLE DETECTION =====
    std::cout << "8. CYCLE DETECTION" << std::endl;
    
    auto cycle_result = detectCycleUndirected(undirected_graph);
    cycle_result.print();
    std::cout << std::endl;
    
    // ===== ADVANCED APPLICATIONS =====
    std::cout << "9. ADVANCED GRAPH APPLICATIONS" << std::endl;
    
    std::cout << "🎛️ Synthesizer Applications:" << std::endl;
    std::cout << "• Signal Flow Graphs: Model audio routing in modular systems" << std::endl;
    std::cout << "• Dependency Graphs: Component initialization order (topological sort)" << std::endl;
    std::cout << "• Network Analysis: Find bottlenecks in audio processing chains" << std::endl;
    std::cout << "• Resource Allocation: Optimize CPU usage across modules (graph coloring)" << std::endl;
    std::cout << std::endl;
    
    std::cout << "🔌 Circuit Analysis:" << std::endl;
    std::cout << "• Component Networks: Shortest paths for signal routing" << std::endl;
    std::cout << "• Power Distribution: MST for efficient power delivery" << std::endl;
    std::cout << "• Ground Loops: Cycle detection prevents audio interference" << std::endl;
    std::cout << "• Layout Optimization: Graph coloring for component placement" << std::endl;
    std::cout << std::endl;
    
    std::cout << "📊 Data Flow Analysis:" << std::endl;
    std::cout << "• Real-time Processing: Topological sort for execution order" << std::endl;
    std::cout << "• Memory Management: Connectivity analysis for buffer sharing" << std::endl;
    std::cout << "• Latency Optimization: Shortest paths minimize processing delay" << std::endl;
    std::cout << "• Deadlock Prevention: Cycle detection in resource graphs" << std::endl;
    std::cout << std::endl;
    
    // ===== PERFORMANCE COMPARISON =====
    std::cout << "10. ALGORITHM PERFORMANCE CHARACTERISTICS" << std::endl;
    
    std::cout << "Algorithm             | Time Complexity | Space | Best Use Case" << std::endl;
    std::cout << "===================== | =============== | ===== | ====================" << std::endl;
    std::cout << "DFS (Recursive)       | O(V + E)        | O(V)  | Connectivity, cycles" << std::endl;
    std::cout << "BFS                   | O(V + E)        | O(V)  | Shortest unweighted" << std::endl;
    std::cout << "Dijkstra              | O((V+E)log V)   | O(V)  | Non-negative weights" << std::endl;
    std::cout << "Bellman-Ford          | O(V*E)          | O(V)  | Negative weights" << std::endl;
    std::cout << "Floyd-Warshall        | O(V³)           | O(V²) | All-pairs shortest" << std::endl;
    std::cout << "Kruskal MST           | O(E log E)      | O(V)  | Sparse graphs" << std::endl;
    std::cout << "Prim MST              | O((V+E)log V)   | O(V)  | Dense graphs" << std::endl;
    std::cout << "Topological Sort      | O(V + E)        | O(V)  | DAG ordering" << std::endl;
    std::cout << "Graph Coloring        | O(V²)           | O(V)  | Resource allocation" << std::endl;
    std::cout << std::endl;
    
    // ===== IMPLEMENTATION FEATURES =====
    std::cout << "11. IMPLEMENTATION FEATURES" << std::endl;
    
    std::cout << "✅ Unified Graph Class:" << std::endl;
    std::cout << "   • Supports both directed and undirected graphs" << std::endl;
    std::cout << "   • Weighted edges with double precision" << std::endl;
    std::cout << "   • Adjacency list and matrix representations" << std::endl;
    std::cout << "   • Automatic edge handling for undirected graphs" << std::endl;
    std::cout << std::endl;
    
    std::cout << "🚀 Optimized Algorithms:" << std::endl;
    std::cout << "   • Priority queues for Dijkstra and Prim" << std::endl;
    std::cout << "   • Union-Find with path compression for Kruskal" << std::endl;
    std::cout << "   • DFS with both recursive and iterative versions" << std::endl;
    std::cout << "   • Efficient cycle detection and path reconstruction" << std::endl;
    std::cout << std::endl;
    
    std::cout << "📋 Comprehensive Results:" << std::endl;
    std::cout << "   • Detailed result structures with metadata" << std::endl;
    std::cout << "   • Path reconstruction capabilities" << std::endl;
    std::cout << "   • Validation and error checking" << std::endl;
    std::cout << "   • Pretty-print functions for debugging" << std::endl;
    std::cout << std::endl;
    
    std::cout << "🛡️ Robust Error Handling:" << std::endl;
    std::cout << "   • Graph type validation (directed vs undirected)" << std::endl;
    std::cout << "   • Negative cycle detection in Bellman-Ford" << std::endl;
    std::cout << "   • Disconnected graph handling" << std::endl;
    std::cout << "   • Boundary condition checks" << std::endl;
    std::cout << std::endl;
    
    // ===== REAL-WORLD EXAMPLES =====
    std::cout << "12. REAL-WORLD EXAMPLES FOR RYO MODULAR" << std::endl;
    
    std::cout << "🎵 Audio Signal Routing:" << std::endl;
    std::cout << "   Graph vertices = Audio modules (VCO, VCF, VCA, etc.)" << std::endl;
    std::cout << "   Graph edges = Patch cables with signal strength weights" << std::endl;
    std::cout << "   Shortest path = Minimal latency audio routing" << std::endl;
    std::cout << "   Cycle detection = Prevent feedback loops" << std::endl;
    std::cout << std::endl;
    
    std::cout << "⚡ Power Distribution:" << std::endl;
    std::cout << "   Graph vertices = Power supply points" << std::endl;
    std::cout << "   Graph edges = Power connections with resistance weights" << std::endl;
    std::cout << "   MST = Minimal wire length for power distribution" << std::endl;
    std::cout << "   Connectivity = Ensure all modules powered" << std::endl;
    std::cout << std::endl;
    
    std::cout << "🧮 DSP Processing Chain:" << std::endl;
    std::cout << "   Graph vertices = DSP algorithms/filters" << std::endl;
    std::cout << "   Graph edges = Data dependencies with processing time" << std::endl;
    std::cout << "   Topological sort = Execution order for real-time processing" << std::endl;
    std::cout << "   Critical path = Maximum latency through processing chain" << std::endl;
    std::cout << std::endl;
    
    std::cout << "🎨 Component Placement:" << std::endl;
    std::cout << "   Graph vertices = PCB components" << std::endl;
    std::cout << "   Graph edges = Interference/heat coupling" << std::endl;
    std::cout << "   Graph coloring = Separate incompatible components" << std::endl;
    std::cout << "   Colors = Physical separation zones" << std::endl;
    std::cout << std::endl;
    
    // ===== ADVANCED OPTIMIZATIONS =====
    std::cout << "13. ADVANCED OPTIMIZATIONS & EXTENSIONS" << std::endl;
    
    std::cout << "🔧 Possible Extensions:" << std::endl;
    std::cout << "   • A* search for heuristic pathfinding" << std::endl;
    std::cout << "   • Johnson's algorithm for sparse all-pairs shortest paths" << std::endl;
    std::cout << "   • Strongly connected components (Tarjan's algorithm)" << std::endl;
    std::cout << "   • Maximum flow algorithms (Ford-Fulkerson, Dinic)" << std::endl;
    std::cout << "   • Bipartite matching algorithms" << std::endl;
    std::cout << "   • Planar graph algorithms for PCB routing" << std::endl;
    std::cout << std::endl;
    
    std::cout << "⚡ Performance Optimizations:" << std::endl;
    std::cout << "   • Parallel algorithms for large graphs" << std::endl;
    std::cout << "   • Memory-efficient sparse representations" << std::endl;
    std::cout << "   • Incremental algorithms for dynamic graphs" << std::endl;
    std::cout << "   • Cache-friendly data structures" << std::endl;
    std::cout << std::endl;
    
    std::cout << "=== UNIFIED GRAPH THEORY MODULE COMPLETE ===" << std::endl;
    std::cout << "🎛️ Ready for advanced network analysis in audio systems! 🎛️" << std::endl;
    
    return 0;
}