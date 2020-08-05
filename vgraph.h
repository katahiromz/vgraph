#include <vector>
#include <stack>
#include <queue>
#include <cstring>
#include <limits>
#include <cassert>
#include <iostream>

namespace vgraph
{
    //////////////////////////////////////////////////////////////////////////
    // Vertex

    typedef std::size_t size_t;
    typedef size_t vertex_size_t;

    inline size_t identity(size_t value)
    {
        return value;
    }

    //////////////////////////////////////////////////////////////////////////
    // Cost

    typedef int cost_t;
    static const cost_t max_cost = std::numeric_limits<cost_t>::max() / 4;

    template <typename T>
    inline cost_t cost_of_zero(const T&)
    {
        return 0;
    }

    inline cost_t cost_of_zero(vertex_size_t i, vertex_size_t j)
    {
        return 0;
    }

    template <typename T>
    inline cost_t cost_of_one(const T&)
    {
        return 1;
    }

    inline cost_t cost_of_one(vertex_size_t i, vertex_size_t j)
    {
        return 1;
    }

    template <typename T>
    inline cost_t cost_of_max(const T&)
    {
        return max_cost;
    }

    inline cost_t cost_of_max(vertex_size_t i, vertex_size_t j)
    {
        return max_cost;
    }

    //////////////////////////////////////////////////////////////////////////
    // Edge

    typedef size_t edge_size_t;

    struct edge
    {
        edge_size_t from, to;
    };

    struct cost_edge
    {
        edge_size_t from, to;
        cost_t cost = 1;
    };

    struct from_centric_edge
    {
        edge_size_t to;
    };

    struct from_centric_cost_edge
    {
        edge_size_t to;
        cost_t cost = 1;
    };

    //////////////////////////////////////////////////////////////////////////
    // Cost 2

    inline cost_t cost_of(const edge& e)
    {
        return 1;
    }

    inline cost_t cost_of(const cost_edge& e)
    {
        return e.cost;
    }

    inline cost_t cost_of(const from_centric_edge& e)
    {
        return 1;
    }

    inline cost_t cost_of(const from_centric_cost_edge& e)
    {
        return e.cost;
    }

    //////////////////////////////////////////////////////////////////////////
    // Count

    template <typename T, size_t t_num>
    inline size_t count_of(const T (&t)[t_num])
    {
        return t_num;
    }

    template <typename T>
    inline size_t count_of(const std::vector<T>& values)
    {
        return values.size();
    }

    //////////////////////////////////////////////////////////////////////////
    // Graph

    template <size_t t_size>
    struct array_graph
    {
        typedef size_t array_t[t_size][t_size];
        array_t m_array;

        array_graph()
        {
        }

        array_graph(const array_t& a)
        {
            std::memcpy(&m_array, a, sizeof(a));
        }

        size_t adjacent(vertex_size_t from, vertex_size_t to) const
        {
            assert(from < t_size);
            assert(to < t_size);
            return m_array[from][to];
        }

        void add_edge(vertex_size_t from, vertex_size_t to, cost_t cost = 1)
        {
            assert(from < t_size);
            assert(to < t_size);
            m_array[from][to]++;
        }

        void clear()
        {
            std::memset(m_array, 0, sizeof(m_array));
        }
    };

    template <size_t t_size>
    struct array_view_graph
    {
        typedef size_t array_t[t_size][t_size];
        array_t *m_parray;

        array_view_graph(array_t *pa) : m_parray(pa)
        {
        }

        size_t adjacent(vertex_size_t from, vertex_size_t to) const
        {
            assert(from < t_size);
            assert(to < t_size);
            return (*m_parray)[from][to];
        }

        void add_edge(vertex_size_t from, vertex_size_t to, cost_t cost = 1)
        {
            assert(from < t_size);
            assert(to < t_size);
            (*m_parray)[from][to]++;
        }

        void clear()
        {
            std::memset(*m_parray, 0, sizeof(*m_parray));
        }
    };

    template <vertex_size_t t_size>
    struct from_centric_graph
    {
        std::vector<from_centric_edge> m_from_edges[t_size];

        size_t adjacent(vertex_size_t from, vertex_size_t to) const
        {
            assert(from < t_size);
            assert(to < t_size);
            size_t num = 0;
            for (auto& edge : m_from_edges[from])
            {
                if (edge.to == to)
                    ++num;
            }
            return num;
        }

        void add_edge(vertex_size_t from, vertex_size_t to, cost_t cost = 1)
        {
            assert(from < t_size);
            assert(to < t_size);
            m_from_edges[from].push_back(from_centric_edge { to });
        }

        void clear()
        {
            for (size_t i = 0; i < t_size; ++i)
            {
                m_from_edges[i].clear();
            }
        }
    };

    template <vertex_size_t t_size>
    struct from_centric_cost_graph
    {
        std::vector<from_centric_cost_edge> m_from_edges[t_size];

        size_t adjacent(vertex_size_t from, vertex_size_t to) const
        {
            assert(from < t_size);
            assert(to < t_size);
            size_t num = 0;
            for (auto& edge : m_from_edges[from])
            {
                if (edge.to == to)
                    ++num;
            }
            return num;
        }

        void add_edge(vertex_size_t from, vertex_size_t to, cost_t cost = 1)
        {
            assert(from < t_size);
            assert(to < t_size);
            m_from_edges[from].push_back(from_centric_cost_edge { to, cost });
        }

        void clear()
        {
            for (vertex_size_t i = 0; i < t_size; ++i)
            {
                m_from_edges[i].clear();
            }
        }
    };

    //////////////////////////////////////////////////////////////////////////
    // Adjacent

    template <size_t t_size>
    inline size_t adjacent(const array_graph<t_size>& graph,
                           vertex_size_t from, vertex_size_t to)
    {
        return graph.adjacent(from, to);
    }

    template <size_t t_size>
    inline size_t&
    adjacent(const array_view_graph<t_size>& graph, vertex_size_t from, vertex_size_t to)
    {
        return graph.adjacent(from, to);
    }

    template <size_t t_size>
    inline size_t
    adjacent(const from_centric_graph<t_size>& graph,
             vertex_size_t from, vertex_size_t to)
    {
        return graph.adjacent(from, to);
    }

    template <vertex_size_t t_size>
    inline size_t
    adjacent(const from_centric_cost_graph<t_size>& graph,
             vertex_size_t from, vertex_size_t to)
    {
        return graph.adjacent(from, to);
    }

    //////////////////////////////////////////////////////////////////////////
    // Count 2

    template <size_t t_num>
    inline size_t count_of(const array_graph<t_num>&)
    {
        return t_num;
    }

    template <size_t t_num>
    inline size_t count_of(const array_view_graph<t_num>&)
    {
        return t_num;
    }

    template <size_t t_num>
    inline size_t count_of(const from_centric_graph<t_num>&)
    {
        return t_num;
    }

    template <size_t t_num>
    inline size_t count_of(const from_centric_cost_graph<t_num>&)
    {
        return t_num;
    }

    //////////////////////////////////////////////////////////////////////////
    // Copy Graph

    template <size_t t_size>
    inline void copy_graph(array_graph<t_size>& g1, const array_graph<t_size>& g2)
    {
        g1 = g2;
    }

    template <size_t t_size>
    inline void copy_graph(array_view_graph<t_size>& g1, const array_view_graph<t_size>& g2)
    {
        g1.m_parray = g2.m_parray;
    }

    template <size_t t_size>
    inline void copy_graph(array_graph<t_size>& g1, const array_view_graph<t_size>& g2)
    {
        g1.m_array = *g2.m_parray;
    }

    template <size_t t_size>
    inline void copy_graph(from_centric_graph<t_size>& g1,
                           const from_centric_graph<t_size>& g2)
    {
        g1 = g2;
    }

    template <size_t t_size>
    inline void copy_graph(array_graph<t_size>& g1,
                           const from_centric_graph<t_size>& g2)
    {
        for (vertex_size_t from = 0; from < t_size; ++from)
        {
            for (vertex_size_t to = 0; to < t_size; ++to)
            {
                g1.m_array[from][to] = adjacent(g2, from, to);
            }
        }
    }

    template <size_t t_size>
    inline void copy_graph(array_view_graph<t_size>& g1,
                           const from_centric_graph<t_size>& g2)
    {
        for (vertex_size_t from = 0; from < t_size; ++from)
        {
            for (vertex_size_t to = 0; to < t_size; ++to)
            {
                (*g1.m_parray)[from][to] = adjacent(g2, from, to);
            }
        }
    }

    template <size_t t_size>
    inline void copy_graph(from_centric_cost_graph<t_size>& g1,
                           const from_centric_cost_graph<t_size>& g2)
    {
        g1 = g2;
    }

    template <size_t t_size>
    inline void copy_graph(array_graph<t_size>& g1,
                           const from_centric_cost_graph<t_size>& g2)
    {
        for (vertex_size_t from = 0; from < t_size; ++from)
        {
            for (vertex_size_t to = 0; to < t_size; ++to)
            {
                g1.m_array[from][to] = adjacent(g2, from, to);
            }
        }
    }

    template <size_t t_size>
    inline void copy_graph(array_view_graph<t_size>& g1,
                           const from_centric_cost_graph<t_size>& g2)
    {
        for (vertex_size_t from = 0; from < t_size; ++from)
        {
            for (vertex_size_t to = 0; to < t_size; ++to)
            {
                (*g1.m_parray)[from][to] = adjacent(g2, from, to);
            }
        }
    }

    template <typename T_GRAPH, typename T_EDGE>
    inline void
    copy_graph(T_GRAPH& g1, const std::vector<T_EDGE>& E)
    {
        g1.clear();
        for (auto& edge : E)
        {
            g1.add_edge(edge.from, edge.to, cost_of(edge));
        }
    }

    //////////////////////////////////////////////////////////////////////////
    // Algorithm

    template <typename T_SEE_FLAGS, typename T_GRAPH>
    inline void dfs(T_SEE_FLAGS& seen, vertex_size_t iv, const T_GRAPH& graph)
    {
        assert(count_of(seen) == count_of(graph));
        assert(iv < count_of(seen));
        seen[iv] = 1;
        for (vertex_size_t jv = 0; jv < count_of(seen); jv++)
        {
            if (!seen[jv] && adjacent(graph, iv, jv))
                dfs(seen, jv, graph);
        }
    }

    template <typename T_SEE_FLAGS, typename T_GRAPH>
    inline void dfs2(T_SEE_FLAGS& seen, vertex_size_t iv, const T_GRAPH& graph)
    {
        assert(count_of(seen) == count_of(graph));
        assert(iv < count_of(seen));
        std::stack<edge_size_t> stack;
        seen[iv] = 1;
        stack.push(iv);
        while (!stack.empty())
        {
            vertex_size_t v = stack.top();
            stack.pop();
            for (size_t jv = 0; jv < count_of(seen); jv++)
            {
                if (!seen[jv] && adjacent(graph, v, jv))
                {
                    seen[jv] = 1;
                    stack.push(jv);
                }
            }
        }
    }

    template <typename T_GRAPH>
    bool is_connected_graph(const T_GRAPH& graph)
    {
        std::vector<edge_size_t> seen;
        seen.resize(count_of(graph));
        dfs(seen, 0, graph);
        for (auto& item : seen)
        {
            if (!item)
                return false;
        }
        return true;
    }

    template <typename T_GRAPH>
    bool is_connected_graph2(const T_GRAPH& graph)
    {
        std::vector<edge_size_t> seen;
        seen.resize(count_of(graph));
        dfs2(seen, 0, graph);
        for (auto& item : seen)
        {
            if (!item)
                return false;
        }
        return true;
    }

    template <typename T_DISTANCES, typename T_EDGE>
    inline void
    bellman_ford(vertex_size_t start, T_DISTANCES& dists,
                 const std::vector<vertex_size_t>& V,
                 const std::vector<T_EDGE>& E)
    {
        assert(count_of(dists) == count_of(V));
        std::fill(&dists[0], &dists[count_of(dists)], max_cost);
        dists[start] = 0;
        for (vertex_size_t i = 0; i < count_of(V); ++i)
        {
            for (edge_size_t j = 0; j < count_of(E); ++j)
            {
                auto& e = E[j];
                auto& old_dist = dists[e.to];
                auto new_dist = dists[e.from] + cost_of(e);
                if (old_dist > new_dist)
                {
                    old_dist = new_dist;
                    assert(i != count_of(V) - 1);
                }
            }
        }
    }

    template <typename T_DISTANCES, vertex_size_t t_size>
    void dijkstra(size_t start, T_DISTANCES& dists,
                  const from_centric_cost_graph<t_size>& graph)
    {
        assert(count_of(dists) == t_size);
        typedef std::pair<typename T_DISTANCES::value_type, vertex_size_t> PAIR;

        std::fill(&dists[0], &dists[count_of(dists)], max_cost);
        dists[start] = 0;

        std::priority_queue<PAIR, std::vector<PAIR>, std::greater<PAIR> > queue;
        queue.push(std::make_pair(0, start));

        while (!queue.empty())
        {
            PAIR pair = queue.top();
            queue.pop();

            if (dists[pair.second] < pair.first)
                continue;

            for (size_t i = 0; i < graph.m_from_edges[pair.second].size(); ++i)
            {
                const from_centric_cost_edge& e = graph.m_from_edges[pair.second][i];
                auto& old_dist = dists[e.to];
                auto new_dist = dists[pair.second] + cost_of(e);
                if (old_dist > new_dist)
                {
                    old_dist = new_dist;
                    queue.push(std::make_pair(dists[e.to], e.to));
                }
            }
        }
    }

    template <vertex_size_t t_size, typename T_GRAPH, typename T_COST_OF_EDGE>
    void warshall_floyd(cost_t (&dp)[t_size][t_size], const T_GRAPH& graph,
                        const T_COST_OF_EDGE& cost_of)
    {
        for (vertex_size_t i = 0; i < t_size; ++i)
        {
            for (vertex_size_t j = 0; j < t_size; ++j)
            {
                if (i == j)
                {
                    dp[i][j] = 0;
                }
                else if (adjacent(graph, i, j))
                {
                    dp[i][j] = cost_of(i, j);
                }
                else
                {
                    dp[i][j] = max_cost;
                }
            }
        }

        for (vertex_size_t k = 0; k < t_size; ++k)
        {
            for (vertex_size_t i = 0; i < t_size; ++i)
            {
                for (vertex_size_t j = 0; j < t_size; ++j)
                {
                    dp[i][j] = std::min(dp[i][j], dp[i][k] + dp[k][j]);
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////
    // Unit test

    inline void unittest()
    {
        static size_t array1[2][2] = { 1, 0, 0, 1};
        array_graph<2> g1(array1);
        assert(!is_connected_graph(g1));
        assert(!is_connected_graph2(g1));

        static size_t array2[3][3] = { {1, 0, 1}, {0, 1, 1}, {1, 1, 0}};
        array_graph<3> g2(array2);
        assert(is_connected_graph(g2));
        assert(is_connected_graph2(g2));

        std::vector<cost_t> dists;
        dists.resize(5);
        std::vector<vertex_size_t> V = { 0, 1, 2, 3, 4 };
        std::vector<cost_edge> E = { { 0, 1, 5 }, {1, 2, 6}, {2, 3, 2}, {0, 2, 15} };
        bellman_ford(0, dists, V, E);
        assert(dists[0] == 0);
        assert(dists[1] == 5);
        assert(dists[2] == 11);
        assert(dists[3] == 13);
        assert(dists[4] == max_cost);

        from_centric_cost_graph<5> g3;
        copy_graph(g3, E);
        assert(!is_connected_graph(g3));
        assert(!is_connected_graph2(g3));
        dijkstra(0, dists, g3);
        assert(dists[0] == 0);
        assert(dists[1] == 5);
        assert(dists[2] == 11);
        assert(dists[3] == 13);
        assert(dists[4] == max_cost);

        array_graph<5> g4;
        copy_graph(g4, g3);
        assert(!is_connected_graph(g4));
        assert(!is_connected_graph2(g4));

        std::vector<cost_edge> E2 = { { 0, 1, 5 }, {1, 2, 6}, {2, 3, 2}, {3, 4, 15} };
        copy_graph(g3, E2);
        assert(is_connected_graph(g3));
        assert(is_connected_graph2(g3));
    }
} // namespace vgraph
