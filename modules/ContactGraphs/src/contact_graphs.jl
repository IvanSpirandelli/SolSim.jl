function construct_graph_from_configuration(centers; threshold = 1.0)
    n = length(centers)
    edges = [(i,j) for i in 1:n for j in i:n if are_balls_intersecting(centers[i], centers[j], threshold, threshold) && i != j]
    g = Graph(n)
    for edge in edges
        add_edge!(g,edge[1], edge[2])
    end
    g
end

function get_elastic_graph_with_target_edge_number(centers, min_edge_dist, max_edge_dist, target_edges; epsilon = 0.01)
    n = length(centers)
    last_edge_dist = min_edge_dist
    last_edge_num = 0
    while last_edge_dist <= max_edge_dist
        potential_edges = [(i,j) for i in 1:n for j in i+1:n if euclidean(centers[i], centers[j]) <= last_edge_dist]
        if length(potential_edges) < target_edges
            last_edge_num = length(potential_edges) 
            last_edge_dist += epsilon
        elseif length(potential_edges) == target_edges
            return true, last_edge_dist, construct_graph_from_configuration(centers; threshold = last_edge_dist / 2.0)
        elseif length(potential_edges) > target_edges
            return false, last_edge_dist, construct_graph_from_configuration(centers; threshold = 0.0)
        end
    end
    return false, last_edge_dist, construct_graph_from_configuration(centers; threshold = 0.0)
end