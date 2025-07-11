function get_line(p1, p2)
    x1, y1 = p1
    x2, y2 = p2
    m = (y2 - y1) / (x2 - x1)
    n = y1 - m * x1

    m, n
end

function get_position_compared_to_line(nodes, edges, line)
    if length(edges) > 1
        return 0
    end

    edge = edges[1]
    m, n = line

    s, t = edge


    s_side = sign(m * lon(nodes[s]) + n - lat(nodes[s]))
    t_side = sign(m * lon(nodes[t]) + n - lat(nodes[t]))
    if s_side != t_side
        if (s_side == 0 && t_side > 0) || (t_side == 0 && s_side > 0)
            return 1
        elseif (s_side == 0 && t_side < 0) || (t_side == 0 && s_side < 0)
            return 2
        end
        return 0
    elseif s_side > 0
        return 1
    end

    return 2
end

function geographical_seed(nodes, fp_edges, line_s, line_t, popSize)
    line = get_line((lon(nodes[line_s]), lat(nodes[line_s])), (lon(nodes[line_t]), lat(nodes[line_t])))
    line_sets = [get_position_compared_to_line(nodes, x, line) for x in fp_edges]

    [deepcopy(line_sets) for _ in 1:popSize]
end