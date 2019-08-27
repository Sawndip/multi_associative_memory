function information= calculate_information(number_vertices, pattern_size, number_active_pattern, number_active_whole_graph)
    p_active=number_active_whole_graph/number_vertices;
    p10=(pattern_size- number_active_pattern)/pattern_size;
    p01=(number_active_whole_graph-number_active_pattern)/(number_vertices-pattern_size);
    information= number_vertices * entropy( p_active)-pattern_size* entropy(p10)-(number_vertices-pattern_size)*entropy(p01);
end




%information= calculate_information(10000000, 16, 16, 16)