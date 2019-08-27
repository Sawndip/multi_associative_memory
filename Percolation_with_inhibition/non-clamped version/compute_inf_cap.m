function inf_cap = compute_inf_cap( p_rec, p_aff, maximal_number_patterns_under_fidelity, number_trials_new_graph, number_vertices, pattern_size, number_associations)
    
   average_information=0;
   for i=1:number_trials_new_graph
      number_patterns=maximal_number_patterns_under_fidelity(i,1);
      number_active=maximal_number_patterns_under_fidelity(i,2);
      number_active_inside=maximal_number_patterns_under_fidelity(i,3);
      total_information = calculate_information(number_vertices, pattern_size, number_active_inside, number_active)*number_patterns;
      %afferent_probability_turn_strong = 1 - (1 - pattern_size/number_vertices*(1 - (1-pattern_size/number_vertices)^number_associations))^(number_patterns);
      %recurrent_probability_turn_strong = 1 - (1 - pattern_size*(pattern_size-1)/(number_vertices*(number_vertices-1)))^(number_patterns);
      %number_of_bits_to_store= number_vertices^2 *p_aff* entropy(afferent_probability_turn_strong) +  number_vertices*(number_vertices-1) *p_rec* entropy(recurrent_probability_turn_strong);
      number_synapses= number_vertices^2 *p_aff +  number_vertices*(number_vertices-1) *p_rec;
      average_information= average_information + total_information/number_synapses/ number_trials_new_graph;
   end
   inf_cap = average_information;
   
end