%each file percolation_i contains data for specific fraction_rec_aff and
%specific number_associations
% it opens such a file and searches for the number of patterns , Mmax, we can
% reveal until fidelity_outside many vertices turn active in the whole
% graph.
% We plot Mmax over d for different s

clear
%fidelity_outside=2*pattern_size;
number_files=40;

a='percolation_';
b=num2str(1);
filename=[a,b];
eval(filename);

number_lines = 4; %Check how many z values were used -> not currently stored on the files
number_z_values = number_files/number_lines;

number_lines_values = zeros(number_lines,1);
z_values=zeros(number_z_values,1);
information_per_synapse_for_plot = zeros(number_z_values,number_lines);
%legend_strings= cell(number_files,1);
information_per_synapse_mean_for_plot = zeros(number_z_values,number_lines);
information_per_synapse_mean_lessstd_for_plot = zeros(number_z_values,number_lines);
information_per_synapse_mean__plusstd_for_plot = zeros(number_z_values,number_lines);

    for z_iterator=1:number_z_values
display(z_iterator);

        for iterator = 1:number_lines
            j = z_iterator + (iterator-1)*number_z_values;
            a='percolation_';
            b=num2str(j);
            filename=[a,b];
            eval(filename);
            z_values(z_iterator) = act_threshold;
            average_total_information = 0;
            number_synapses = number_vertices^2*(p_aff+p_rec);
            for graph_trial_iterator = 1:number_trials_new_graph
                total_information_graph_trial = 0;
                number_patterns= maximal_number_patterns_under_fidelity(graph_trial_iterator);
                for number_pattern_iterator = 1:number_patterns
                    p10 = (pattern_size -inside_activity(number_pattern_iterator, graph_trial_iterator))/pattern_size;
                    p01 = outside_activity(number_pattern_iterator, graph_trial_iterator)/(number_vertices - pattern_size);
                    %max_number_pattern(z_iterator, iterator) = number_patterns;
                    p_active= (inside_activity(number_pattern_iterator, graph_trial_iterator) + outside_activity(number_pattern_iterator, graph_trial_iterator))/number_vertices;
                    entropy_p_active = -p_active*log(p_active) - (1-p_active)*log(1-p_active);
                    if p_active == 0
                        entropy_p_active = 0;
                    elseif p_active == 1
                        entrop20y_p_active = 0;
                    end
                    entropy_p_10 = -p10*log(p10) - (1-p10)*log(1 - p10);
                    if p10 == 0
                        entropy_p_10 = 0;
                    elseif p10 == 1
                        entropy_p_10 = 0;
                    end
                    entropy_p_01 = -p01*log(p01) - (1-p01)*log(1 - p01);
                    if p01 == 0
                        entropy_p_01 = 0;
                    elseif p01 == 1
                        entropy_p_01 = 0;
                    end
                    total_information_graph_trial = total_information_graph_trial + number_vertices * entropy_p_active -pattern_size* entropy_p_10-(number_vertices-pattern_size)*entropy_p_01;
                end
                average_total_information = average_total_information + total_information_graph_trial/number_synapses/number_trials_new_graph;
            end
            information_per_synapse_for_plot(z_iterator, iterator)= average_total_information/8;
            
            %information_per_pattern_vector = zeros(number_trials_new_graph,1);
            %for graph_iterator = 1:number_trials_new_graph
            %    number_patterns= maximal_number_patterns_under_fidelity(graph_iterator);
            %    max_number_pattern(z_iterator, iterator) =number_patterns;
            %    error_prob_ins= inside_error_probability(graph_iterator);
            %    error_prob_out=outside_error_probability(graph_iterator);
            %    p10=error_prob_ins;
            %    p01=error_prob_out;
            %    p_active= (pattern_size*(1-p10)+ (number_vertices-pattern_size)*p01)/number_vertices;
            %    entropy_p_active = -p_active*log(p_active) - (1-p_active)*log(1-p_active);
            %    entropy_p_10 = -p10*log(p10) - (1-p10)*log(1 - p10);
            %    entropy_p_01 = -p01*log(p01) - (1-p01)*log(1 - p01);
            %    information_per_pattern_vector(graph_iterator)= number_vertices * entropy_p_active -pattern_size* entropy_p_10-(number_vertices-pattern_size)*entropy_p_01;
            %    number_synapses=number_vertices^2*(p_aff+p_rec);
            %end
            %information_per_synapse_mean_for_plot(z_iterator, iterator)= mean(information_per_pattern_vector)* number_patterns/ number_synapses;
            %information_per_synapse_mean_lessstd_for_plot(z_iterator, iterator)= (mean(information_per_pattern_vector)-std(information_per_pattern_vector))* number_patterns/ number_synapses;
            %information_per_synapse_mean__plusstd_for_plot(z_iterator, iterator)= (mean(information_per_pattern_vector)+std(information_per_pattern_vector))* number_patterns/ number_synapses;

            %define varibels
            %number_patterns=max_number_pattern(z_iterator, iterator);
            %p_rec=best_rec_density_values(z_iterator, iterator);
            %p_aff=best_aff_density_values(z_iterator, iterator);
            %max_number_patterns_divided_{'--'; '-'}by_density(z_iterator, iterator)=max_number_pattern(j, iterator)/(p_aff+p_rec);
            %compute strong density
            %number_multi_patterns=number_patterns/ number_associations;
            %p_aff_get_strong=1- (1 - (1-((number_vertices-pattern_size)/number_vertices)^number_associations)*pattern_size/number_vertices )^number_multi_patterns;
            %p_aff_strong= p_aff*p_aff_get_strong;
            %p_rec_get_strong=1- (1-pattern_size*(pattern_size-1)/number_vertices/(number_vertices-1))^number_multi_patterns; 
            %p_rec_strong= p_rec*p_rec_get_strong;
            %end compute strong density
            %max_number_patterns_divided_by_strong_density(z_iterator, iterator)=max_number_pattern(j, iterator)/(number_vertices*number_vertices*(p_aff_strong+p_rec_strong));
        end
    end

clf

Y=information_per_synapse_for_plot;
X=z_values;


H1=plot(X,Y);
AX=gca;
h_legend=legend('Sparse Willshaw, $k=71$', 'Sparse Willshaw, $k=142$','With recurrent edges, $k=71$' , 'With recurrent edges, $k=142$' ,  'Location', 'SouthEast');

set(H1, {'color'},{'r'; 'b'; 'r';'b'})
set(H1, {'marker'},{'*'; '*'; 'o'; 'o'})
set(H1, {'markers'}, {10;10; 10;10})
set(H1, {'LineStyle'},{'--'; '--'; '-'; '-'})
set(H1, {'LineWidth'},{2; 2; 2; 2})
title('Information capacity, $s=8$','Interpreter','latex', 'FontSize', 50)
set(h_legend,'Interpreter','latex','FontSize', 40);
%, 'Mmax / (p_aff+p_rec)'
xlabel('Activation threshold $z$','Interpreter','latex', 'FontSize', 40);
ylabel( 'Information capacity','Interpreter','latex', 'FontSize', 40);
set(AX, 'FontSize', 40)
set(AX, 'ylim', [0,12*10^(-3)])

%set(AX, 'YTick',[0, 4, 8, 60, 71])

%print('fig_maximal_number_of_patterns_for_INI_talk', '-djpeg') 


Y1=information_per_synapse_mean_for_plot;
Y2=information_per_synapse_mean_lessstd_for_plot;
Y3=information_per_synapse_mean__plusstd_for_plot;
X=z_values;


H1=plot(X,Y1,X,Y2, X, Y3);
AX=gca;
h_legend=legend('Sparse Willshaw z=10', 'Sparse Willshaw z=20','With recurrent edges z=10' , 'With recurrent edges z=20' ,  'Location', 'NorthEast');

set(H1, {'color'},{'r'; 'b'; 'r';'b'})
set(H1, {'marker'},{'*'; '*'; 'o'; 'o'})
set(H1, {'markers'}, {10;10; 10;10})
set(H1, {'LineStyle'},{'--'; '--'; '-'; '-'})
set(H1, {'LineWidth'},{2; 2; 2; 2})
title('Optimal Densities s=8', 'FontSize', 25)
set(h_legend,'FontSize',20);
%, 'Mmax / (p_aff+p_rec)'
xlabel('Pattern Size', 'FontSize', 20);
ylabel( 'Density', 'FontSize', 20);
set(AX, 'FontSize', 20)
%print('fig_maximal_number_of_patterns_for_INI_talk', '-djpeg') 



