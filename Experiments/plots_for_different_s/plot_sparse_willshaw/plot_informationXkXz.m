%each file percolation_i contains data for specific fraction_rec_aff and
%specific number_associations
% it opens such a file and searches for the number of patterns , Mmax, we can
% reveal until fidelity_outside many vertices turn active in the whole
% graph.
% We plot Mmax over d for different s

clear
%fidelity_outside=2*pattern_size;
number_files=24;

a='percolation_';
b=num2str(1);
filename=[a,b];
eval(filename);

number_lines = 4; %Check how many z values were used -> not currently stored on the files
number_pattern_sizes = number_files/number_lines;

number_lines_values = zeros(number_lines,1);
number_pattern_size_values=zeros(number_pattern_sizes,1);
information_per_synapse_for_plot = zeros(number_pattern_sizes,number_lines);
%legend_strings= cell(number_files,1);


    for pattern_size_iterator=1:number_pattern_sizes


        for iterator = 1:number_lines
            j = pattern_size_iterator + (iterator-1)*number_pattern_sizes;
            a='percolation_';
            b=num2str(j);
            filename=[a,b];
            eval(filename);
            number_pattern_size_values(pattern_size_iterator) = pattern_size;
           
            number_patterns= mean(maximal_number_patterns_under_fidelity);
            max_number_pattern(pattern_size_iterator, iterator) =number_patterns;
            error_prob_ins= mean(inside_error_probability);
            error_prob_out=mean(outside_error_probability);
            p10=error_prob_ins;
            p01=error_prob_out;
            addpath('../../../../../../matlab/new_expectation_calculation/expectation_calc_plots/information_asymptotic_willshaw')
            p_active= (pattern_size*(1-p10)+ (number_vertices-pattern_size)*p01)/number_vertices;
            information_per_pattern= number_vertices * entropy( p_active)-pattern_size* entropy(p10)-(number_vertices-pattern_size)*entropy(p01);
            number_synapses=number_vertices^2*(p_aff+p_rec);
            information_per_synapse_for_plot(pattern_size_iterator, iterator)= information_per_pattern* number_patterns/ number_synapses;
            
            rmpath('../../../../../../matlab/new_expectation_calculation/expectation_calc_plots/information_asymptotic_willshaw')

            %define varibels
            %number_patterns=max_number_pattern(pattern_size_iterator, iterator);
            %p_rec=best_rec_density_values(pattern_size_iterator, iterator);
            %p_aff=best_aff_density_values(pattern_size_iterator, iterator);
            %max_number_patterns_divided_{'--'; '-'}by_density(pattern_size_iterator, iterator)=max_number_pattern(j, iterator)/(p_aff+p_rec);
            %compute strong density
            %number_multi_patterns=number_patterns/ number_associations;
            %p_aff_get_strong=1- (1 - (1-((number_vertices-pattern_size)/number_vertices)^number_associations)*pattern_size/number_vertices )^number_multi_patterns;
            %p_aff_strong= p_aff*p_aff_get_strong;
            %p_rec_get_strong=1- (1-pattern_size*(pattern_size-1)/number_vertices/(number_vertices-1))^number_multi_patterns; 
            %p_rec_strong= p_rec*p_rec_get_strong;
            %end compute strong density
            %max_number_patterns_divided_by_strong_density(pattern_size_iterator, iterator)=max_number_pattern(j, iterator)/(number_vertices*number_vertices*(p_aff_strong+p_rec_strong));
        end
    end

clf

Y=information_per_synapse_for_plot;
X=number_pattern_size_values;


H1=plot(X,Y);
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



