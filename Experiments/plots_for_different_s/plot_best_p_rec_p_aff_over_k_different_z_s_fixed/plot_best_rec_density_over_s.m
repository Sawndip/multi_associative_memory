%each file percolation_i contains data for specific fraction_rec_aff and
%specific number_associations
% it opens such a file and searches for the number of patterns , Mmax, we can
% reveal until fidelity_outside many vertices turn active in the whole
% graph.
% We plot Mmax over d for different s

clear
%fidelity_outside=2*pattern_size;
number_files=12;

a='percolation_';
b=num2str(1);
filename=[a,b];
eval(filename);

number_z = 2; %Check how many z values were used -> not currently stored on the files
number_pattern_sizes = number_files/number_z;

number_z_values = zeros(number_z,1);
number_pattern_size_values=zeros(number_pattern_sizes,1);
best_rec_density_values=zeros(number_pattern_sizes,number_z);
best_aff_density_values=zeros(number_pattern_sizes,number_z);
max_number_pattern=zeros(number_pattern_sizes,number_z);
max_number_patterns_divided_by_density=zeros(number_pattern_sizes,number_z);
max_number_patterns_divided_by_strong_density=zeros(number_pattern_sizes,number_z);
information_number_patterns = zeros(number_pattern_sizes,number_z);
information_sparse_willshaw_values= zeros(number_pattern_sizes,  number_z);
information_complete_willshaw_values= zeros(number_pattern_sizes,  1);
%legend_strings= cell(number_files,1);


for pattern_size_iterator=1:number_pattern_sizes


    for iterator = 1:number_z
        j = pattern_size_iterator + (iterator-1)*number_pattern_sizes;
        a='percolation_';
        b=num2str(j);
        filename=[a,b];
        eval(filename);
        number_pattern_size_values(pattern_size_iterator) = pattern_size;
        number_z_values(iterator) = act_threshold;
        best_rec_density_values(pattern_size_iterator, iterator) = mat_s_prec_paff_mmax(2);
        best_aff_density_values(pattern_size_iterator, iterator) = mat_s_prec_paff_mmax(3);
        max_number_pattern(pattern_size_iterator, iterator) = mat_s_prec_paff_mmax(4);
        information_number_patterns(pattern_size_iterator, iterator) = max_number_pattern(pattern_size_iterator, iterator)*number_pattern_size_values(pattern_size_iterator)*(log2(number_vertices/pattern_size))/(number_associations*number_vertices*number_vertices*(best_aff_density_values(pattern_size_iterator,iterator)+best_rec_density_values(pattern_size_iterator,iterator)));
        
        addpath('../../../../../../matlab/new_expectation_calculation/expectation_calc_plots/information_asymptotic_willshaw')
        fid_ins=fidelity_inside;
        fid_out= (fidelity_outside-pattern_size)/ (number_vertices-pattern_size); %information_sparse_willshaw has this fidelity definition
        information_sparse_willshaw_values(pattern_size_iterator, iterator)= information_sparse_willshaw(number_vertices, pattern_size, act_threshold, number_associations,  fid_ins, fid_out);
        information_complete_willshaw_values(pattern_size_iterator, 1)= information_willshaw(number_vertices, pattern_size, number_associations,  fid_ins, fid_out);
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

Y=[];
X=number_pattern_size_values;
for iterator=1:number_z % the 2 is for the ini talk
    Y=[Y,  best_rec_density_values(:,iterator)];    
end
for iterator=1:number_z % the 2 is for the ini talk
    Y=[Y,  best_aff_density_values(:,iterator)];    
end
H1=plot(X,Y);
AX=gca;
h_legend=legend('Optimal $p_r$, $z=10$' , 'Optimal $p_r$, $z=20$' , 'Optimal $p_a$, $z=10$', 'Optimal $p_a$, $z=20$', 'Location', 'NorthEast');

set(H1, {'color'},{'r'; 'b'; 'r';'b'})
set(H1, {'marker'},{'*'; '*'; 'o'; 'o'})
set(H1, {'markers'}, {10;10; 10;10})
set(H1, {'LineStyle'},{'--'; '--'; '-'; '-'})
set(H1, {'LineWidth'},{2; 2; 2; 2})
title('Densities optimized for $\mathcal{M}$, $s=8$','Interpreter','latex', 'FontSize', 50)
set(h_legend,'Interpreter','latex','FontSize',40);
%, 'Mmax / (p_aff+p_rec)'
xlabel('Pattern size $k$','Interpreter','latex', 'FontSize', 40);
ylabel( 'Density','Interpreter','latex', 'FontSize', 40);
set(AX, 'FontSize', 40)
%print('fig_maximal_number_of_patterns_for_INI_talk', '-djpeg') 

%information plot:

clf

Y=[];
X=number_pattern_size_values;
for iterator=1:number_z % the 2 is for the ini talk
    Y=[Y,  information_number_patterns(:,iterator)];    
end
for iterator=1:number_z % the 2 is for the ini talk
    Y=[Y,  information_sparse_willshaw_values(:,iterator)];    
end
Y=[Y,  information_complete_willshaw_values(:,1)];
%for iterator=1:number_z % the 2 is for the ini talk
 %   Y=[Y,  best_aff_density_values(:,iterator)];    
%end
H1=plot(X,Y);
AX=gca;
h_legend=legend('Information capacity z=10' , 'Information capacity z=20', 'Sparse Willshaw z=10', 'Sparse Willshaw z=20', 'Complete Willshaw z=k','Location', 'NorthEast');

set(H1, {'color'},{'r'; 'b'; 'r';'b'; 'g'})
set(H1, {'marker'},{'*'; '*'; 'o'; 'o';'x'})
set(H1, {'markers'}, {10;10;10;10;10})
set(H1, {'LineStyle'},{'--'; '--'; '-'; '-'; '-'})
set(H1, {'LineWidth'},{2; 2;2;2;2})
title('Optimal Densities s=8', 'FontSize', 25)
set(h_legend,'FontSize',20);
%, 'Mmax / (p_aff+p_rec)'
xlabel('Pattern Size', 'FontSize', 20);
ylabel( 'Information per synapse', 'FontSize', 20);
set(AX, 'FontSize', 20)
%print('fig_maximal_number_of_patterns_for_INI_talk', '-djpeg') 



