clear
%fidelity_outside=2*pattern_size;
%number_files=10;
number_files=777;

a='percolation_';
b=num2str(01);
filename = [a,b];
eval(filename);

rec_density_values=zeros(number_files,1);
aff_density_values=zeros(number_files,1);
max_number_pattern=zeros(number_files,1);
max_number_patterns_divided_by_density=zeros(number_files,1);
max_number_patterns_divided_by_strong_density=zeros(number_files,1);
information_per_strong_synapse = zeros(number_files,1);
information_compressed_strong_synapse = zeros(number_files,1);
strong_synapse_willshaw_benchmark = zeros(number_files,1);
compressed_benchmark = zeros(number_files,1);
compressed_model_bits_values = zeros(number_files,1);

%paffstrong=1-(1-(pattern_size/number_vertices)*(1-(1-pattern_size/number_vertices)^number_associations))^(mat_s_prec_paff_mmax(4)/number_associations)
%paffstrong=mat_s_prec_paff_mmax(3)*paffstrong
%precstrong=1-(1-(pattern_size/number_vertices)*((pattern_size-1)/(number_vertices-1)))^(mat_s_prec_paff_mmax(4)/number_associations)
%precstrong=mat_s_prec_paff_mmax(2)*precstrong

%patternperstrongdensity = mat_s_prec_paff_mmax(4)/(number_vertices*number_vertices*(paffstrong+precstrong))
%legend_strings= cell(number_files,1);

for j=1:number_files
    a='percolation_';
    b=num2str(j-1);
   
    filename = [a,b];
    %c='_simple';
    %filename=[a,b,c];
    eval(filename);
    
    %p_aff = 1;
    %p_rec = 0.1*(j-1);
    rec_density_values(j) = p_rec;
    aff_density_values(j) = p_aff;
    max_number_pattern(j) = mean(maximal_number_patterns_under_fidelity)*number_associations;
    afferent_probability_turn_strong = 1 - (1 - pattern_size/number_vertices*(1 - (1-pattern_size/number_vertices)^(number_associations)))^(max_number_pattern(j)/number_associations);
    recurrent_probability_turn_strong = 1 - (1 - pattern_size*(pattern_size-1)/(number_vertices*(number_vertices-1)))^(max_number_pattern(j)/number_associations);
    max_number_patterns_divided_by_strong_density(j) = max_number_pattern(j)/(number_vertices*number_vertices*(p_aff*afferent_probability_turn_strong+p_rec*recurrent_probability_turn_strong));
    information_per_strong_synapse(j) = max_number_patterns_divided_by_strong_density(j)*pattern_size*log2(number_vertices/pattern_size)/number_associations;
    entropy_afferent = -afferent_probability_turn_strong*log2(afferent_probability_turn_strong) - (1-afferent_probability_turn_strong)*log2(1-afferent_probability_turn_strong);
    if (recurrent_probability_turn_strong == 0) entropy_recurrent = 0;
    elseif (recurrent_probability_turn_strong == 1) entropy_recurrent = 0;
    else entropy_recurrent = -recurrent_probability_turn_strong*log2(recurrent_probability_turn_strong) - (1-recurrent_probability_turn_strong)*log2(1-recurrent_probability_turn_strong);
    end
    compressed_model_bits = number_vertices*number_vertices*(rec_density_values(j)*entropy_recurrent + aff_density_values(j)*entropy_afferent);
    compressed_model_bits_values(j) = compressed_model_bits;
    information_compressed_strong_synapse(j) = max_number_pattern(j)*pattern_size*log2(number_vertices/pattern_size)/compressed_model_bits/number_associations;
    compressed_benchmark(j) = log(2)/number_associations;
    %compressed_real_benchmark(j) = log(2)*n^2/k^2
    strong_synapse_willshaw_benchmark(j) = 2*log(2)/number_associations;
    %define varibels
    %number_patterns=max_number_pattern(j);
    %p_rec=best_rec_density_values(j);
    %p_aff=best_aff_density_values(j);
    %max_number_patterns_divided_by_density(j)=max_number_pattern(j)/(p_aff+p_rec);
    %compute strong density
    %   number_multi_patterns=number_patterns/ number_associations;
    %   p_aff_get_strong=1- (1 - (1-((number_vertices-pattern_size)/number_vertices)^number_associations)*pattern_size/number_vertices )^number_multi_patterns;
    %   p_aff_strong= p_aff*p_aff_get_strong;
    %   p_rec_get_strong=1- (1-pattern_size*(pattern_size-1)/number_vertices/(number_vertices-1))^number_multi_patterns; 
    %   p_rec_strong= p_rec*p_rec_get_strong;
    %end compute strong density
    %max_number_patterns_divided_by_strong_density(j)=max_number_pattern(j)/(number_vertices*number_vertices*(p_aff_strong+p_rec_strong));

end

clf
x=number_associations_values;
Y=[ best_rec_density_values,best_aff_density_values ];
Z=[max_number_pattern, max_number_patterns_divided_by_strong_density*5*10^5]; % the 10 is only to make it mor visible
%, max_number_patterns_divided_by_density
[AX, H1, H2]=plotyy(x, Y, x, Z)
set(AX,{'ycolor'},{'r';'b'})
set(H1, {'color'},{'r'; 'r'})
set(H1, {'marker'},{'*'; 'o'})
set(H1, {'markers'}, {10;10})
set(H1, {'LineStyle'},{'-'; '--'})
set(H1, {'LineWidth'},{2; 2})
set(H2, {'color'},{'b'; 'b'})
set(H2, {'marker'},{'*'; 'o'})
set(H2, {'markers'}, {10;10})
set(H2, {'LineStyle'},{'-'; '--'})
set(H2, {'LineWidth'},{2; 2})
title('Densities optimized for $\mathcal{M}_{\overline{\textrm{syn}}}$', 'Interpreter', 'latex','FontSize', 50)
h_legend=legend('Optimal $p_r$' , 'Optimal $p_a$' , '$\mathcal{M}$', '$\mathcal{M}_{\overline{\textrm{syn}}}  \cdot 5 \cdot 10^5$', 'Location', 'SouthEast');
set(h_legend,'Interpreter', 'latex','FontSize',40);
%, 'Mmax / (p_aff+p_rec)'
xlabel('$s$','Interpreter','latex', 'FontSize', 40);
ylabel(AX(1), 'Density', 'FontSize', 40);
set(AX(1), 'FontSize', 40)
set(AX(1), 'YTick', 0:0.05:0.2)
ylabel(AX(2), '$\mathcal{M}$','Interpreter','latex', 'FontSize', 40);
set(AX(2), 'FontSize', 40)
ylim(AX(2),[ 0, 2000]);
set(AX(2), 'YTick', 0:500:2000)
set(AX(1),'Position', [0.13 0.14 0.730 0.76]);


%print('fig_best_rec_density_over_number_associations_for_INI_talk', '-djpeg') 
%print('fig_best_rec_density_over_number_associations', '-djpeg') 
%savefig('fig_best_rec_density_over_number_associations') 


%%%%Information per strong synapse
clf
x=number_associations_values;
Y=[ best_rec_density_values,best_aff_density_values ];
Z=[information_per_strong_synapse, strong_synapse_willshaw_benchmark]; % the 10 is only to make it mor visible
%, max_number_patterns_divided_by_density
[AX, H1, H2]=plotyy(x, Y, x, Z)
set(AX,{'ycolor'},{'r';'b'})
set(H1, {'color'},{'r'; 'r'})
set(H1, {'marker'},{'*'; 'o'})
set(H1, {'markers'}, {10;10})
set(H1, {'LineStyle'},{'-'; '--'})
set(H1, {'LineWidth'},{2; 2})
set(H2, {'color'},{'b'; 'b'})
set(H2, {'marker'},{'*'; 'o'})
set(H2, {'markers'}, {10;10})
set(H2, {'LineStyle'},{'-'; '--'})
set(H2, {'LineWidth'},{2; 2})
title('Densities optimized for $\mathcal{M}_{\overline{\textrm{syn}}}$', 'Interpreter', 'latex','FontSize', 50)
h_legend=legend('Optimal $p_r$' , 'Optimal $p_a$' , '$\mathcal{M}$', '$\mathcal{M}_{\overline{\textrm{syn}}}  \cdot 5 \cdot 10^5$', 'Location', 'NorthEast');
set(h_legend,'Interpreter', 'latex','FontSize',40);
%, 'Mmax / (p_aff+p_rec)'
xlabel('$s$','Interpreter','latex', 'FontSize', 40);
ylabel(AX(1), 'Density', 'FontSize', 40);
set(AX(1), 'FontSize', 40)
set(AX(1), 'YTick', 0:0.05:0.2)
ylabel(AX(2), '$\mathcal{M}$','Interpreter','latex', 'FontSize', 40);
set(AX(2), 'FontSize', 40)
%ylim(AX(2),[ 0, 2000]);
%set(AX(2), 'YTick', 0:500:2000)
set(AX(1),'Position', [0.13 0.14 0.730 0.76]);



%%%%%Information compressed graph
clf
x=number_associations_values;
Y=[ best_rec_density_values,best_aff_density_values ];
Z=[information_compressed_strong_synapse, compressed_benchmark]; % the 10 is only to make it mor visible
%, max_number_patterns_divided_by_density
[AX, H1, H2]=plotyy(x, Y, x, Z)
set(AX,{'ycolor'},{'r';'b'})
set(H1, {'color'},{'r'; 'r'})
set(H1, {'marker'},{'*'; 'o'})
set(H1, {'markers'}, {10;10})
set(H1, {'LineStyle'},{'-'; '--'})
set(H1, {'LineWidth'},{2; 2})
set(H2, {'color'},{'b'; 'b'})
set(H2, {'marker'},{'*'; 'o'})
set(H2, {'markers'}, {10;10})
set(H2, {'LineStyle'},{'-'; '--'})
set(H2, {'LineWidth'},{2; 2})
title('Densities optimized for $\mathcal{M}_{\overline{\textrm{syn}}}$', 'Interpreter', 'latex','FontSize', 50)
h_legend=legend('Optimal $p_r$' , 'Optimal $p_a$' , '$\mathcal{M}$', '$\mathcal{M}_{\overline{\textrm{syn}}}  \cdot 5 \cdot 10^5$', 'Location', 'NorthEast');
set(h_legend,'Interpreter', 'latex','FontSize',40);
%, 'Mmax / (p_aff+p_rec)'
xlabel('$s$','Interpreter','latex', 'FontSize', 40);
ylabel(AX(1), 'Density', 'FontSize', 40);
set(AX(1), 'FontSize', 40)
set(AX(1), 'YTick', 0:0.05:0.2)
ylabel(AX(2), '$\mathcal{M}$','Interpreter','latex', 'FontSize', 40);
set(AX(2), 'FontSize', 40)
%ylim(AX(2),[ 0, 2000]);
%set(AX(2), 'YTick', 0:500:2000)
set(AX(1),'Position', [0.13 0.14 0.730 0.76]);