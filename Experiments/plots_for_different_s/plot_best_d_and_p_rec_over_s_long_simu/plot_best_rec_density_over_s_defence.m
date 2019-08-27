%each file percolation_i contains data for specific fraction_rec_aff and
%specific number_associations
% it opens such a file and searches for the number of patterns , Mmax, we can
% reveal until fidelity_outside many vertices turn active in the whole
% graph.
% We plot Mmax over d for different s

clear
%fidelity_outside=2*pattern_size;
number_files=10;

a='percolation_';
b=num2str(1);
filename=[a,b];
eval(filename);

number_associations_values=zeros(number_files,1);
best_rec_density_values=zeros(number_files,1);
best_aff_density_values=zeros(number_files,1);
max_number_pattern=zeros(number_files,1);
max_number_patterns_divided_by_density=zeros(number_files,1);
max_number_patterns_divided_by_strong_density=zeros(number_files,1);



%legend_strings= cell(number_files,1);


for j=1:number_files
    a='percolation_';
    b=num2str(j);
    filename=[a,b];
    eval(filename);
    
    
    number_associations_values(j)=mat_s_prec_paff_mmax(1);
    best_rec_density_values(j)=mat_s_prec_paff_mmax(2);
    best_aff_density_values(j)=mat_s_prec_paff_mmax(3);
    max_number_pattern(j)=mat_s_prec_paff_mmax(4);
    
    %define varibels
    number_patterns=max_number_pattern(j);
    p_rec=best_rec_density_values(j);
    p_aff=best_aff_density_values(j);
    max_number_patterns_divided_by_density(j)=max_number_pattern(j)/(number_vertices*number_vertices*(p_aff+p_rec));
    %compute strong density
       number_multi_patterns=number_patterns/ number_associations;
       p_aff_get_strong=1- (1 - (1-((number_vertices-pattern_size)/number_vertices)^number_associations)*pattern_size/number_vertices )^number_multi_patterns;
       p_aff_strong= p_aff*p_aff_get_strong;
       p_rec_get_strong=1- (1-pattern_size*(pattern_size-1)/number_vertices/(number_vertices-1))^number_multi_patterns; 
       p_rec_strong= p_rec*p_rec_get_strong;
    %end compute strong density
    max_number_patterns_divided_by_strong_density(j)=max_number_pattern(j)/(number_vertices*number_vertices*(p_aff_strong+p_rec_strong));

end
    

clf
x=number_associations_values;
Y=[ best_rec_density_values,best_aff_density_values ];
%, max_number_patterns_divided_by_density
H1 = plot(x, Y)
%set(AX,{'ycolor'},{'r'})
set(H1, {'color'},{[0 0.5 0]; [0 0 0.5]})
set(H1, {'marker'},{'*'; 'o'})
set(H1, {'markers'}, {10;10})
set(H1, {'LineStyle'},{'-'; '--'})
set(H1, {'LineWidth'},{3; 3})
%set(H2, {'color'},{'b'; 'b'; 'b'})
%set(H2, {'marker'},{'*'; 'o'; '+'})
%set(H2, {'markers'}, {10;10;10})
%set(H2, {'LineStyle'},{'-'; '--';'-.'})
%set(H2, {'LineWidth'},{2; 2; 2})
title('Densities optimized for $\mathcal{M}$', 'Interpreter', 'latex','FontSize', 50)
h_legend=legend('Optimal $p_r$' , 'Optimal $p_a$', 'Location', 'SouthEast');
set(h_legend,'Interpreter', 'latex','FontSize',35);
%, 'Mmax / (p_aff+p_rec)'
xlabel('Number associations','Interpreter','latex', 'FontSize', 40);
ylabel('Density', 'FontSize', 40);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 40);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 40);
xlim([1,10]);


%ylabel(AX(2), '$\mathcal{M}$','Interpreter','latex', 'FontSize', 40);
%set(AX(2), 'FontSize', 40)
%ylim(AX(2),[ 0, 2000]);
%set(AX(2), 'YTick', 0:500:2000)
%set(H1,'Position', [0.13 0.14 0.730 0.76]);



print('fig_best_rec_density_over_number_associations_for_defence', '-djpeg') 
%print('fig_best_rec_density_over_number_associations', '-djpeg') 
%savefig('fig_best_rec_density_over_number_associations') 

