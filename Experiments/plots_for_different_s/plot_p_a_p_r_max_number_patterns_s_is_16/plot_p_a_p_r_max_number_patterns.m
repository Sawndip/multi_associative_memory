%each file percolation_i contains data for specific fraction_rec_aff and
%specific number_associations
%we plot the p_a, p_r values and the M obtained in each case

clear 

number_files=50;

a='percolation_';
b=num2str(0);
filename=[a,b];
eval(filename);

num_num_association_values=1;

p_r_values_plot = zeros(num_p_r_values, num_num_association_values);
p_a_values_plot = zeros(num_p_r_values, num_num_association_values);

Mmax_values_plot = zeros(num_p_r_values, num_num_association_values);
Mmax_values_plus_std_plot = zeros(num_p_r_values, num_num_association_values);
Mmax_values_minus_std_plot = zeros(num_p_r_values, num_num_association_values);


Mmax_per_synapse_values_plot = zeros(num_p_r_values, num_num_association_values);
Mmax_per_synapse_values_plus_std_plot = zeros(num_p_r_values, num_num_association_values);
Mmax_per_synapse_values_minus_std_plot = zeros(num_p_r_values, num_num_association_values);

Mmax_per_strong_synapse_values_plot = zeros(num_p_r_values, num_num_association_values);
Mmax_per_strong_synapse_values_plus_std_plot = zeros(num_p_r_values, num_num_association_values);
Mmax_per_strong_synapse_values_minus_std_plot = zeros(num_p_r_values, num_num_association_values);

legend_strings= cell(ceil(num_num_association_values/2),1);

for j=1:number_files
    
    a='percolation_';
    b=num2str(j-1);
    filename=[a,b];
    eval(filename);
    
    p_r_values_plot(j,1) = p_rec;
    p_a_values_plot(j,1) = p_aff;
    
    Mmax_values_plot(j,1) = mean(maximal_number_patterns_under_fidelity);
    Mmax_values_plus_std_plot(j,1) = mean(maximal_number_patterns_under_fidelity) + std(maximal_number_patterns_under_fidelity);
    Mmax_values_minus_std_plot(j,1) = mean(maximal_number_patterns_under_fidelity) - std(maximal_number_patterns_under_fidelity);

    Mmax_per_synapse_values_plot(j,1) = Mmax_values_plot(j,1)/(number_vertices*number_vertices*(p_aff+p_rec));
    Mmax_per_synapse_values_plus_std_plot(j,1) = Mmax_values_plus_std_plot(j,1)/(number_vertices*number_vertices*(p_aff+p_rec));
    Mmax_per_synapse_values_minus_std_plot(j,1) = Mmax_values_minus_std_plot(j,1)/(number_vertices*number_vertices*(p_aff+p_rec));

end


clf


Y=[Mmax_values_plot, Mmax_per_synapse_values_plot*5*10^6];
X=p_r_values_plot;
Z=p_a_values_plot;
E=[[Mmax_values_minus_std_plot; transpose(fliplr(transpose(Mmax_values_plus_std_plot)))], [Mmax_per_synapse_values_minus_std_plot*5*10^6; transpose(fliplr(transpose(Mmax_per_synapse_values_plus_std_plot*5*10^6)))]];


[AX, H1, H2]=plotyy(X,Z,X,Y);
set(AX,{'ycolor'},{'r';'b'})
set(H1, {'color'},{'r'})
set(H1, {'marker'},{'*'})
set(H1, {'markers'}, {10})
set(H1, {'LineStyle'},{'-'})
set(H1, {'LineWidth'},{2})
set(H2, {'color'},{'b'; 'b'})
set(H2, {'marker'},{'o';'*'})
set(H2, {'markers'}, {10;10})
set(H2, {'LineStyle'},{'--';':'})
set(H2, {'LineWidth'},{2;2})
%title('Afferent/Recurrent Density tradeoff', 'Interpreter', 'latex', 'FontSize', 50)
h_legend=legend( '$p_a$', '$\mathcal{M}$', '$\mathcal{M}_{\textrm{syn}}  \cdot 5 \cdot 10^5$','Location', 'NorthEast');
set(h_legend,'Interpreter','latex','FontSize',40);
%, 'Mmax / (p_aff+p_rec)'
xlabel('$p_r$','Interpreter','latex', 'FontSize', 40);
ylabel(AX(1), '$p_a$','Interpreter','latex', 'FontSize', 40);
set(AX(1), 'FontSize', 40)
set(AX(1), 'YTick', 0:0.05:0.2)
set(AX(1), 'YLim', [0 0.2])
ylabel(AX(2), '$\mathcal{M}$','Interpreter','latex', 'FontSize', 40);
set(AX(2), 'FontSize', 40)
set(AX(2), 'YTick', 0:500:2000)
set(AX(2), 'YLim', [0 2000])

set(AX(1),'Position', [0.11 0.14 0.7775 0.75]);
%set(AX(2),'Position', [0.13 0.13 0.775 0.815]);

%print('fig_best_rec_density_over_number_associations_for_INI_talk', '-djpeg') 
%print('fig_best_rec_density_over_number_associations', '-djpeg') 
%savefig('fig_best_rec_density_over_number_associations') 




set(gcf, 'CurrentAxes', AX(2))
patch([X;transpose(fliplr(transpose(X)))], E(:,1), [0.8 0.8 1], 'EdgeColor', 'none')
patch([X;transpose(fliplr(transpose(X)))], E(:,2), [0.8 0.8 1], 'EdgeColor', 'none')

set(gca,'children',flipud(get(gca,'children')))

set(gcf, 'CurrentAxes', AX(1))
set(gca,'children',flipud(get(gca,'children')))
