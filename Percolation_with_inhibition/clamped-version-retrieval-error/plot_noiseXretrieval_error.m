clear

a='percolation_';
b=num2str(1);
filename = [a,b];
eval(filename);


false_activation_noise_values = zeros(number_false_activation_noise_values,1);

p_rec_values = zeros(number_false_activation_noise_values,1);
p_aff_values = zeros(number_false_activation_noise_values,1);
pattern_stored_values = zeros(number_patterns_values,1);
outside_retrieval_error_with_rec_edges = zeros(number_false_activation_noise_values,number_patterns_values);
outside_retrieval_error_sparse_willshaw = zeros(number_false_activation_noise_values,number_patterns_values);
%opt_inf_cap_willshaw_clamped_version = zeros(number_false_activation_noise_values,1);
%average_act_threshold_clamped_version = zeros(number_false_activation_noise_values,1);
%average_act_threshold_clamped_version = zeros(number_false_activation_noise_values,1);
%optimal_files_clamped_version = zeros(number_false_activation_noise_values,1);

number_act_threshold_values = 1; % should have been on the c++ file!!!!
number_p_rec_values = 1;
for i = 1 : number_false_activation_noise_values
    for j = 1 : number_patterns_values
        file_number = j-1 + (i-1)*number_patterns_values;
        a='percolation_';
        b=num2str(file_number);
        filename = [a,b];
        eval(filename);
        false_activation_noise_values(i,1) = false_activation_noise;
        pattern_stored_values(j,1) = maximal_number_patterns_under_fidelity_average;
        outside_retrieval_error_sparse_willshaw(i,j) = mean(maximal_number_patterns_under_fidelity(:,2)-maximal_number_patterns_under_fidelity(:,3))/number_vertices;
        file_number = j-1 + (i-1 + 6)*number_patterns_values;
        a='percolation_';
        b=num2str(file_number);
        filename = [a,b];
        eval(filename);
        false_activation_noise_values(i,1) = false_activation_noise;
        outside_retrieval_error_with_rec_edges(i,j) = mean(maximal_number_patterns_under_fidelity(:,2)-maximal_number_patterns_under_fidelity(:,3))/number_vertices;
    end
end

clf

Y=[outside_retrieval_error_with_rec_edges(:,2), outside_retrieval_error_sparse_willshaw(:,2),outside_retrieval_error_with_rec_edges(:,4), outside_retrieval_error_sparse_willshaw(:,4)];
X=[false_activation_noise_values, false_activation_noise_values, false_activation_noise_values, false_activation_noise_values];


H1=plot(X,Y);
AX=gca;
h_legend=legend('$\mathcal{M}_{\textrm{syn}}=6\times 10^{-5}, p_r = 0.15$','$\mathcal{M}_{\textrm{syn}}=6\times 10^{-5}, p_r=0$' , '$\mathcal{M}_{\textrm{syn}}=12\times 10^{-5}, p_r=0.15$','$\mathcal{M}_{\textrm{syn}}=12\times 10^{-5}, p_r=0$', 'Location', 'NorthWest');

set(H1, {'color'},{'r'; 'b'; 'r'; 'b'})
set(H1, {'marker'},{'*';'*'; '+'; '+'})
set(H1, {'markers'}, {10;10;10;10})
set(H1, {'LineStyle'},{'-';'-';'--';'--'})
set(H1, {'LineWidth'},{2;2;2;2})
title('Retrieval error $p_{10}$','Interpreter','latex', 'FontSize', 50)
set(h_legend,'Interpreter','latex','FontSize', 40);
%, 'Mmax / (p_aff+p_rec)'
xlabel('activation noise $\kappa$','Interpreter','latex', 'FontSize', 40);
set(AX, 'xlim', [0,1])
ylabel( '$p_{10}$','Interpreter','latex', 'FontSize', 40);
set(AX, 'FontSize', 40)
set(AX, 'ylim', [0,0.1])

%set(AX, 'YTick',[0, 4, 8, 60, 71])

%%%%%%%%%
clf

y = pattern_stored_values;
x = false_activation_noise_values;
[X,Y] = meshgrid(x,y);

zlevels=[10, 20, 40, 100, 200, 300, 400, 500, 700, 900];
figure
H1=contour(X,Y,transpose(outside_retrieval_error_sparse_willshaw),zlevels, 'LineWidth', 2, 'LabelSpacing',300);
AX=gca;
%h_legend=legend(' Number patterns $\mathcal{M}$' , 'Location', 'NorthEast');
%set(h_legend,'Interpreter','latex','FontSize', 40);
title(' Outside retrieval error sparse Willshaw ','Interpreter','latex', 'FontSize', 50)
xlabel('False activation noise','Interpreter','latex', 'FontSize', 40);
%set(AX, 'xlim', [0,1])
ylabel( 'Number of patterns stored','Interpreter','latex', 'FontSize', 40);
set(AX, 'FontSize', 40)
clabel(H1,'Fontsize', 30);
%set(AX, 'ylim', [0.0,3])


clf

y = pattern_stored_values;
x = false_activation_noise_values;
[X,Y] = meshgrid(x,y);

zlevels=[1, 2, 5, 10, 50, 200, 500, 1000, 2000, 3000];
figure
H1=contour(X,Y,transpose(outside_retrieval_error_with_rec_edges),zlevels, 'LineWidth', 2, 'LabelSpacing',300);
AX=gca;
%h_legend=legend(' Number patterns $\mathcal{M}$' , 'Location', 'NorthEast');
%set(h_legend,'Interpreter','latex','FontSize', 40);
title(' Outside retrieval error with rec edges ','Interpreter','latex', 'FontSize', 50)
xlabel('False activation noise','Interpreter','latex', 'FontSize', 40);
%set(AX, 'xlim', [0,1])
ylabel( 'Number of patterns stored','Interpreter','latex', 'FontSize', 40);
set(AX, 'FontSize', 40)
clabel(H1,'Fontsize', 30);
%set(AX, 'ylim', [0.0,3])
