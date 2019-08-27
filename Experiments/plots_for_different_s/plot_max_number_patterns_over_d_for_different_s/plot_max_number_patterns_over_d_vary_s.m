%each file percolation_i contains data for specific fraction_rec_aff and
%specific number_associations
% it opens such a file and searches for the number of patterns , Mmax, we can
% reveal until fidelity_outside many vertices turn active in the whole
% graph.
% We plot Mmax over d for different s

clear 
addpath('../plot_max_number_patterns_over_d_for_different_s_c++output_only_mmax');

number_files=10500;

a='percolation_';
b=num2str(1);
filename=[a,b];
eval(filename);

number_trials_graph = 30;

d_values_plot= zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_values_plot = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_error_bars = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_minus_std_values_plot = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_plus_std_values_plot = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_values_graph_trials = zeros(number_trials_graph,1);

Mmax_per_synapse_values_plot = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_per_synapse_error_bars = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_per_synapse_minus_std_values_plot = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_per_synapse_plus_std_values_plot = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_per_synapse_values_graph_trials = zeros(number_trials_graph,1);

Mmax_per_strong_synapse_values_plot = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_per_strong_synapse_error_bars = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_per_strong_synapse_minus_std_values_plot = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_per_strong_synapse_plus_std_values_plot = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_per_strong_synapse_values_graph_trials = zeros(number_trials_graph,1);

legend_strings= cell(ceil(num_num_association_values/2),1);

for j=1:number_files
    
    a='percolation_';
    b=num2str(j-1);
    filename=[a,b];
    eval(filename);
    
    if mod(j, number_trials_graph)==0
        Mmax_values_graph_trials(30) = maximal_number_patterns_under_fidelity;
        
        Mmax_per_synapse_values_graph_trials(30) = maximal_number_patterns_under_fidelity/(number_vertices*number_vertices*(p_aff+p_rec));
        
        p_aff_get_strong=1-(1-(pattern_size/number_vertices)*(1-(1-pattern_size/number_vertices)^number_associations))^(maximal_number_patterns_under_fidelity/number_associations);
        p_aff_strong = p_aff_get_strong*p_aff;
        p_rec_get_strong=1-(1-(pattern_size/number_vertices)*((pattern_size-1)/(number_vertices-1)))^(maximal_number_patterns_under_fidelity/number_associations);
        p_rec_strong = p_rec_get_strong*p_rec;
        Mmax_per_strong_synapse_values_graph_trials(30) = maximal_number_patterns_under_fidelity/(number_vertices*number_vertices*(p_aff_strong+p_rec_strong));
    else
        Mmax_values_graph_trials(mod(j, number_trials_graph)) = maximal_number_patterns_under_fidelity;
        
        Mmax_per_synapse_values_graph_trials(mod(j, number_trials_graph)) = maximal_number_patterns_under_fidelity/(number_vertices*number_vertices*(p_aff+p_rec));
        
        p_aff_get_strong=1-(1-(pattern_size/number_vertices)*(1-(1-pattern_size/number_vertices)^number_associations))^(maximal_number_patterns_under_fidelity/number_associations);
        p_aff_strong = p_aff_get_strong*p_aff;
        p_rec_get_strong=1-(1-(pattern_size/number_vertices)*((pattern_size-1)/(number_vertices-1)))^(maximal_number_patterns_under_fidelity/number_associations);
        p_rec_strong = p_rec_get_strong*p_rec;
        Mmax_per_strong_synapse_values_graph_trials(mod(j, number_trials_graph)) = maximal_number_patterns_under_fidelity/(number_vertices*number_vertices*(p_aff_strong+p_rec_strong));
    end
    
    if mod(s_iterator,2)==1 && mod(j, number_trials_graph)==0

        
        d_values_plot(d_iterator, s_iterator)= fraction_rec_aff;
        Mmax_values_plot(d_iterator,s_iterator)= mean(Mmax_values_graph_trials); % mean of each individual graph 
        Mmax_standard_deviation=std(Mmax_values_graph_trials);
        Mmax_error_bars( d_iterator, s_iterator)=Mmax_standard_deviation;
        Mmax_minus_std_values_plot( d_iterator, s_iterator)=Mmax_values_plot(d_iterator,s_iterator)-Mmax_standard_deviation;
        Mmax_plus_std_values_plot( d_iterator, s_iterator)=Mmax_values_plot(d_iterator,s_iterator)+Mmax_standard_deviation;
        
        Mmax_per_synapse_values_plot(d_iterator,s_iterator)= mean(Mmax_per_synapse_values_graph_trials);
        Mmax_standard_deviation=std(Mmax_per_synapse_values_graph_trials);
        Mmax_per_synapse_error_bars( d_iterator, s_iterator)=Mmax_standard_deviation;
        Mmax_per_synapse_minus_std_values_plot( d_iterator, s_iterator)=Mmax_per_synapse_values_plot(d_iterator,s_iterator)-Mmax_standard_deviation;
        Mmax_per_synapse_plus_std_values_plot( d_iterator, s_iterator)=Mmax_per_synapse_values_plot(d_iterator,s_iterator)+Mmax_standard_deviation;

        Mmax_per_strong_synapse_values_plot(d_iterator,s_iterator)= mean(Mmax_per_strong_synapse_values_graph_trials);
        Mmax_standard_deviation=std(Mmax_per_strong_synapse_values_graph_trials);
        Mmax_per_strong_synapse_error_bars( d_iterator, s_iterator)=Mmax_standard_deviation;
        Mmax_per_strong_synapse_minus_std_values_plot( d_iterator, s_iterator)=Mmax_per_strong_synapse_values_plot(d_iterator,s_iterator)-Mmax_standard_deviation;
        Mmax_per_strong_synapse_plus_std_values_plot( d_iterator, s_iterator)=Mmax_per_strong_synapse_values_plot(d_iterator,s_iterator)+Mmax_standard_deviation;



        if (d_iterator==1 && mod(j,number_trials_graph)==0)

            str1 = 'Multi-ass. factor s= ';
            str2=num2str(number_associations);
            new_legend_string=[str1,str2];
            legend_strings{ceil(s_iterator/2)}=new_legend_string;
            %legend_strings{(s_iterator-1)*3+1}=new_legend_string;

            %str1 = 'average Mmax minus std for Multiassociation s= ';
            %str2=num2str(number_associations);
            %new_legend_string=[str1,str2];
            %legend_strings{(s_iterator-1)*3+2}=new_legend_string;


            %str1 = 'average Mmax plus std for Multiassociation s= ';
            %str2=num2str(number_associations);
            %new_legend_string=[str1,str2];
            %legend_strings{(s_iterator-1)*3+3}=new_legend_string;  

        end
    end
end
for j=1:number_files
    a='percolation_';
    b=num2str(j-1);
    filename=[a,b];
    eval(filename);
  
    if mod(s_iterator,2)==1 && mod(j, number_trials_graph)==0

   

        if (d_iterator==1 && mod(j,number_trials_graph)==0)

            str1 = '$s = ';
            str2=num2str(number_associations);
            str3='$';
            new_legend_string=[str1,str2,str3];
            
            legend_strings{ceil(s_iterator/2)}=new_legend_string;


        end
    end
end

rmpath('../plot_max_number_patterns_over_d_for_different_s_c++output_only_mmax')


clf

Y=[];
E=[];
X=[];
Z=[];
for i=1:2:num_num_association_values % the 2 is for the ini talk
    Y=[Y,  Mmax_values_plot(:,i)];
    E=[E,  Mmax_error_bars(:,i)];
    X=[X, d_values_plot(:,1)];
    Z=[Z, [Mmax_minus_std_values_plot(:,i); transpose(fliplr(transpose(Mmax_plus_std_values_plot(:,i))))]];
    %plot(d_values_plot(:,i), Mmax_values_plot(:,i))%, d_values_plot(:,i), Mmax_minus_std_values_plot(:,i), d_values_plot(:,i), Mmax_plus_std_values_plot(:,i))
        
end
H1=plot(X,Y);
AX=gca;
h_legend=legend(legend_strings);

set(H1, {'color'},{'r'; 'k'; 'b';'g'})
set(H1, {'marker'},{'*'; 'o'; 'd'; 'x'})
set(H1, {'markers'}, {10;10; 10;10})
set(H1, {'LineStyle'},{'-'; '--'; ':'; '-.'})
set(H1, {'LineWidth'},{2; 2; 2; 2})
%title('Maximal number of patterns $\mathcal{M}$','Interpreter', 'latex', 'FontSize', 50)
set(h_legend,'Interpreter', 'latex','FontSize',40);
%, 'Mmax / (p_aff+p_rec)'
xlabel('$p_a / p_r$','Interpreter', 'latex', 'FontSize', 40);
ylabel( '$\mathcal{M}$','Interpreter', 'latex', 'FontSize', 40);
set(AX, 'FontSize', 40)
%print('fig_maximal_number_of_patterns_for_INI_talk', '-djpeg') 


% To highlight the region between 0 and 1 in the x-axis.
% Add lines
%h1 = line([0 0],[1 1800]);
%h2 = line([1 1],[1 1800]);

% Set properties of lines
%set([h1 h2],'Color','k')
%set(Mmax_minus_std_values_plot Mmax_plus_std_values_plot, 'Color', 'k')
% Add a patch
%patch([0 1 1 0],[1 1 1800 1800], [1 1 0.7])


%colors = [[1 0.5 0.5];[0.5 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
%colors = [0;0;0;0];
patch([X(:,1);transpose(fliplr(transpose(X(:,1))))], Z(:,1), [1 0.8 0.8], 'EdgeColor', 'none')
patch([X(:,2);transpose(fliplr(transpose(X(:,2))))], Z(:,2), [0.8 0.8 0.8], 'EdgeColor', 'none')
patch([X(:,3);transpose(fliplr(transpose(X(:,3))))], Z(:,3), [0.8 0.8 1], 'EdgeColor', 'none')
patch([X(:,4);transpose(fliplr(transpose(X(:,4))))], Z(:,4), [0.8 1 0.8], 'EdgeColor', 'none')

% The order of the "children" of the plot determines which one appears on top.
% I need to flip it here.
set(gca,'children',flipud(get(gca,'children')))





% Plot for Mmax_per_synapse
clf

Y=[];
E=[];
X=[];
Z=[];
for i=1:2:num_num_association_values % the 2 is for the ini talk
    Y=[Y,  Mmax_per_synapse_values_plot(:,i)];
    E=[E, Mmax_per_synapse_error_bars(:,i)];
    X=[X, d_values_plot(:,1)];
    Z=[Z, [Mmax_per_synapse_minus_std_values_plot(:,i); transpose(fliplr(transpose(Mmax_per_synapse_plus_std_values_plot(:,i))))]];
    %plot(d_values_plot(:,i), Mmax_values_plot(:,i))%, d_values_plot(:,i), Mmax_minus_std_values_plot(:,i), d_values_plot(:,i), Mmax_plus_std_values_plot(:,i))
    
    
end
H1=plot(X,Y);
AX=gca;
h_legend=legend(legend_strings);

set(H1, {'color'},{'r'; 'k'; 'b';'g'})
set(H1, {'marker'},{'*'; 'o'; 'd'; 'x'})
set(H1, {'markers'}, {10;10; 10;10})
set(H1, {'LineStyle'},{'-'; '--'; ':'; '-.'})
set(H1, {'LineWidth'},{2; 2; 2; 2})
%title('\textbf{Maximal number of associations} $\mathcal{M}$','Interpreter', 'latex', 'FontSize', 50)
set(h_legend,'Interpreter', 'latex','FontSize',40);
%, 'Mmax / (p_aff+p_rec)'
xlabel('$p_a / p_r$','Interpreter', 'latex', 'FontSize', 40);
ylabel( '$\mathcal{M}_{\textrm{syn}}$','Interpreter', 'latex', 'FontSize', 40);
set(AX, 'FontSize', 40)
%print('fig_maximal_number_of_patterns_for_INI_talk', '-djpeg') 

% To highlight the region between 0 and 1 in the x-axis.
% Add lines
%h1 = line([0 0],[1 1800]);
%h2 = line([1 1],[1 1800]);

% Set properties of lines
%set([h1 h2],'Color','k')
%set(Mmax_minus_std_values_plot Mmax_plus_std_values_plot, 'Color', 'k')
% Add a patch
%patch([0 1 1 0],[1 1 1800 1800], [1 1 0.7])


%colors = [[1 0.5 0.5];[0.5 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
%colors = [0;0;0;0];
patch([X(:,1);transpose(fliplr(transpose(X(:,1))))], Z(:,1), [1 0.8 0.8], 'EdgeColor', 'none')
patch([X(:,2);transpose(fliplr(transpose(X(:,2))))], Z(:,2), [0.8 0.8 0.8], 'EdgeColor', 'none')
patch([X(:,3);transpose(fliplr(transpose(X(:,3))))], Z(:,3), [0.8 0.8 1], 'EdgeColor', 'none')
patch([X(:,4);transpose(fliplr(transpose(X(:,4))))], Z(:,4), [0.8 1 0.8], 'EdgeColor', 'none')

% The order of the "children" of the plot determines which one appears on top.
% I need to flip it here.
set(gca,'children',flipud(get(gca,'children')))


% Plot for Mmax_per_strong_synapse
clf

Y=[];
E=[];
X=[];
Z=[];

for i=1:2:num_num_association_values % the 2 is for the ini talk
    Y=[Y,  Mmax_per_strong_synapse_values_plot(:,i)];
    E=[E, Mmax_per_strong_synapse_error_bars(:,i)];
    X=[X, d_values_plot(:,1)];
    Z=[Z, [Mmax_per_strong_synapse_minus_std_values_plot(:,i); transpose(fliplr(transpose(Mmax_per_strong_synapse_plus_std_values_plot(:,i))))]];
    %plot(d_values_plot(:,i), Mmax_values_plot(:,i))%, d_values_plot(:,i), Mmax_minus_std_values_plot(:,i), d_values_plot(:,i), Mmax_plus_std_values_plot(:,i))
    
end
H1=plot(X,Y);
AX=gca;
h_legend=legend(legend_strings);

set(H1, {'color'},{'r'; 'k'; 'b';'g'})
set(H1, {'marker'},{'*'; 'o'; 'd'; 'x'})
set(H1, {'markers'}, {10;10; 10;10})
set(H1, {'LineStyle'},{'-'; '--'; ':'; '-.'})
set(H1, {'LineWidth'},{2; 2; 2; 2})
%title('\textbf{Maximal number of patterns} $\mathcal{M}$','Interpreter', 'latex', 'FontSize', 50)
set(h_legend,'Interpreter', 'latex','FontSize',40);
%, 'Mmax / (p_aff+p_rec)'
xlabel('$p_a / p_r$','Interpreter', 'latex', 'FontSize', 40);
ylabel( '$\mathcal{M}_{\overline{\textrm{syn}}}$','Interpreter', 'latex', 'FontSize', 40);set(AX, 'FontSize', 40)
%print('fig_maximal_number_of_patterns_for_INI_talk', '-djpeg') 


% To highlight the region between 0 and 1 in the x-axis.
% Add lines
%h1 = line([0 0],[1 1800]);
%h2 = line([1 1],[1 1800]);

% Set properties of lines
%set([h1 h2],'Color','k')
%set(Mmax_minus_std_values_plot Mmax_plus_std_values_plot, 'Color', 'k')
% Add a patch
%patch([0 1 1 0],[1 1 1800 1800], [1 1 0.7])


%colors = [[1 0.5 0.5];[0.5 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
%colors = [0;0;0;0];
patch([X(:,1);transpose(fliplr(transpose(X(:,1))))], Z(:,1), [1 0.8 0.8], 'EdgeColor', 'none')
patch([X(:,2);transpose(fliplr(transpose(X(:,2))))], Z(:,2), [0.8 0.8 0.8], 'EdgeColor', 'none')
patch([X(:,3);transpose(fliplr(transpose(X(:,3))))], Z(:,3), [0.8 0.8 1], 'EdgeColor', 'none')
patch([X(:,4);transpose(fliplr(transpose(X(:,4))))], Z(:,4), [0.8 1 0.8], 'EdgeColor', 'none')

% The order of the "children" of the plot determines which one appears on top.
% I need to flip it here.
set(gca,'children',flipud(get(gca,'children')))


%savefig('fig_number_patterns_over_d_vary_s')
%print('fig_number_patterns_over_d_vary_s', '-djpeg') 

    %plot(p_rec_values(:,1), p_aff_values(:,1), p_rec_values(:,2), p_aff_values(:,2),        p_rec_values(:,3), p_aff_values(:,3),        p_rec_values(:,4), p_aff_values(:,4),        p_rec_values(:,5), p_aff_values(:,5),        p_rec_values(:,6), p_aff_values(:,6),     p_rec_values(:,7), p_aff_values(:,7),        p_rec_values(:,8), p_aff_values(:,8) ) 
    
    %legend(legend_strings(1),        legend_strings(2),        legend_strings(3),        legend_strings(4),        legend_strings(5),        legend_strings(6),        legend_strings(7),        legend_strings(8)    )
    % 'average minus std',
    %  X, average_active_plus_std,
    %'average plus std' ,

