%each file percolation_i contains data for specific fraction_rec_aff and
%specific number_associations
% it opens such a file and searches for the number of patterns , Mmax, we can
% reveal until fidelity_outside many vertices turn active in the whole
% graph.
% We plot Mmax over d for different s

clear 


number_files=56;

a='percolation_';
b=num2str(1);
filename=[a,b];
eval(filename);

fidelity_outside=2*pattern_size;
d_values_plot= zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_values_plot = zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_minus_std_values_plot =  zeros(num_fraction_rec_aff_values, num_num_association_values);
Mmax_plus_std_values_plot =  zeros(num_fraction_rec_aff_values, num_num_association_values);

legend_strings= cell(ceil(num_num_association_values/2),1);


for j=1:number_files
    
    a='percolation_';
    b=num2str(j);
    filename=[a,b];
    eval(filename);
    if mod(s_iterator,2)==1


        d_values_plot(d_iterator, s_iterator)= fraction_rec_aff;
        Mmax_values_plot(d_iterator,s_iterator)= mean(maximal_number_patterns_under_fidelity); % mean of each individual graph 
        Mmax_standard_deviation=std(maximal_number_patterns_under_fidelity);
        Mmax_minus_std_values_plot( d_iterator, s_iterator)=Mmax_values_plot(d_iterator,s_iterator)-Mmax_standard_deviation;
        Mmax_plus_std_values_plot( d_iterator, s_iterator)=Mmax_values_plot(d_iterator,s_iterator)+Mmax_standard_deviation;




        if (d_iterator==1)

            str1 = 'Multi-Association s= ';
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
    

clf
Y=[];
X=d_values_plot(:,1);
for i=1:2:num_num_association_values % the 2 is for the ini talk
    Y=[Y,  Mmax_values_plot(:,i)];
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
title('Maximal Number of Patterns', 'FontSize', 25)
set(h_legend,'FontSize',20);
%, 'Mmax / (p_aff+p_rec)'
xlabel('Afferent Density divided by Recurrent Density', 'FontSize', 20);
ylabel( 'Number of Patterns', 'FontSize', 20);
set(AX, 'FontSize', 20)

print('fig_maximal_number_of_patterns_for_INI_talk', '-djpeg') 






%savefig('fig_number_patterns_over_d_vary_s')
%print('fig_number_patterns_over_d_vary_s', '-djpeg') 

    %plot(p_rec_values(:,1), p_aff_values(:,1), p_rec_values(:,2), p_aff_values(:,2),        p_rec_values(:,3), p_aff_values(:,3),        p_rec_values(:,4), p_aff_values(:,4),        p_rec_values(:,5), p_aff_values(:,5),        p_rec_values(:,6), p_aff_values(:,6),     p_rec_values(:,7), p_aff_values(:,7),        p_rec_values(:,8), p_aff_values(:,8) ) 
    
    %legend(legend_strings(1),        legend_strings(2),        legend_strings(3),        legend_strings(4),        legend_strings(5),        legend_strings(6),        legend_strings(7),        legend_strings(8)    )
    % 'average minus std',
    %  X, average_active_plus_std,
    %'average plus std' ,

