%each file percolation_i contains data for specific fraction_rec_aff and
%specific number_associations
% it opens such a file and searches for the number of patterns , Mmax, we can
% reveal until fidelity_outside many vertices turn active in the whole
% graph.
% We plot Mmax over d for different s

clear 


number_files=707;

a='percolation_';
b=num2str(1);
filename=[a,b];
eval(filename);

fidelity_outside=2*pattern_size;
p_rec_values_plot= zeros(num_p_rec_values, num_num_association_values);
Mmax_values_plot = zeros(num_p_rec_values, num_num_association_values);
Mmax_minus_stp_rec_values_plot =  zeros(num_p_rec_values, num_num_association_values);
Mmax_plus_stp_rec_values_plot =  zeros(num_p_rec_values, num_num_association_values);

legend_strings= cell(num_num_association_values,1);


for j=1:number_files
    a='percolation_';
    b=num2str(j);
    filename=[a,b];
    eval(filename);
    
    
    p_rec_values_plot(p_rec_iterator, s_iterator)= p_rec;
    Mmax_values_plot(p_rec_iterator,s_iterator)= mean(maximal_number_patterns_under_fidelity); % mean of each individual graph 
    Mmax_standard_deviation=std(maximal_number_patterns_under_fidelity);
    Mmax_minus_stp_rec_values_plot( p_rec_iterator, s_iterator)=Mmax_values_plot(p_rec_iterator,s_iterator)-Mmax_standard_deviation;
    Mmax_plus_stp_rec_values_plot( p_rec_iterator, s_iterator)=Mmax_values_plot(p_rec_iterator,s_iterator)+Mmax_standard_deviation;

  
    

    if (p_rec_iterator==1)
        str1 = 'average Mmax for Multiassociation s= ';
        str2=num2str(number_associations);
        new_legend_string=[str1,str2];
        legend_strings{s_iterator}=new_legend_string;
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
    

clf
for i=1:num_num_association_values
    plot(p_rec_values_plot(:,i), Mmax_values_plot(:,i))%, p_rec_values_plot(:,i), Mmax_minus_stp_rec_values_plot(:,i), p_rec_values_plot(:,i), Mmax_plus_stp_rec_values_plot(:,i))
    
    hold on
end
legend(legend_strings)

xlabel('recurrent density');
ylabel('Number of Patterns');

print('fig_number_patterns_over_p_rec_vary_s', '-djpeg') 
    %plot(p_rec_values(:,1), p_aff_values(:,1), p_rec_values(:,2), p_aff_values(:,2),        p_rec_values(:,3), p_aff_values(:,3),        p_rec_values(:,4), p_aff_values(:,4),        p_rec_values(:,5), p_aff_values(:,5),        p_rec_values(:,6), p_aff_values(:,6),     p_rec_values(:,7), p_aff_values(:,7),        p_rec_values(:,8), p_aff_values(:,8) ) 
    
    %legend(legend_strings(1),        legend_strings(2),        legend_strings(3),        legend_strings(4),        legend_strings(5),        legend_strings(6),        legend_strings(7),        legend_strings(8)    )
    % 'average minus std',
    %  X, average_active_plus_std,
    %'average plus std' ,

