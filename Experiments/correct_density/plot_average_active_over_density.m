% this file plots number of active vertices i inside pattern over the
% recurrent density: for exp calt, simulation average, simulation average+-
% std deviation and max min of simulation
% This file imports a file correct_density_i with raw simulation data
% This file imports a file containing expectation calculation

correct_density_1;

number_p_rec_values= size(mat_dens_rec_aff_numbactive,1)/number_samples;
p_rec_values=zeros(number_p_rec_values,1);
p_aff_values=zeros(number_p_rec_values,1);
exp_calc_num_active=zeros(number_p_rec_values,1);
average_active=zeros(number_p_rec_values,1);
average_active_plus_std=zeros(number_p_rec_values,1);
average_active_minus_std=zeros(number_p_rec_values,1);
min_active=zeros(number_p_rec_values,1);
max_active=zeros(number_p_rec_values,1);
quantile02=zeros(number_p_rec_values,1);
quantile08=zeros(number_p_rec_values,1);

current_data=zeros(number_samples, 1);
for i=1:number_p_rec_values

    cd_start= (i-1)*number_samples+1;
    cd_end= i*number_samples;
    current_data=mat_dens_rec_aff_numbactive(cd_start:cd_end, 3);
    p_rec_values(i)=mat_dens_rec_aff_numbactive(cd_start,1);
    p_aff_values(i)=mat_dens_rec_aff_numbactive(cd_start,2);
    average= mean(current_data);
    average_active(i)=average;
    standard_deviation=std(current_data);
    average_active_plus_std(i)=average + standard_deviation;
    average_active_minus_std(i)=average - standard_deviation;
    min_active(i)=min(current_data);
    max_active(i)=max(current_data);    
    quantile02(i)=quantile(current_data,0.2);
    quantile08(i)=quantile(current_data,0.8);
    
    %Calculate the expectation exp_calc
    addpath('/home/flo/Documents/SVN/csa/florian_M/large_pattern_learning/matlab/new_expectation_calculation')
    [active_in, active_out]= compute_number_active( number_vertices, pattern_size, p_aff_values(i), p_rec_values(i), act_threshold, 0, 1);
    rmpath('/home/flo/Documents/SVN/csa/florian_M/large_pattern_learning/matlab/new_expectation_calculation')
    exp_calc_num_active(i)=active_in;
end

X=p_rec_values;
plot(X, min_active, X, quantile02, X, average_active, X, quantile08, X, max_active,X, exp_calc_num_active)
legend('min active', '0.2 quantile', 'average', '0.8 quantile',  'max active', 'calculated expectation')
% X, average_active_minus_std,
% 'average minus std',
%  X, average_active_plus_std,
%'average plus std' ,

