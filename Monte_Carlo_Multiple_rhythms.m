% Original Author: Fedor Shmarov

%  Edited by: Matthew Lilley
% This script is the Third iteration of the stochastic model of melatonin
% dysregulation.
% This Script Runs two rhythm "c" amount of times defined by variable "C"
% on  line 123.
% so that the intuiton behind Monte Carlo analysis can be used to determine
% the rhythms of patients that exert ambiguity.
%
% This script has all of the same functionality of
% "Artificial_Flare_Interpolation.m" but recreating patient flares is done
% stochastically defined by the code lines 117-191.
% The Gausian function is defined at the end of this script.
% This Script uses the inferred rhythmic stimulus of TNF-A and uses the
% RMSD of two Gausiians representing different melatonin rhtyms in order to
% apply an amplitude onto the rhythmic stimuli.
% In order to Run this script you will need to:
%  - Observe the arrays on lines:
%       117 - 119
%  - Within the Dysregulated Gaussians on lines:
%       164 & 167
%  - Select the index of the array that corresponds to the rhythm you are
%  using.
% - The caption is statically typed, meaning the rhythm will have to be
% named
% IMPORTANT
% In order to run this script it is Required to have the UVB Parameter
% already known. An array of all patietns and their UV sensitivitiy values
% are listed in the text file "UV_Values".
% Change the UVB parameter on line 233.


% SBML import of the model
m1 = sbmlimport('../models/psor_v8_5.xml');


% times when PASI values are recorded
time_pasis = [];
cur_pasis_time = 0;
for i=1:13
    time_pasis = [time_pasis cur_pasis_time];
    cur_pasis_time = cur_pasis_time + 7;
end

% setting the times for when the UVB doses are administered
time_doses = [];
% UVB admission pattern (days) for every week (0 = Monday, 1 = Tuesday, ...)
cur_doses_time = [0 2 4];
% generating the doses days for 11 weeks of therapy
for i=1:11
    time_doses = [time_doses cur_doses_time];
    cur_doses_time = cur_doses_time + 7;
end

% importing the data from an XLSX file
data = readtable("../data/data_matlab.xlsx");

% extracting only the patients whose IDs are in a particular range
data = data((data.ID > 0 & data.ID <= 100),:);

data_ids = [54];

% extracting PASI trajectories        
pasis = [data.PASI_PRE_TREATMENT ...
            data.PASI_END_WEEK_1 ...
            data.PASI_END_WEEK_2 ...
            data.PASI_END_WEEK_3 ...
            data.PASI_END_WEEK_4 ...
            data.PASI_END_WEEK_5 ...
            data.PASI_END_WEEK_6 ...
            data.PASI_END_WEEK_7 ...
            data.PASI_END_WEEK_8 ...
            data.PASI_END_WEEK_9 ...
            data.PASI_END_WEEK_10 ...
            data.PASI_END_WEEK_11 ...
            data.PASI_END_WEEK_12];
% extracting UVB doses        
doses = [data.UVB_DOSE_1 data.UVB_DOSE_2 data.UVB_DOSE_3 ...
            data.UVB_DOSE_4 data.UVB_DOSE_5 data.UVB_DOSE_6 ...
            data.UVB_DOSE_7 data.UVB_DOSE_8 data.UVB_DOSE_9 ...
            data.UVB_DOSE_10 data.UVB_DOSE_11 data.UVB_DOSE_12 ...
            data.UVB_DOSE_13 data.UVB_DOSE_14 data.UVB_DOSE_15 ...
            data.UVB_DOSE_16 data.UVB_DOSE_17 data.UVB_DOSE_18 ...
            data.UVB_DOSE_19 data.UVB_DOSE_20 data.UVB_DOSE_21 ...
            data.UVB_DOSE_22 data.UVB_DOSE_23 data.UVB_DOSE_24 ...
            data.UVB_DOSE_25 data.UVB_DOSE_26 data.UVB_DOSE_27 ...
            data.UVB_DOSE_28 data.UVB_DOSE_29 data.UVB_DOSE_30 ...
            data.UVB_DOSE_31 data.UVB_DOSE_32 data.UVB_DOSE_33]; 

% extracting Follow Up Pasi Values
FU_Pasi = data.LAST_FU_PASI;

% extracting Last follow up months
FU_Pasi_Month = data.LAST_FU_MONTH;
        
% obtaining scaled PASI trajectories: every PASI value is divided by the
% baseline PASI
pasis_scaled = pasis;
for i = 1:length(pasis_scaled)
%     pasis_scaled(i,:) = totC_h + (totC_p - totC_h)*(pasis_scaled(i,:)/max(pasis_scaled(i,:)));
    pasis_scaled(i,:) = pasis_scaled(i,:)/pasis_scaled(i,1);
end
% model simulation upper time bound
%stop_time = 735;
stop_time = 1000;
% list for the UVB efficacy values
uv_eff = [];
% number of PASI points used for parameter fitting
num_pasi_points = 13;
% iterating through all the IDs in the data

% generating a figure
fig = figure('Units','normalized','OuterPosition',[10 10 30 90],'Visible','off');
% setting the font name
set(gca, 'FontName', 'Arial');
for k=1:length(data.ID)
    for b=1:1
        Normal = {1, 1.5, 19, 22};
        Moderate = {1.5, 2, 16, 24};
        Severe = {2, 4, 8, 29};

        for c=1:50
            data_row = data(data.ID == data.ID(k), :);
            p_id = data.ID(k);
            % only running it for the ids in data_ids
            if(~ismember(p_id, data_ids))
                continue;
            end
             % obtaining the index of the row for the patient
            for i=1:length(data.ID)
                if(p_id == data.ID(i))
                    index = i;
                    break;
                end
            end
            
            %%%%
            delete(m1.Events);
        
            % Parameters
            mu = 19; % Mean
            sigma = 1; % Standard deviation
            num_days = 365; % Number of days to repea
            % Generate x values
            x = linspace(0, num_days * 24, 1000);
            
            % Compute the reoccurring Gaussian function for healthy and dysregulated
            y_dysregulated = zeros(size(x));
            y_healthy = zeros(size(x));
            mse_values = [];
            rmse_values = [];
            
            % for loop to generate dysregulated values and compute SSE
            for n = 1:stop_time
                % generating random values for dysregulated production
                % min val = [normal, moderate, severe]
                min_sigma = 1; % standard Deviation maximum decreses due to increasingly disturbed rhythms
                max_sigma = 1.25;
                min_mean = 19; % mean
                max_mean = 20;
                % Rhythm 1
                if b == 1
                    y_dysregulated = y_dysregulated + gaussian(x - (n - 1) * 24, randi([Normal{3}, Normal{4}]), Normal{1} + (Normal{2} - Normal{1}) * rand());
                % Rhythm 2
                elseif b == 2
                    y_dysregulated = y_dysregulated + gaussian(x - (n - 1) * 24, randi([Moderate{3}, Moderate{4}]), Moderate{1} + (Moderate{2} - Moderate{1}) * rand());                
                end
                y_healthy = y_healthy + gaussian(x - (n - 1) * 24, randi([min_mean, max_mean]), min_sigma + (max_sigma - min_sigma) * rand());
                
                % calculate sum of squared errors
                start_index = floor((n - 1) * length(x) / stop_time) + 1;
                end_index = floor(n * length(x) / stop_time);
                errors_day = y_dysregulated(start_index:end_index) - y_healthy(start_index:end_index);
                
                % Calculate MSE and RMSE for the current day
                mse = mean(errors_day .^ 2);
                rmse = sqrt(mse);
                
        
                % making all rmse vals positive
                disturbance = rmse*1000;
                
        
                
                daily_stimValues = sum([0.08, 0.09, 0.1, 0.11, 0.12, 0.11, 0.1, 0.09, ...
                                    0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, ...
                                    0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]) * disturbance;
                addevent(m1, ['time>=', num2str(n) ], ['dc_stim=' (num2str(daily_stimValues))]); 
        
            end
            % Helper function for Gaussian calculation (AT THE BOTTOM)
            
        
            %%%%
           
            
            
            % active apoptosis pariod in days
            a_time = 0.99999;
            % adding events for the UVB doses
            for i=1:30
                if ~isnan(doses(index,i))
                    addevent(m1, ['time>' num2str(time_doses(i)+300)], ['uv_dose=' num2str(doses(index,i))]);
                    addevent(m1, ['time>=' num2str(time_doses(i)+a_time+300)], 'uv_dose=0');
                end
            end
            % initialising the best error and the best UVB efficacy values
            best_err = [1e3];
            best_uv_eff = 0;
            % transoposing the time and the PASI trajectories
            Time_full = transpose(time_pasis+300);
            PASI_full = transpose(pasis_scaled(index,:));
            % truncating the PASI trajectories up to the number of points specified
            % by <num_pasi_points>
            Time_full = Time_full(1:num_pasi_points);
            PASI_full = PASI_full(1:num_pasi_points);
            % getting the indeces of the species from the simulation data by the
            % species name
            species_to_plot = ["PASI"];
            plot_index = [];
            for i=1:length(species_to_plot)
                for j=1:length(m1.Species)
                    if species_to_plot(i) == m1.Species(j).Name
                        plot_index = [plot_index j];
                        break;
                    end
                end
            end
            
            
            % simulating the model with the best UVB efficacy parameter
            m1 = sbml_set_parameter_value(m1, "uv_eff", 0.43);
            sim_data = model_sim(m1, stop_time);
            
            % plotting the simulated PASI trajectory
            line1 = plot((sim_data.Time-300)/7, sim_data.Data(:, plot_index(1))*pasis(index,1), 'LineWidth', 8, 'Color', 'black', 'LineStyle', '-');
            if b == 1
                line1.Color = "#000000";
            elseif b == 2
                line1.Color = "#FF00FF";
            end
            line1.Color(4) = 0.3;
            hold on;
            % drawing an empty line to display the legend properly
            line(NaN,NaN,'LineWidth',8,'LineStyle','none','Marker','x','MarkerSize', 100);
            % plotting the PASI data points
            scatter((time_pasis(1:end)+300-300)/7, pasis(index,1:end), 4000, 'x', 'LineWidth', 8, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
            hold on;
            disp(num2str("patient %i", c))      
        end
    end
end
% Converting Follow up Pasi Data to be plotted
        FU_Pasi_converted = cell2mat(FU_Pasi(index,1:end));
        FU_Pasi_Month_converted = FU_Pasi_Month(index,1:end);
        disp(FU_Pasi_Month_converted);
        disp(FU_Pasi_converted);
        
        % plotting Patient Follow up value at the final time they attended
        % Follow up month * 4 to convert to weeks
        % +12 for after treatment
        scatter((FU_Pasi_Month_converted * 4)+12, str2double(FU_Pasi_converted), 'Marker', 'x', 'LineWidth', 30, 'Color', 'k');
% setting the bounds for the plot
xlim([-0.5 100 ]);
ylim([0-0.1*max(pasis(index,:)) max(pasis(index,:))+0.1*max(pasis(index,:))]);
        % labelling the axes
xlabel('Time (weeks)');
ylabel('PASI');
        % changing the font size
set(gca,'FontSize',22);
        % setting the title
title(string(['ID = ' num2str(data.ID(index)) ', UVB sensitivity = ' num2str(0.43, '%.3f'), ', Rhythm: Normal']));
        % adding the legend
legend({'Model simulation', 'Patients data', 'FU Pasi'}, FontSize=8, Location="southeast");
        % saving the generated plot as a PNG file
work_dir = '../img';
saveas(fig, [work_dir 'patient_' num2str(p_id) '_w_' num2str(num_pasi_points-1)], 'png'); 
% saving the data-table with newly estimated "uv_eff" values to a XLSX file
%writetable(data, "../data/data_matlab.xlsx");

% Helper function for Guassian probability distribution calculation
 function y = gaussian(x, mu, sigma)
        y = exp(-((x - mu).^2) / (2 * sigma^2)) / (sqrt(2 * pi) * sigma);
    end
