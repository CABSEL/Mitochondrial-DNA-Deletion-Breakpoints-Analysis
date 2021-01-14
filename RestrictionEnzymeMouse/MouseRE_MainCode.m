%% Code description and usage instructions
   % Description: Code to perform mixture model analysis 
   % Reference: Lakshmanan et al.,2021, "Thermodynamic analysis of DNA hybridization signatures near mitochondrial DNA deletion breakpoints" (submitted to iScience). 
   % Inputs required for this code include (1) mtDNA deletion breakpoints dataset and (2) mtDNA hybridization partition function data.
   % Code outputs model prediction results as a histogram figure, similar to those in figures 2 and 3 in the manuscript. 
   % Output figures are saved in both matlab figure format (.fig) and in high resolution adobe illustrator format (.eps)  
   % numerical results are saved in excel2003 format (.xls)
   
   % Requirements:  (1) input data, (2) optimization toolbox (MATLAB), (3) ancillary functions used within the code should be present in the same folder

   % Usage instructions: 
   % (1) format input breakpoint data ( 3 columns: [serialnumber, 5' BP, 3' BP], each row is a breakpoint pair )
   %     (redundant breakpoints (two or more deletions with identical 5' and 3' breakpoints), if present, will be included only once in the dataset)
   % (2) specify the input breakpoints datset filename " Breakpoints_FileName "
   % (3) specify the correct directory to load mtDNA hybridization partition function file
   % (4) can specify output filename
   
   % Recommendation: 
   % use dataset with atleast 25 non-redundant breakpoints, within the regions for which partition function values are available
   % model predictions can have large standard deviation and can give errors for really small datasets
   
   % Output file details
   % (1) Component predictions file: 
   %         Row 1: predicted values for each component [0,5,10,15,20,25,50,75,100]
   %         Row 2: SD values for each component (calculated from samples)
   %
   % (2) Grouped predictions file:
   %         Row 1: predicted values for each group, short, medium and long
   %         Row 2: SD values for each group
   %
   % (3) Full results file:
   %         Row 1: for original dataset, predicted values for each component [0,5,10,15,20,25,50,75,100] + LnLikelihood value
   %         Rows 2-101: sampled datasets, predicted values for each component + LnLikeliood value
   %         one row for each sample set, total 100 rows for 100 samples 

   % Last modification date: 14/01/2021

clear; clc; close('all');

tic

format long;

%% Regions in mouse mtDNA for which hybridization partition function values were calculated
% do not modify
% --------------

Ori_LBR_start_position = 4500;          Ori_LBR_stop_position = 10500;

Ori_RBR_start_position = 11000;         Ori_RBR_stop_position = 16200;

% -------------------------------------------------------------------


%% user specified parameters

main_seed_value = 13 * 10000000;       % seed for random number generator, to ensure reproducible results

window_length_val = 100;               % length of sequence analysed around each breakpoint (i.e. left-breakpoint (+/-) 100  and right-breakpoint (+/-) 100 )
 
bin_width = 10;                        % length used for binning (i.e., 10 nt x 10 nt ) 

total_sample_sets = 100;               % number of randomly sampled breakpoint sets to calculate standard deviation

% input deletion breakpoints dataset file name
% --------------------------------------------
    % see example data file in the folder for correct input format 
    % file should be in excel2003 format
    
Breakpoints_FileName = 'MouseRE_Breakpoints_Lakshmanan_et_al.xls'; 


% load mouse mtDNA hybridization partition function data
% provide the correct directory where the file is saved in the computer

load M://'Mixture Model'/Dec28_2020_Data_MixtureModel/Mouse/Mouse_Complete_MD_Propensity_DS.mat 

% output file name
% user can specify their desired filename for "Dataset_Name"
% no need to specify file format

Dataset_Name = [ 'MouseRE','-W',num2str(window_length_val),'-BW', num2str(bin_width),'-SS',num2str(total_sample_sets)]; % name includes details about dataset, window length, bin size, number of samples 

% ----------------------------------------------------------------------------------------------------------------------------------------------------------------


%% Data Analysis 
% do not modify  


Standard_BPdata = [];                       Standard_BPdata = xlsread(Breakpoints_FileName);

% redundant breakpoints filter

NonRedundant_BPdata = [];                   NonRedundant_BPdata = unique(Standard_BPdata, 'rows');

if( ( length(Standard_BPdata(:,1)) - length(NonRedundant_BPdata(:,1)) ) > 0 )

    fprintf('\n number of redundant breakpoints filtered = %d \n', ( length(Standard_BPdata(:,1)) - length(NonRedundant_BPdata(:,1)) ) )

end


% filtering breakpoints present outside two mtDNA regions  

Feasible_BP_Indices = [];                   Feasible_BP_Indices = find( ( NonRedundant_BPdata(:,2) > (Ori_LBR_start_position+window_length_val) ) & ( NonRedundant_BPdata(:,2) < (Ori_LBR_stop_position-window_length_val) ) & ( NonRedundant_BPdata(:,3) > (Ori_RBR_start_position+window_length_val) ) & ( NonRedundant_BPdata(:,3) < (Ori_RBR_stop_position-window_length_val) ) );

Filtered_Breakpoints = [];                  Filtered_Breakpoints = NonRedundant_BPdata(Feasible_BP_Indices,:);


if( ( length(NonRedundant_BPdata(:,1)) - length(Filtered_Breakpoints(:,1)) ) > 0 )
    
    fprintf('\n number of breakpoints (located outside the regions of analysis) filtered = %d \n', ( length(NonRedundant_BPdata(:,1)) - length(Filtered_Breakpoints(:,1)) ) )
    
end


% calculating the boundaries of the analysis window around breakpoints

LBR_LB = 0;                                 LBR_LB = Filtered_Breakpoints(:,2) - window_length_val;

LBR_UB = 0;                                 LBR_UB = Filtered_Breakpoints(:,2) + window_length_val;

RBR_LB = 0;                                 RBR_LB = Filtered_Breakpoints(:,3) - window_length_val;

RBR_UB = 0;                                 RBR_UB = Filtered_Breakpoints(:,3) + window_length_val;

Breakpoints_Data_Matrix = [];               Breakpoints_Data_Matrix = [ Filtered_Breakpoints, LBR_LB, LBR_UB, RBR_LB, RBR_UB ];

total_reported_deletions = 0;               total_reported_deletions = length( Breakpoints_Data_Matrix(:,1) );



% Grouping of components

Groups_DS = struct( 'Group_Number',{}, 'Components_Vector',{}, 'Overall_Component_Indices',{} );

Groups_DS(1).Group_Number = 1;
Groups_DS(1).Components_Vector = [5];
Groups_DS(1).Overall_Component_Indices = [1,2];

Groups_DS(2).Group_Number = 2;
Groups_DS(2).Components_Vector = [10,15,20,25];
Groups_DS(2).Overall_Component_Indices = [3,4,5,6];

Groups_DS(3).Group_Number = 3;
Groups_DS(3).Components_Vector = [ 50,75,100 ];
Groups_DS(3).Overall_Component_Indices = [7,8,9];

total_groups = length( Groups_DS );

MD_Components_Combination = [];

for group_number = 1:1:total_groups
    
    MD_Components_Combination = [ MD_Components_Combination, Groups_DS(group_number).Components_Vector ];   
    
end

total_MD_components = length( MD_Components_Combination );

MI_component_status = true ;





total_mixture_components = 0; 

Labels_CA = { };

if( MI_component_status == true )

    Labels_CA = [ Labels_CA, '0 bp' ];

    total_mixture_components = total_MD_components + 1;

else

    total_mixture_components = total_MD_components;

end

Labels_CA = [ Labels_CA, '5bp','10bp','15bp','20bp','25bp','50bp','75bp','100bp' ];


BoundaryValues_Matrix = Breakpoints_Data_Matrix(:,4:7);

Union_Matrix = union(BoundaryValues_Matrix,BoundaryValues_Matrix,'rows');

total_boundary_combinations = length( Union_Matrix(:,1) );

Bounds_DataStructure = struct('Boundary_Count',{}, 'BoundaryValues',{}, 'MI_PDF_Unbinned',{}, 'MD_DS', {}, 'Bin_Bounds',{}, 'Binned_PDF_DS',{}, 'Deletions_Count', {});

Overall_BP_Probability_DataStructure = struct( 'Component_Number',{}, 'Component_ID',{}, 'Probability_Vector',{} );

Breakpoints_ProbValues_Matrix = zeros(total_reported_deletions,total_mixture_components);


for boundary_combination_loop = 1:1:total_boundary_combinations  
    
        
    Bounds_DataStructure(boundary_combination_loop).Boundary_Count = boundary_combination_loop;
    
    Bounds_DataStructure(boundary_combination_loop).BoundaryValues = Union_Matrix(boundary_combination_loop,:);
    
    
    
    LBR_start_position = 0;         LBR_start_position = Union_Matrix(boundary_combination_loop,1);
    
    LBR_stop_position = 0;          LBR_stop_position = Union_Matrix(boundary_combination_loop,2);
    
    RBR_start_position = 0;         RBR_start_position = Union_Matrix(boundary_combination_loop,3);
    
    RBR_stop_position = 0;          RBR_stop_position = Union_Matrix(boundary_combination_loop,4);

    % MI component PDF calculation
    
    LBR_length = 0;                     LBR_length = LBR_stop_position - LBR_start_position + 1;

    RBR_length = 0;                     RBR_length = RBR_stop_position - RBR_start_position + 1;
    
    MI_probability_value_per_bp = 0;    MI_probability_value_per_bp = 1 / ( LBR_length * RBR_length );
    
    MI_Unbinned_PDF_Marix = [];         MI_Unbinned_PDF_Marix = MI_probability_value_per_bp * ones( LBR_length, RBR_length );
    
    Bounds_DataStructure(boundary_combination_loop).MI_PDF_Unbinned = MI_Unbinned_PDF_Marix;
    
    % MD components PDFs
    
    MDcomponents_DS = struct( 'Component_ID',{}, 'Component_PDF_Matrix',{} );
    
    for MD_component_number = 1:1:total_MD_components
    
        component_ID_val = '';              component_ID_val = [ num2str(MD_Components_Combination(MD_component_number)), 'bp'];
        
        MDcomponents_DS(MD_component_number).Component_ID = component_ID_val;

        component_index_in_reference = 0;
    
        for n_ref = 1:1:length(Mouse_Complete_PropMatrix_DS)
    
            if( strcmp( component_ID_val, Mouse_Complete_PropMatrix_DS(n_ref).MD_Component_ID ) == 1 )
        
                component_index_in_reference = n_ref; break;
            
            end        
        end

        F_MDcomponent_Propensity_Matrix = [];   F_MDcomponent_Propensity_Matrix = Mouse_Complete_PropMatrix_DS(component_index_in_reference).MD_Component_AvgPF_PropensityMatrix;   % select between AvgPF and MaxPF in this step
        
        row_lb = 0;                             row_lb = LBR_start_position - Ori_LBR_start_position + 1;
        
        row_ub = 0;                             row_ub = LBR_stop_position - Ori_LBR_start_position + 1;
        
        col_lb = 0;                             col_lb = RBR_start_position - Ori_RBR_start_position + 1;
    
        col_ub = 0;                             col_ub = RBR_stop_position - Ori_RBR_start_position + 1;
    
        MDcomponent_Propensity_Matrix = [];     MDcomponent_Propensity_Matrix = F_MDcomponent_Propensity_Matrix( (row_lb:row_ub),(col_lb:col_ub) );
        
        normalization_factor = 0;               normalization_factor = sum( sum(MDcomponent_Propensity_Matrix) );
    
        MDcomponents_DS(MD_component_number).Component_PDF_Matrix = MDcomponent_Propensity_Matrix / normalization_factor ;
        
    end
    
    Bounds_DataStructure(boundary_combination_loop).MD_DS = MDcomponents_DS;
    
    % Bin Boundaries Specification
    
    R_bin_count = 0;

    RBR_Bin_Boundary_Matrix = [];

    for RBR_loop = RBR_start_position:bin_width:(RBR_stop_position-bin_width)

        R_bin_count = R_bin_count + 1;

        R_bin_lb = 0;    R_bin_lb = RBR_loop;

        R_bin_ub = 0;    R_bin_ub = R_bin_lb + (bin_width-1);

        RBR_Bin_Boundary_Matrix = [RBR_Bin_Boundary_Matrix; R_bin_lb, R_bin_ub];            

    end

    %-------
    if( ( RBR_stop_position - RBR_Bin_Boundary_Matrix(end,2)) + 1 >= 10 )

        R_bin_count = R_bin_count + 1;

        R_bin_lb = 0;    R_bin_lb = RBR_Bin_Boundary_Matrix(end,2) + 1;

        R_bin_ub = 0;    R_bin_ub = RBR_stop_position;

        RBR_Bin_Boundary_Matrix = [RBR_Bin_Boundary_Matrix; R_bin_lb, R_bin_ub]; 

    else

        RBR_Bin_Boundary_Matrix(end,2) = RBR_stop_position;

    end


    % -------

    %RBR_Bin_Boundary_Matrix(end-2:end,:)
    
    Bin_Boundary_Matrix = [];

    for LBR_loop = LBR_start_position:bin_width:(LBR_stop_position-bin_width)

        L_bin_lb = 0;            L_bin_lb =  LBR_loop;

        L_bin_ub = 0;            L_bin_ub =  L_bin_lb + (bin_width-1);

        Local_BB_Matrix = [];    Local_BB_Matrix = [ [L_bin_lb * ones(R_bin_count,1)], [L_bin_ub * ones(R_bin_count,1)], RBR_Bin_Boundary_Matrix ];

        Bin_Boundary_Matrix = [ Bin_Boundary_Matrix; Local_BB_Matrix];

    end

    %---
    if( (LBR_stop_position - Bin_Boundary_Matrix(end,2)) + 1 >= 10  )

        L_bin_lb = 0;            L_bin_lb = Bin_Boundary_Matrix(end,2) + 1;

        L_bin_ub = 0;            L_bin_ub = LBR_stop_position;

        Local_BB_Matrix = [];    Local_BB_Matrix = [ [L_bin_lb * ones(R_bin_count,1)], [L_bin_ub * ones(R_bin_count,1)], RBR_Bin_Boundary_Matrix ];

        Bin_Boundary_Matrix = [ Bin_Boundary_Matrix; Local_BB_Matrix];

    else

        Bin_Boundary_Matrix( ((end - R_bin_count)+1):(end), 2 ) = LBR_stop_position;

    end

    Bounds_DataStructure(boundary_combination_loop).Bin_Bounds = Bin_Boundary_Matrix;

    %---
    % Bin_Boundary_Matrix( ((end - R_bin_count)+1):(end), 2 ) = LBR_stop_position;

    % Bin_Boundary_Matrix(((end - R_bin_count)):((end - R_bin_count)+1),:)

    % Bin_Boundary_Matrix( (end-2):(end),: )
    
    Components_binnedPDF_DS = struct( 'Component_ID', {}, 'Component_Bins_PDF_Vector', {}, 'Component_Bins_CDF_Vector',{} );

    Components_binnedPDF_DS = Binned_PDF_calculation_function( Bin_Boundary_Matrix, MI_component_status, MI_probability_value_per_bp, MDcomponents_DS, LBR_start_position, RBR_start_position);  
    
    Bounds_DataStructure(boundary_combination_loop).Binned_PDF_DS = Components_binnedPDF_DS;
      
    
    Boundary_Deletion_Indices = [];             Boundary_Deletion_Indices = find( ( Breakpoints_Data_Matrix(:,4) == LBR_start_position ) & ( Breakpoints_Data_Matrix(:,5) == LBR_stop_position ) & ( Breakpoints_Data_Matrix(:,6) == RBR_start_position ) & ( Breakpoints_Data_Matrix(:,7) == RBR_stop_position ) );

    total_deletions_withinBoundary = 0;         total_deletions_withinBoundary = length( Boundary_Deletion_Indices );
    
    Bounds_DataStructure(boundary_combination_loop).Deletions_Count = total_deletions_withinBoundary;
    
    Bounds_Breakpoints_Matrix = [];             Bounds_Breakpoints_Matrix = Breakpoints_Data_Matrix(Boundary_Deletion_Indices,1:3);
    
    BP_Probability_DataStructure = struct( 'Component_Number',{}, 'Component_ID',{}, 'Probability_Vector',{} );

    BP_Probability_DataStructure = Binned_BP_probability_function( Bounds_Breakpoints_Matrix, Components_binnedPDF_DS, Bin_Boundary_Matrix );    
         
    for component_number = 1:1:total_mixture_components
    
        Overall_BP_Probability_DataStructure(component_number).Component_Number = component_number;
    
        Overall_BP_Probability_DataStructure(component_number).Probability_Vector = [Overall_BP_Probability_DataStructure(component_number).Probability_Vector; BP_Probability_DataStructure(component_number).Probability_Vector];
    
    end
    
          

end

% Optimisation

total_start_points_reported = 5;

Initial_Guess_Vector = (1/total_mixture_components) * ones(1,total_mixture_components);

A_ineq = [];

b_ineq = [];

A_eq = ones(1,total_mixture_components);

b_eq = 1;

Variables_LB = zeros(1,total_mixture_components);

Variables_UB = ones(1,total_mixture_components);

NonLCon = [];

Opt = optimoptions('fmincon', 'Display','none');%, 'PlotFcns', @optimplotfval );

Opt_Problem = createOptimProblem('fmincon', 'Aineq',A_ineq, 'bineq',b_ineq, 'Aeq',A_eq, 'beq',b_eq, 'lb',Variables_LB, 'ub',Variables_UB, 'objective', @(x)LnLikelihood_Calculation_Fn(x,Overall_BP_Probability_DataStructure), 'x0',Initial_Guess_Vector, 'options', Opt );

MS = MultiStart;

x_optimum = [];     fval_optimum = 0;

[x_optimum, fval_optimum, exit_flag_val, output_val, solutions_ds] = run( MS, Opt_Problem, total_start_points_reported );


% % % Data_Likelihood_Matrix = [];
% % % 
% % % Normalised_LogLikelihood_Vector = [];
% % % 
% % % [Data_Likelihood_Matrix,Normalised_LogLikelihood_Vector] = LnLikelihood_Analysis_Fn(x_optimum,Overall_BP_Probability_DataStructure);
% % % 
% % % figure;
% % % 
% % % imagesc(Data_Likelihood_Matrix); colorbar;
% % % 
% % % ax1 = gca;
% % % 
% % % ax1.XTick = [1:1:total_mixture_components];
% % % 
% % % ax1.XTickLabel = Labels_CA;
% % % 
% % % title( Dataset_Name,'FontName','Arial','FontSize',16,'FontWeight','demi');
% % % 
% % % set(gcf,'Position',[100,100,1000,600]); grid on; 
% % % 
% % % saveas( gcf, [Dataset_Name,'_HeatMap'], 'fig'); 
% % % 
% % % saveas( gcf, [Dataset_Name,'_HeatMap'], 'jpg'); 




% sample generation

total_samples = 0;                          total_samples = total_reported_deletions;

Preset_Component_Fraction = [];             Preset_Component_Fraction = x_optimum;

Cumulative_Space_Matrix = zeros(total_mixture_components,2);  

Cumulative_Space_Matrix(1,2) = Preset_Component_Fraction(1);

for component_loop = 2:1:total_mixture_components

    Cumulative_Space_Matrix(component_loop,1) = Cumulative_Space_Matrix(component_loop-1,2);

    Cumulative_Space_Matrix(component_loop,2) = Cumulative_Space_Matrix(component_loop,1) + Preset_Component_Fraction(component_loop);

end

Prediction_Results = zeros(total_sample_sets,total_mixture_components); 

Prediction_LikelihoodValues = zeros(total_sample_sets,1);



URN_stream1 = RandStream('mt19937ar','Seed',main_seed_value);

URN_stream_ComponentSelection = RandStream('mt19937ar','Seed',(main_seed_value+10000));

total_start_points_random = 5;

for set_number = 1:1:total_sample_sets
    
    overall_sample_count = 1;
    
    Sampled_Breakpoints_Data = []; 
    
    for boundary_combination_loop = 1:1:total_boundary_combinations  
        
        
        Bounds_Vector = [];             Bounds_Vector = Bounds_DataStructure(boundary_combination_loop).BoundaryValues;
    
        LBR_start_position = 0;         LBR_start_position = Bounds_Vector(1,1);
    
        LBR_stop_position = 0;          LBR_stop_position = Bounds_Vector(1,2);

        RBR_start_position = 0;         RBR_start_position = Bounds_Vector(1,3);

        RBR_stop_position = 0;          RBR_stop_position = Bounds_Vector(1,4);        
        
        MI_Unbinned_PDF_Marix = [];         MI_Unbinned_PDF_Marix = Bounds_DataStructure(boundary_combination_loop).MI_PDF_Unbinned;
        
        clear MDcomponents_DS; 
        
        MDcomponents_DS = Bounds_DataStructure(boundary_combination_loop).MD_DS;
    
        bounds_total_sample_count = 0;      bounds_total_sample_count = Bounds_DataStructure(boundary_combination_loop).Deletions_Count;  
    
        %boundary_combination_loop,bounds_sample_loop
        for bounds_sample_loop = 1:1:bounds_total_sample_count
        
            if( ( boundary_combination_loop == 1 ) & ( bounds_sample_loop == 1 ) ) % the very first sample

                r_l1 = rand(URN_stream_ComponentSelection,1,1);
                
                component_index = 0;        component_index = find( (r_l1 >= Cumulative_Space_Matrix(:,1)) & (r_l1 < Cumulative_Space_Matrix(:,2)) );

                LBR_length = 0;                     LBR_length = LBR_stop_position - LBR_start_position + 1;

                RBR_length = 0;                     RBR_length = RBR_stop_position - RBR_start_position + 1;
        
                Local_Component_Unbinned_PDF_Matrix = [];

                if( component_index == 1 )
                    
                    Local_Component_Unbinned_PDF_Matrix = MI_Unbinned_PDF_Marix;    % MI component                

                elseif ( component_index > 1 )
                    
                    MD_component_index = component_index - 1;
                    
                    Local_Component_Unbinned_PDF_Matrix = MDcomponents_DS(MD_component_index).Component_PDF_Matrix;
                    
                else
                    
                    fprintf('\n ERROR IN COMPONENT INDEX WHILE SAMPLING RANDOM DELETIONS\n'); 
                    
                end
                
                Local_PDF_Tr_ColumnVector = [];         Local_PDF_Tr_ColumnVector = reshape(Local_Component_Unbinned_PDF_Matrix',[],1);

                Local_Cuml_ColumnVector = [];           Local_Cuml_ColumnVector = cumsum(Local_PDF_Tr_ColumnVector);

                local_leftBP_position = 0;

                local_rightBP_position = 0;

                r_l2 = rand(URN_stream1,1,1);

                local_position_index = 0;               local_position_index = min( find( r_l2 <= Local_Cuml_ColumnVector ) );

                local_PDF_Matrix_col_id = 0;             

                local_rem_value = 0;                    local_rem_value =  rem(local_position_index, RBR_length);          

                if( local_rem_value == 0 )

                    local_PDF_Matrix_col_id = RBR_length;

                    local_PDF_Matrix_row_id = 0;            local_PDF_Matrix_row_id = floor( local_position_index /   RBR_length );

                else

                    local_PDF_Matrix_col_id = local_rem_value;

                    local_PDF_Matrix_row_id = 0;            local_PDF_Matrix_row_id = (floor( local_position_index /   RBR_length )) + 1;

                end

                local_leftBP_position = 0;              local_leftBP_position = (LBR_start_position - 1) + local_PDF_Matrix_row_id;

                local_rightBP_position = 0;             local_rightBP_position = (RBR_start_position - 1) + local_PDF_Matrix_col_id;

                Sampled_Breakpoints_Data = [Sampled_Breakpoints_Data; overall_sample_count, local_leftBP_position, local_rightBP_position, Bounds_Vector,boundary_combination_loop,bounds_sample_loop,component_index ]; 
                
                overall_sample_count = overall_sample_count + 1;
                
              
            else
              
                r_l1 = rand(URN_stream_ComponentSelection,1,1);

                component_index = 0;        component_index = find( (r_l1 >= Cumulative_Space_Matrix(:,1)) & (r_l1 < Cumulative_Space_Matrix(:,2)) );

                LBR_length = 0;                     LBR_length = LBR_stop_position - LBR_start_position + 1;

                RBR_length = 0;                     RBR_length = RBR_stop_position - RBR_start_position + 1;

                Local_Component_Unbinned_PDF_Matrix = [];
                
                if( component_index == 1 )
                    
                    Local_Component_Unbinned_PDF_Matrix = MI_Unbinned_PDF_Marix;    % MI component                

                elseif ( component_index > 1 )
                    
                    MD_component_index = component_index - 1;
                    
                    Local_Component_Unbinned_PDF_Matrix = MDcomponents_DS(MD_component_index).Component_PDF_Matrix;
                    
                else
                    
                    fprintf('\n ERROR IN COMPONENT INDEX WHILE SAMPLING RANDOM DELETIONS\n'); 
                    
                end
                
                Local_PDF_Tr_ColumnVector = [];         Local_PDF_Tr_ColumnVector = reshape(Local_Component_Unbinned_PDF_Matrix',[],1);

                Local_Cuml_ColumnVector = [];           Local_Cuml_ColumnVector = cumsum(Local_PDF_Tr_ColumnVector);

                local_leftBP_position = 0;

                local_rightBP_position = 0;

                redundancy_index = 99;

                while( (redundancy_index == 0) || (redundancy_index == 99))

                    r_l2 = rand(URN_stream1,1,1);

                    local_position_index = 0;               local_position_index = min( find( r_l2 <= Local_Cuml_ColumnVector ) );

                    local_PDF_Matrix_col_id = 0;             

                    local_rem_value = 0;                    local_rem_value =  rem(local_position_index, RBR_length);         

                    if( local_rem_value == 0 )

                        local_PDF_Matrix_col_id = RBR_length;

                        local_PDF_Matrix_row_id = 0;            local_PDF_Matrix_row_id = floor( local_position_index /   RBR_length );

                    else

                        local_PDF_Matrix_col_id = local_rem_value;

                        local_PDF_Matrix_row_id = 0;            local_PDF_Matrix_row_id = (floor( local_position_index /   RBR_length )) + 1;

                    end


                    local_leftBP_position = 0;              local_leftBP_position = (LBR_start_position - 1) + local_PDF_Matrix_row_id;

                    local_rightBP_position = 0;             local_rightBP_position = (RBR_start_position - 1) + local_PDF_Matrix_col_id;

                    redundancy_index = isempty(  find( (Sampled_Breakpoints_Data(:,2) == local_leftBP_position) & (Sampled_Breakpoints_Data(:,3) == local_rightBP_position) )  );

                end

                Sampled_Breakpoints_Data = [Sampled_Breakpoints_Data; overall_sample_count, local_leftBP_position, local_rightBP_position, Bounds_Vector,boundary_combination_loop,bounds_sample_loop,component_index ];
                
                overall_sample_count = overall_sample_count + 1;
                
                               
            end
                
            
        end
    
    
    end
       
    Sampled_BP_Probability_DataStructure = struct( 'Component_Number',{}, 'Component_ID',{}, 'Probability_Vector',{} );
    
    Sampled_BP_Probability_DataStructure = SampledBP_Probability_Calculation_Function(Sampled_Breakpoints_Data, Bounds_DataStructure, total_mixture_components);
    
    x_L = []; fval_L = 0;

    Opt_Problem_S = createOptimProblem('fmincon', 'Aineq',A_ineq, 'bineq',b_ineq, 'Aeq',A_eq, 'beq',b_eq, 'lb',Variables_LB, 'ub',Variables_UB, 'objective', @(x)LnLikelihood_Calculation_Fn(x,Sampled_BP_Probability_DataStructure), 'x0',Initial_Guess_Vector, 'options', Opt );

    [x_L,fval_L,exitflag_L,output_L,solutions_L] = run( MS, Opt_Problem_S, total_start_points_random );

   % [x_L,fval_L,exitflag_L,output_L,solutions_L] = fmincon( @(x)LnLikelihood_Calculation_Fn(x,Sampled_BP_Probability_DataStructure), Initial_Guess_Vector, A_ineq, b_ineq, A_eq, b_eq, Variables_LB, Variables_UB, NonLCon, Opt);  

    Prediction_Results(set_number,:) =  x_L;

    Prediction_LikelihoodValues(set_number,:) =  -fval_L ;
    
       
    
    
    
end


Predictions_MeanValue_Vector = [];

Predictions_StdValue_Vector = [];

for component_number = 1:1:total_mixture_components

    Predictions_MeanValue_Vector = [ Predictions_MeanValue_Vector; mean( Prediction_Results(:,component_number) ) ];

    Predictions_StdValue_Vector = [ Predictions_StdValue_Vector; std( Prediction_Results(:,component_number) ) ];       

end

xlswrite( [Dataset_Name, '_ComponentPredictions'], [x_optimum;(Predictions_StdValue_Vector')] );

xlswrite( [Dataset_Name, '_FullResults'], [ [x_optimum, -fval_L];[Prediction_Results,Prediction_LikelihoodValues]] );



Samples_Average_Lk_Vector = [];                     Samples_Average_Lk_Vector = mean( Prediction_LikelihoodValues );

Samples_Deviation_Lk_Vector = [];                   Samples_Deviation_Lk_Vector = std( Prediction_LikelihoodValues );

zt_h = 99;      zt_p = 99;

[zt_h,zt_p] = ztest( -fval_optimum, Samples_Average_Lk_Vector, Samples_Deviation_Lk_Vector, 'Tail', 'left' );

y_axis_min = 0;             y_axis_min = min( [Prediction_LikelihoodValues;  -fval_optimum ] );            y_axis_min = y_axis_min + ( 0.1 * y_axis_min );

y_axis_max = 0;             y_axis_max = max( [Prediction_LikelihoodValues;  -fval_optimum ] );            y_axis_max = y_axis_max - (0.1 * y_axis_max );

figure;

boxplot( Prediction_LikelihoodValues, 'color','k', 'Symbol', 'k+' ); hold on;

plot( 1, -fval_optimum, 'Marker','p', 'MarkerSize', 10, 'Color','k', 'LineWidth',2, 'LineStyle',':' ); hold off;

ylim([y_axis_min, y_axis_max]);  

ylabel('Ln (likelihood)','FontName','Arial','FontSize',16,'FontWeight','demi');

title( Dataset_Name,'FontName','Arial','FontSize',16,'FontWeight','demi');

grid on;

text( 1.1, median(Prediction_LikelihoodValues) + 4, ['Z-test'] ,'FontName','Arial','FontSize',16,'FontWeight','demi');

text( 1.1, median(Prediction_LikelihoodValues), ['p-value: ', num2str(zt_p)] ,'FontName','Arial','FontSize',16,'FontWeight','demi');

set(gca,'FontName','Arial', 'FontSize',16, 'FontWeight','demi', 'XTick',[1], 'XTickLabel', {});  

saveas(gcf, [Dataset_Name,'_Likelihood'], 'epsc'); 
        
saveas(gcf, [Dataset_Name, '_Likelihood'], 'fig');  clear gcf;



Grouped_SamplePredictions = zeros(total_sample_sets,total_groups);

for set_number = 1:1:total_sample_sets

    for group_number = 1:1:total_groups

        Group_Indices_Vector = [];              Group_Indices_Vector = Groups_DS(group_number).Overall_Component_Indices;
        
        Grouped_SamplePredictions(set_number,group_number) = sum( Prediction_Results(set_number, Group_Indices_Vector) );       
        
    end

end

GroupPredictions_StdValue_Vector = zeros(total_groups,1);

for group_number = 1:1:total_groups
    
   GroupPredictions_StdValue_Vector(group_number,1) = std( Grouped_SamplePredictions(:,group_number) );    
    
end

ReportedData_GroupPredictions = zeros(total_groups,1);

for group_number = 1:1:total_groups
    
    Group_Indices_Vector = [];              Group_Indices_Vector = Groups_DS(group_number).Overall_Component_Indices;
    
    ReportedData_GroupPredictions(group_number,1) = sum( x_optimum(Group_Indices_Vector) );
    
end

xlswrite( [Dataset_Name, '_GroupedPredictions'], [(ReportedData_GroupPredictions'); (GroupPredictions_StdValue_Vector')] );


Color_Triples_Matrix = [    0.50, 0.50, 0.00;
                            0.00, 0.50, 0.50;
                            0.50, 0.00, 0.50;
                            0.75, 0.75, 0.75;
                            1.00, 0.00, 0.00 ];  %  [ olive, teal, purple,silver,red ];


figure;

additional_width = 1;

for group_number = 1:1:total_groups
    
% % %     if(group_number == 1)
% % %         
% % %        x_lb = -additional_width; 
% % %        
% % %        if(isempty(Groups_DS(group_number).Components_Vector) == 1)
% % % 
% % %             x_mid_value = 0;
% % %         
% % %        else
% % %         
% % %             x_mid_value = 2.5;
% % %     
% % %        end
% % %        
% % %        
% % %         
% % %     else
% % %         
% % %         x_lb = ( Groups_DS(group_number).Components_Vector(1) ) - additional_width;
% % %         
% % %         x_mid_value = ( ( Groups_DS(group_number).Components_Vector(1) ) + (Groups_DS(group_number).Components_Vector(end)) ) / 2;
% % %         
% % %     end
% % %     
% % %     if(group_number == total_groups)
% % %         
% % %         x_mid_value = 45;
% % %         
% % %     end
% % %     
% % %     if(isempty(Groups_DS(group_number).Components_Vector) == 1)
% % % 
% % %         x_ub = additional_width;
% % %         
% % %     else
% % %         
% % %         x_ub = (Groups_DS(group_number).Components_Vector(end)) + additional_width;
% % %     
% % %     end


    if(group_number == 1)
    
        x_lb = 0;
        
        x_ub = 5;
        
        x_mid_value = 2.5;
        
    end
    
    if(group_number == 2)
       
        x_lb = 10;
        
        x_ub = 25;
        
        x_mid_value = 17.5;
    end
    
    if(group_number == 3)
    
        x_lb = 30;
        
        x_ub = 40;
    
        x_mid_value = 35;
    
    end
        
    y_ub = ReportedData_GroupPredictions(group_number,1);
    
    X_Vector = [];          X_Vector = [ x_lb, x_lb, x_ub, x_ub, x_lb ];

    Y_Vector = [];          Y_Vector = [ 0, y_ub, y_ub, 0, 0 ];
    
    group_SD_value = 0;     group_SD_value = GroupPredictions_StdValue_Vector(group_number,1);

    patch( X_Vector,Y_Vector, Color_Triples_Matrix(group_number,:), 'EdgeColor','none'); hold on;
    
    plot( (x_mid_value * ones(3,1)), [(y_ub-group_SD_value),(y_ub),(y_ub+group_SD_value)], 'LineStyle','-','LineWidth',3, 'Marker','none', 'Color',Color_Triples_Matrix(group_number,:)) ; hold on;
    
    

end

alpha(0.25);

plot([0:5:40], x_optimum, 'Marker','o',  'Color','k',  'LineWidth',1,  'LineStyle',':' );  hold on;

errorbar( [0:5:40], x_optimum, Predictions_StdValue_Vector, 'Marker', 'none', 'Color','k', 'LineStyle','none', 'LineWidth',1 ); hold off;

set(gca,'XTick', [0:5:40], 'XTickLabel', [0,MD_Components_Combination] , 'FontName','Arial', 'FontSize',14, 'FontWeight','demi');  ylim([0,1]); xlim([-5,45]);

%xlabel('components', 'FontName','Arial', 'FontSize',16, 'FontWeight','demi'); 

ylabel('predicted fraction', 'FontName','Arial', 'FontSize',16, 'FontWeight','demi');

title([Dataset_Name, ' (n: ', num2str(total_reported_deletions), ' deletions)'], 'FontName','Arial', 'FontSize',16, 'FontWeight','demi');

set(gcf,'Position',[100,100,1000,600]); grid on; 

saveas( gcf, [Dataset_Name,'_CombinedPredictions'], 'fig'); 
 
saveas( gcf, [Dataset_Name,'_CombinedPredictions'], 'epsc'); 

[ ReportedData_GroupPredictions, GroupPredictions_StdValue_Vector ]

toc




