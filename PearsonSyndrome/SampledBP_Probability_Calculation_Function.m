function[Overall_BP_Probability_DataStructure] = SampledBP_Probability_Calculation_Function(Sampled_Breakpoints_Data, Bounds_DataStructure, total_mixture_components)

    
    Breakpoints_Data_Matrix = Sampled_Breakpoints_Data;

    total_boundary_combinations = length( Bounds_DataStructure );
    
    Overall_BP_Probability_DataStructure = struct( 'Component_Number',{}, 'Component_ID',{}, 'Probability_Vector',{} );
    
    for boundary_combination_loop = 1:1:total_boundary_combinations  
    
        
        LBR_start_position = 0;         LBR_start_position = Bounds_DataStructure(boundary_combination_loop).BoundaryValues(1,1);

        LBR_stop_position = 0;          LBR_stop_position = Bounds_DataStructure(boundary_combination_loop).BoundaryValues(1,2);

        RBR_start_position = 0;         RBR_start_position = Bounds_DataStructure(boundary_combination_loop).BoundaryValues(1,3);

        RBR_stop_position = 0;          RBR_stop_position = Bounds_DataStructure(boundary_combination_loop).BoundaryValues(1,4);
        
                
        
        Bin_Boundary_Matrix = Bounds_DataStructure(boundary_combination_loop).Bin_Bounds;

        Components_binnedPDF_DS = Bounds_DataStructure(boundary_combination_loop).Binned_PDF_DS;

        Boundary_Deletion_Indices = [];             Boundary_Deletion_Indices = find( ( Breakpoints_Data_Matrix(:,4) == LBR_start_position ) & ( Breakpoints_Data_Matrix(:,5) == LBR_stop_position ) & ( Breakpoints_Data_Matrix(:,6) == RBR_start_position ) & ( Breakpoints_Data_Matrix(:,7) == RBR_stop_position ) );

        Bounds_Breakpoints_Matrix = [];             Bounds_Breakpoints_Matrix = Breakpoints_Data_Matrix(Boundary_Deletion_Indices,1:3);

        BP_Probability_DataStructure = struct( 'Component_Number',{}, 'Component_ID',{}, 'Probability_Vector',{} );

        BP_Probability_DataStructure = Binned_BP_probability_function( Bounds_Breakpoints_Matrix, Components_binnedPDF_DS, Bin_Boundary_Matrix );    

        for component_number = 1:1:total_mixture_components

            Overall_BP_Probability_DataStructure(component_number).Component_Number = component_number;

            Overall_BP_Probability_DataStructure(component_number).Probability_Vector = [Overall_BP_Probability_DataStructure(component_number).Probability_Vector; BP_Probability_DataStructure(component_number).Probability_Vector];

        end     



    end



end