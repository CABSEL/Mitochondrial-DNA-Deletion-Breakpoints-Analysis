function[BP_Probability_DataStructure] = Binned_BP_probability_function( Deletion_Breakpoints_Matrix, Components_binnedPDF_DS, Bin_Boundary_Matrix )

    BP_Probability_DataStructure = struct( 'Component_Number',{}, 'Component_ID',{}, 'Probability_Vector',{} );
    
    total_deletions_count = 0;                  total_deletions_count = length( Deletion_Breakpoints_Matrix(:,1) );

    total_mixture_components = 0;               total_mixture_components = length( Components_binnedPDF_DS );

    for component_number = 1:1:total_mixture_components
        
        %component_number
        
        BP_Probability_DataStructure(component_number).Component_ID = Components_binnedPDF_DS(component_number).Component_ID;
        
        Component_BP_Probability_Vector = zeros(total_deletions_count,1);
        
        Component_PDF_Vector = [];              Component_PDF_Vector = Components_binnedPDF_DS(component_number).Component_Bins_PDF_Vector;

        for deletion_number = 1:1:total_deletions_count
            
                       
            %[ deletion_number, Deletion_Breakpoints_Matrix(deletion_number,2),  Deletion_Breakpoints_Matrix(deletion_number,3) ]

            bin_index = 0;                       bin_index = find( ( Deletion_Breakpoints_Matrix(deletion_number,2) >=  Bin_Boundary_Matrix(:,1) ) & ( Deletion_Breakpoints_Matrix(deletion_number,2) <= Bin_Boundary_Matrix(:,2) ) & ( Deletion_Breakpoints_Matrix(deletion_number,3) >= Bin_Boundary_Matrix(:,3) ) & ( Deletion_Breakpoints_Matrix(deletion_number,3) <= Bin_Boundary_Matrix(:,4) ) );
            
            Component_BP_Probability_Vector(deletion_number,1) = Component_PDF_Vector(bin_index,1);
            
        end
        
        BP_Probability_DataStructure(component_number).Probability_Vector = Component_BP_Probability_Vector;

    end



end