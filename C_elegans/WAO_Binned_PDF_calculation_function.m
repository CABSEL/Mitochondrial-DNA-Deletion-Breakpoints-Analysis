function[Components_binnedPDF_DS] = WAO_Binned_PDF_calculation_function( Bin_Boundary_Matrix, MI_component_status, MI_Probability_Matrix, MDcomponents_DS, LBR_start_position, RBR_start_position)

    Components_binnedPDF_DS = struct( 'Component_ID', {}, 'Component_Bins_PDF_Vector', {}, 'Component_Bins_CDF_Vector',{} );
    
    total_bins = 0;                         total_bins = length( Bin_Boundary_Matrix(:,1) );
    
    total_component_count = 0;
    
    
    
    if( MI_component_status == true )
        
        total_component_count = total_component_count + 1;
        
        Components_binnedPDF_DS(total_component_count).Component_ID = 'MI';        
        
        Component_PDFmatrix = [];               Component_PDFmatrix = MI_Probability_Matrix;
        
        Component_Binned_PDF_Vector = zeros(total_bins,1);
        
        for bin_number = 1:1:total_bins
        
            L_bin_lb = 0;                       L_bin_lb = Bin_Boundary_Matrix(bin_number,1);           row_lb = 0;         row_lb = L_bin_lb - (LBR_start_position - 1);
        
            L_bin_ub = 0;                       L_bin_ub = Bin_Boundary_Matrix(bin_number,2);           row_ub = 0;         row_ub = L_bin_ub - (LBR_start_position - 1);
        
            R_bin_lb = 0;                       R_bin_lb = Bin_Boundary_Matrix(bin_number,3);           col_lb = 0;         col_lb = R_bin_lb - (RBR_start_position - 1);
            
            R_bin_ub = 0;                       R_bin_ub = Bin_Boundary_Matrix(bin_number,4);           col_ub = 0;         col_ub = R_bin_ub - (RBR_start_position - 1);
            
            Bin_Submatrix = [];                 Bin_Submatrix = Component_PDFmatrix( (row_lb:row_ub),(col_lb:col_ub) );   
            
            bin_p_value = 0;                    bin_p_value = sum( sum(Bin_Submatrix) );
        
            Component_Binned_PDF_Vector(bin_number,1) = bin_p_value;
        
        end
        
        Components_binnedPDF_DS(total_component_count).Component_Bins_PDF_Vector = Component_Binned_PDF_Vector;
        
        Component_Binned_CDF_Vector = [];       Component_Binned_CDF_Vector = cumsum(Component_Binned_PDF_Vector);
        
        Components_binnedPDF_DS(total_component_count).Component_Bins_CDF_Vector = Component_Binned_CDF_Vector;
        
        
        
    end  
    
    
    total_MD_components = 0;            total_MD_components = length( MDcomponents_DS );  % % % % %  MDcomponents_DS  struct( 'Component_ID',{}, 'Component_PDF_Matrix',{} );
    
    for MD_component_count = 1:1:total_MD_components

        total_component_count = total_component_count + 1;
        
        Components_binnedPDF_DS(total_component_count).Component_ID = MDcomponents_DS(MD_component_count).Component_ID;
        
        Component_PDFmatrix = [];               Component_PDFmatrix = MDcomponents_DS(MD_component_count).Component_PDF_Matrix;
        
       % sum( sum(Component_PDFmatrix) )
        
        Component_Binned_PDF_Vector = zeros(total_bins,1);
        
        for bin_number = 1:1:total_bins
        
            L_bin_lb = 0;                       L_bin_lb = Bin_Boundary_Matrix(bin_number,1);           row_lb = 0;         row_lb = L_bin_lb - (LBR_start_position - 1);
        
            L_bin_ub = 0;                       L_bin_ub = Bin_Boundary_Matrix(bin_number,2);           row_ub = 0;         row_ub = L_bin_ub - (LBR_start_position - 1);
        
            R_bin_lb = 0;                       R_bin_lb = Bin_Boundary_Matrix(bin_number,3);           col_lb = 0;         col_lb = R_bin_lb - (RBR_start_position - 1);
            
            R_bin_ub = 0;                       R_bin_ub = Bin_Boundary_Matrix(bin_number,4);           col_ub = 0;         col_ub = R_bin_ub - (RBR_start_position - 1);
            
            Bin_Submatrix = [];                 Bin_Submatrix = Component_PDFmatrix( (row_lb:row_ub),(col_lb:col_ub) );   
            
            bin_p_value = 0;                    bin_p_value = sum( sum(Bin_Submatrix) );
        
            Component_Binned_PDF_Vector(bin_number,1) = bin_p_value;
        
        end
        
        Components_binnedPDF_DS(total_component_count).Component_Bins_PDF_Vector = Component_Binned_PDF_Vector;
        
        Component_Binned_CDF_Vector = [];       Component_Binned_CDF_Vector = cumsum(Component_Binned_PDF_Vector);
        
        Components_binnedPDF_DS(total_component_count).Component_Bins_CDF_Vector = Component_Binned_CDF_Vector;
        
    end


end