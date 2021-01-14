% Automated filtering of redundant deletions
% primer positions can be present after 3rd column
function[NonRedundant_BP_Matrix] = Redundant_Breakpoints_Filtering_Function(Breakpoints_Matrix)


    total_initial_deletions = 0;            total_initial_deletions = length( Breakpoints_Matrix(:,1) );    

    NonRedundant_BP_Matrix = [];
    
    total_NR_deletions = 0;

    for deletion_count = 1:1:total_initial_deletions

        if( deletion_count == 1)
            
            total_NR_deletions = total_NR_deletions + 1;

            NonRedundant_BP_Matrix = [ NonRedundant_BP_Matrix; [total_NR_deletions, Breakpoints_Matrix(deletion_count,2:3)] ];

        else

            comp_leftBP = 0;            comp_leftBP = Breakpoints_Matrix(deletion_count,2);

            comp_rightBP = 0;           comp_rightBP = Breakpoints_Matrix(deletion_count,3);

            Comparison_Indices = [];    Comparison_Indices = find( ( NonRedundant_BP_Matrix(:,2) == comp_leftBP ) & ( NonRedundant_BP_Matrix(:,3) == comp_rightBP ) );

            if( isempty(Comparison_Indices) == 1 )
                
                total_NR_deletions = total_NR_deletions + 1;

                NonRedundant_BP_Matrix = [ NonRedundant_BP_Matrix; total_NR_deletions, Breakpoints_Matrix(deletion_count,2:3) ];          

            end

       end

    end
    
    
end