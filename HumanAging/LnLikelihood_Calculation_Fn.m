function f = LnLikelihood_Calculation_Fn(x,BP_Probability_DataStructure)

    total_mixture_components = 0;           total_mixture_components = length( BP_Probability_DataStructure );
    
    total_deletions = 0;                    total_deletions = length( BP_Probability_DataStructure(1).Probability_Vector );
    
    Data_Likelihood_Vector = zeros(total_deletions,1);
    
    for component_number = 1:1:total_mixture_components

        Data_Likelihood_Vector = Data_Likelihood_Vector + ( ( x(component_number) ) * ( BP_Probability_DataStructure(component_number).Probability_Vector ) ); 

    end
    
    Data_LogLikelihood_Vector = [];         Data_LogLikelihood_Vector = log( Data_Likelihood_Vector );
       
    Non_Inf_Indices = [];                   Non_Inf_Indices = find( isinf(Data_LogLikelihood_Vector) == 0 );
    
    Non_Inf_ValuesVector = [];              Non_Inf_ValuesVector = Data_LogLikelihood_Vector(Non_Inf_Indices);
    
    if( length(Non_Inf_Indices) ~= total_deletions )
        
        length(Non_Inf_Indices)
        
    end
    
    total_loglikelihood_value = 0;          total_loglikelihood_value = (sum( Non_Inf_ValuesVector )) ;% / (length(Non_Inf_Indices)) ;
    
    

    f = -total_loglikelihood_value;
    
end