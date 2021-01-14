function [Data_Likelihood_Matrix,Normalised_LogLikelihood_Vector] = LnLikelihood_Analysis_Fn(x,BP_Probability_DataStructure)

    total_mixture_components = 0;           total_mixture_components = length( BP_Probability_DataStructure );
    
    total_deletions = 0;                    total_deletions = length( BP_Probability_DataStructure(1).Probability_Vector );
    
    Data_Likelihood_Vector = zeros(total_deletions,1);
    
    Data_Likelihood_Matrix = zeros(total_deletions,total_mixture_components);
    
    for component_number = 1:1:total_mixture_components
        
        Data_Likelihood_Matrix(:,component_number) = ( ( x(component_number) ) * ( BP_Probability_DataStructure(component_number).Probability_Vector ) ); 

        Data_Likelihood_Vector = Data_Likelihood_Vector + ( ( x(component_number) ) * ( BP_Probability_DataStructure(component_number).Probability_Vector ) ); 

    end
    
    for deletion_number = 1:1:total_deletions
        
        Data_Likelihood_Matrix(deletion_number,:) = Data_Likelihood_Matrix(deletion_number,:) / (sum(Data_Likelihood_Matrix(deletion_number,:)));
        
    end
    
    Data_LogLikelihood_Vector = [];         Data_LogLikelihood_Vector = log( Data_Likelihood_Vector );
    
    Normalised_LogLikelihood_Vector = [];   Normalised_LogLikelihood_Vector = Data_LogLikelihood_Vector / (sum(Data_LogLikelihood_Vector));
       
%     Non_Inf_Indices = [];                   Non_Inf_Indices = find( isinf(Data_LogLikelihood_Vector) == 0 );
%     
%     Non_Inf_ValuesVector = [];              Non_Inf_ValuesVector = Data_LogLikelihood_Vector(Non_Inf_Indices);
%     
%     if( length(Non_Inf_Indices) ~= total_deletions )
%         
%         length(Non_Inf_Indices)
%         
%     end
%     
%     total_loglikelihood_value = 0;          total_loglikelihood_value = (sum( Non_Inf_ValuesVector )) ;% / (length(Non_Inf_Indices)) ;
%     
%     
% 
%     f = -total_loglikelihood_value;
    
end