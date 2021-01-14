% local bounda data generation

clear; clc; close('all');

Ori_LBR_start_position = 5000;          Ori_LBR_stop_position = 11300;

Ori_RBR_start_position = 11500;         Ori_RBR_stop_position = 16465;


%Standard_BPdata = xlsread('HumanAging_PrimerBounds_BreakpointsData.xls');

Standard_BPdata = xlsread('HumanAging_BreakpointsData.xls');

local_window_length = 100; 

output_file_name = ['HumanAging_ClassII_LocalBounds_', num2str(local_window_length), '_BreakpointsData']; % 'MouseRE_LocalBounds_250_BreakpointsData';

% local_window_length = 250; minimum_resolution = 25; 

% local_window_length = 500; minimum_resolution = 50; 


LBR_LB = Standard_BPdata(:,2) - local_window_length;

LBR_UB = Standard_BPdata(:,2) + local_window_length;

RBR_LB = Standard_BPdata(:,3) - local_window_length;

RBR_UB = Standard_BPdata(:,3) + local_window_length;


% % % LBR_LB
% % 
% % LBR_LB = floor( ( Standard_BPdata(:,2) - (local_window_length-1) )/minimum_resolution ) * minimum_resolution;
% % 
% % LBR_LB_Indices = find( LBR_LB < Ori_LBR_start_position );
% % 
% % if( isempty(LBR_LB_Indices) ~= 1 )
% % 
% %     LBR_LB(LBR_LB_Indices) = Ori_LBR_start_position;
% % 
% % end
% % 
% % % LBR_UB
% % 
% % LBR_UB = floor(( Standard_BPdata(:,2) + local_window_length )/minimum_resolution ) * minimum_resolution;
% % 
% % LBR_UB_Indices = find( LBR_UB > Ori_LBR_stop_position );
% % 
% % if( isempty(LBR_UB_Indices) ~= 1 )
% %     
% %     LBR_UB(LBR_UB_Indices) = Ori_LBR_stop_position;
% %     
% % end
% % 
% % % RBR_LB
% % 
% % RBR_LB = floor( ( Standard_BPdata(:,3) - (local_window_length-1) )/minimum_resolution ) * minimum_resolution;
% % 
% % RBR_LB_Indices = find( RBR_LB < Ori_RBR_start_position );
% % 
% % if( isempty(RBR_LB_Indices) ~= 1 )
% %     
% %     RBR_LB(RBR_LB_Indices) = Ori_RBR_start_position;
% % 
% % end
% % 
% % 
% % % RBR_UB
% % 
% % RBR_UB = floor(( Standard_BPdata(:,3) + local_window_length )/minimum_resolution ) * minimum_resolution;
% % 
% % RBR_UB_Indices = find( RBR_UB > Ori_RBR_stop_position );
% % 
% % if( isempty(RBR_UB_Indices) ~= 1 )
% %     
% %     RBR_UB(RBR_UB_Indices) = Ori_RBR_stop_position;
% %     
% % end

% Bounds_Matrix 

Bounds_Matrix = [ LBR_LB, Standard_BPdata(:,2), LBR_UB, RBR_LB, Standard_BPdata(:,3), RBR_UB ];

LocalBounds_BPdata = [ Standard_BPdata(:,1:3), LBR_LB, LBR_UB, RBR_LB, RBR_UB ];

xlswrite(output_file_name,LocalBounds_BPdata);
