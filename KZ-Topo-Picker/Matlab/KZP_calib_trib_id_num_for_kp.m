function [all_trib_ids, stream_id, ksn_SA_list, ks_chiplot_list,...
    m_n_basinwide_list, all_ks_trib_list] = ...
    KZP_calib_trib_id_num_for_kp(beta_all_tribs_all_basins, kp_elev_stored_lips,theta_bf,...
    ks_chiplot,ksn_SA_store)

% ^^ inputs: measured reference profile steepness for each tributary. One
% of the cell arrays the contained the measurement for one of the knickzone
% properties (used to figure out how many knickzones were selected from
% each tributary), list of m/n used for chiplot construction, ks derived
% from chiplot, ksn derived from SA plot

% outputs: list containing the tributary that each knickzone was selected
% from. list containing the basin # that each knickzone was selected from.
% list containing ksn from SA plot for each basin, list containing ks from
% chiplot for each basin (ksn if ref concavity used), list containing m/n
% ratio used for chiplot construction (either ref concavity or best fit m/n
% depending on what was specified in parameters script)
% list containing the steepness of the reference profile for each tributary

% ^^ these get organized into a table with row = knickzone # an column =
% knickzone attribute

% get the number of basins that contained knickzones
num_basins2 = length(beta_all_tribs_all_basins);


for i = 1:num_basins2 
    % for all stream basins that contained knickzones
    
    current_basin_elevation = kp_elev_stored_lips{i};
    % need the length of these to find out how many knickzones were
    % selected in each tributary
    
    beta_list_current_basin = beta_all_tribs_all_basins{i};
    % unpack the steepness of the reference profile for each tributary
    
    for j= 1:length(current_basin_elevation);
        % for the number of tributaries  in the current basin
        trib_number= j;
        % record the id tributary number
        trib_id_matrix = [];
        % create empty matrix to store tributary id numbers for each
        % knickzone in the tributary
        beta_matrix = [];
        % create empty matrix to store the steepness of the refrence
        % profile for the tributary (for each knickzone in the tributary)
        
        number_of_kps = length(current_basin_elevation{j}) ; 
        % organize cell array so each tributary records the same id number  
        % (how many knickpoints are in this current tributary)
        
        current_beta = beta_list_current_basin{j}; 
        % get the beta used for the current tributary
        
        for k = 1:number_of_kps ;
            % of knickzones in the tributary
            trib_id_matrix(k) = j;
            % give each knickzone the current tributary number
            beta_matrix(k) = current_beta; 
            % give each knickzone the current reference profile steepness
        end
        trib_array{j} = transpose(trib_id_matrix);  
        % these are the trib id's for the current tributary (repeat number
        % if multiple knickzones in the same tributary)
        beta_array{j} = transpose(beta_matrix);
        % these are the reference profile steepness's for the current
        % tributary
    end
    
    all_trib_ids_each_basin = [];
    all_beta_each_basin = [];
    % vertically concatinate for each tributary in the basin
    for j= 1:length(current_basin_elevation);
        all_trib_ids_each_basin = vertcat(all_trib_ids_each_basin,trib_array{j}); 
        % concatinate each trib id within the basin
        all_beta_each_basin = vertcat(all_beta_each_basin,beta_array{j}); 
        %same with slope of the reference profile
    end
    all_trib_ids_all_basins{i} = all_trib_ids_each_basin;  
    % store each trib id list for each basin in a different cell
    all_beta_all_basins{i}= all_beta_each_basin;


    % get the stream id number for each knickpoint too (which stream
    % did this come from), also assign measured Ksn from SA plot and from
    % chiplot
    
    % attributes for each basin
    stream_id_matrix= [];
    ksn_SA_matrix = [];
    ks_chiplot_matrix = []; 
    % ^^ note if fixed m/n is used, then this is also a measure of ksn
    m_n_basin_matrix = [];
    

    for k =1:length(all_trib_ids_all_basins{i}); 
        % for the number of knickpoints found in that stream
        stream_id_matrix(k) = i ;  % we want to record the current stream too
        
        ksn_SA_matrix(k) = ksn_SA_store(i);
        % save slope area derived ksn for basin (always uses theta_ref  <-- from parameters script)
        
        ks_chiplot_matrix(k) = ks_chiplot(i);
        % save chiplot derived ks for basin
        
         m_n_basin_matrix(k) = theta_bf(i);
         % save the m/n ratio used to calculate ks in chispace
        
    end
    stream_id_array{i} = transpose(stream_id_matrix);  
    % save the stream id matrix
    ksn_SA_array{i} = transpose(ksn_SA_matrix); 
    ks_chiplot_array{i} = transpose(ks_chiplot_matrix); 
    m_n_basin_array{i} = transpose(m_n_basin_matrix); 
    
end

all_trib_ids = [];
stream_id = [];
all_ks_trib_list = [];
ksn_SA_list = [];
ks_chiplot_list = [];
m_n_basinwide_list = [];

% vertically concatinate for each basin (to get full list)
for i = 1:num_basins2 % % for all stream basins
    all_trib_ids = vertcat(all_trib_ids,all_trib_ids_all_basins{i});  
% vertcally concatinate all the trib ids from each basin stored in the cell array above
    stream_id = vertcat(stream_id,stream_id_array{i});
    ksn_SA_list = vertcat(ksn_SA_list,ksn_SA_array{i});
    ks_chiplot_list = vertcat(ks_chiplot_list,ks_chiplot_array{i});
    m_n_basinwide_list = vertcat(m_n_basinwide_list,m_n_basin_array{i});
    all_ks_trib_list = vertcat(all_ks_trib_list,all_beta_all_basins{i});
end