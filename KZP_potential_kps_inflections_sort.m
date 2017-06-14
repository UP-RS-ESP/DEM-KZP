function [base_cells_current_trib] = ...
    KZP_potential_kps_inflections_sort(kps_cells_current_trib,kps_lips_length,...
    kps_base_length, base_cells_current_trib,stream_mouth_cell)
% Reorganize knickzone base vector to make sure upstream lip is compared
% to downstream base!

% CASE 1
% downstream        <------              Upstream
%  knickpoint - base - knickpoint - base
% not what we want.. make adjustments
if ~isempty(kps_cells_current_trib)
    if kps_lips_length == kps_base_length && min(kps_cells_current_trib)> min(base_cells_current_trib) ;
        % if the stream starts out with the first "0" detrend gradient as a
        % base (near headwaters)^^ and then ends with the
        % last "0" detrend gradient as a knickpoint (near mouth
        % of stream)

        base_cells_current_trib = base_cells_current_trib(2:kps_base_length,1);
        % ^^ then we need to remove the first base and add the mouth of the stream so all the
        % knickpoint lips will be compared to their respective
        % knickpoint bases, and the last knickpoint near the
        % mouth of the stream will be compared to the mouth of the stream

        base_cells_current_trib = [base_cells_current_trib;stream_mouth_cell];
        % Lastly, we must add the mouth of the stream, so the last
        % knickpoint will be compared to the mouth of the stream "n"
    end
end

% End with:
% downstream        <------              Upstream
%  mouth stream - knickpoint - base - knickpoint
% ^^ that will work

%
% CASE 2
% downstream        <------              Upstream
%  knickpoint - base - knickpoint - base - knickpoint
% this needs to be fixed too, the knickpoint at the end of the stream needs
% a reference base

if kps_lips_length> kps_base_length ;
    base_cells_current_trib = [base_cells_current_trib;stream_mouth_cell];
end
% if there are more knickpoints than bases, this loop
% adds the last cell
% point (mouth of the stream 'n') to the end of the bases vector to make them
% the same size. So we compare the elevation of the last knickpoint
% (furthest downstream) to the elevation of the mouth of the stream.

% End with:
% downstream        <------              Upstream
%  mouth stream - knickpoint - base - knickpoint
%^^ that will work

%
% CASE 3
% downstream        <------              Upstream
%  base - knickpoint - base - knickpoint - base
% ^ this won't work because the first inflection is a base

if kps_lips_length < kps_base_length;
    base_cells_current_trib = base_cells_current_trib(2:kps_base_length(1),1);
end

% if there are more base than knickpoints, this loop removes the first
% base (the extra base that occurs near the headwaters of the stream)

% End with:
% downstream        <------              Upstream
%  base - knickpoint - base - knickpoint
% ^ that will work

% Case 4
% downstream        <------              Upstream
%                     base
% ^ here there is only 1 base and no knickpoints. Just remove the
% base because we aren't interested in that

if kps_base_length == 1 && kps_lips_length == 0
    base_cells_current_trib = [];
end
%
% ALWAYS WANT EITHER:
% downstream        <------              Upstream
%  base - knickpoint - base - knickpoint
% or
% downstream        <------              Upstream
%  mouth stream - knickpoint - base - knickpoint

%  All of the knickpoint cells are organized so the upstream knickpoint 
% lips are compared to a respective downstream knickpoint base
