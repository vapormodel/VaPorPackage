function terminations = findTerms(stroke_loc,vessels)
%FINDTERM Summary of this function goes here
%   Detailed explanation goes here
%to_search = [stroke_loc];
to_search = zeros(2e6,1);
to_search_marker = 1;
to_search_insert_point = 1;
to_search(to_search_insert_point) = stroke_loc;
to_search_insert_point = to_search_insert_point + 1;
% terminations = [];
terminations = zeros(2e6,1);
term_insert_point = 1;
% Find segments which have this as the parent.

%while ~isempty(to_search)
while to_search(to_search_marker) ~= 0
    this_search = to_search(to_search_marker);
    %to_search = to_search(1:end-1);
    children = find(vessels(:,7) == this_search);
    if isempty(children)
        % This node has no children and is a terminating node.
        %terminations = [terminations; this_search];
        terminations(term_insert_point) = this_search;
        term_insert_point = term_insert_point + 1;
    else
        % Push children to heap.
        %to_search = [to_search; children];
        to_search(to_search_insert_point:to_search_insert_point+length(children)-1) = children;
        to_search_insert_point = to_search_insert_point + length(children);
    end
    to_search_marker = to_search_marker + 1;
end

% Trim the terminations vector. 
terminations = terminations(logical(terminations));

end

