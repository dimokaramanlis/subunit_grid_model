function  plotSubunitGrid(subcents, subwts)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

subcents = double(gather(subcents));
subwts   = double(gather(subwts));


opts = {'QJ', 'Pp'};
[v,c] = voronoin(subcents, opts); 
gridXc = mean(subcents(:,1));
gridYc = mean(subcents(:,2));

% maxd = max(range(subcents(:,1)), range(subcents(:,2)));
% 
% badVerts = any(abs(v - [gridXc, gridYc]) > maxd/2, 2);

badVerts = find(sqrt(sum((v - [gridXc, gridYc]).^2, 2)) > 5e5);
badFaces = cellfun(@(x) any(ismember(x, badVerts)), c);

f = CelltoMatUE(c(~badFaces));

patch('Faces', f, 'Vertices', v, 'FaceVertexCData', subwts(~badFaces),...
    'FaceColor','flat', 'EdgeAlpha', 0.1);


end

