function [subg, count, GY] = gspan (G, minsup, size_req, ...
	boostY, boostWeights, boostTau, boostN, boostMax, boostType);

% gSpan frequent graph substructure algorithm.
%   and
% gSpan based graph boosting.
%
% (This is a Matlab wrapper around the modified gspan implementation of Taku
%  Kudo.)
%
% Author: Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
% Date: 15th November 2006
%
% Find all subgraphs in the set of graphs G that appear at least minsup times.
%
% Note: For graph boosting, this function should normally not be used
%    directly, but through the gboost1d5 or gboost2 functions.
%
%%%
% Input
%   G: (1,n) cellarray of n graph structures with this layout
%     g.nodelabels: (n,1) discrete integer labels [L_1 ; L_2 ; ... ; L_n];
%     g.edges: (m,2) edges, [from to] at each line:
%       [e_1_{from} e_1_{to} edgelabel_1 ; ... ; e_m_{from} e_m_{to} edgelabel_m]
%       The node indices go from 1 to n.  (They will be converted to 0-(n-1)
%       internally and converted back.)
%
%   minsup: The minimum number of times a frequent pattern has to appear in G.
%
%   (optional) size_req: Either [] for no limit, or (1,2) real col-vector
%      [min_node_count max_node_count], which limits the returned graphs to
%      have a number of nodes p, such that
%         min_node_count <= p <= max_node_count.
%      If max_node_count < min_node_count, only the lower bound is used.
%      (default: [0 0])
%
% For graph boosting:
%   boostY: (n,1) array of labels in {-1, 1}.
%
%   boostWeights: (n,1) array of reals in [0 ; 1].
%
%   boostTau: tau parameter (normally set to 0).  This is the start gain it
%      has to beat.
%
%   boostN: integer > 0.  From the best patterns, only the best boostN patterns
%      are returned.
%
%   boostMax: The maximum number of exploration steps to take or zero for
%      infinity (1e6 is a good choice here).
%
%   boostType: 1 for 1.5-class weak learners with outputs {0, 1}, 2 for
%      2-class weak learners with outputs {-1, 1}.
%
%%%
% Output
%   subg: (1,p) cellarray of p graph structures that appear frequently in G.
%
%   (optional) count: (1,p) uint32 array giving at element i the number of times
%     graphs{i} appears in G (duplicates within one graph only counted once).
%
%   (in case of graph boosting count is)
%   count: (1,p) cellarray of { -1, 1 }, giving the classifiers h_{<subg>,<subgY>}.
%
%   (optional) GY: (n,p) double array giving individual counts for each
%       subgraph/graph pair.  At (i,j) the number of counts of the j'th
%       subgraph found in the i'th graph is recorded.  This takes no
%       additional computation time to record.
%

if nargin < 2 || nargin > 10 || nargout > 3
	error(['Invalid number of input or output parameters.']);
end

d_size_req = [0 0];
if nargin >= 3
	if isempty(size_req)
		size_req = [0 0];
	end

	if size(size_req,1) ~= 1 || size(size_req,2) ~= 2
		error(['size_req parameter invalid.']);
	end
	d_size_req = size_req;
end

for i=1:length(G)
	[v, r] = verifygraph (G{i});
	if v ~= 1
		error(['Graph ', num2str(i), ': ', r]);
	end
end
if size(G,2) == 1
	G=G';	% put it as row.
end

if nargin > 3
	disp(['Starting gspan-boost run...']);
	if boostTau < 0.0
		boostTau = 0.0;
	end

	if nargout == 3
		[subg, count, GY] = mexgspan (G, minsup, d_size_req, 0, ...
			boostY, boostWeights, boostTau, boostN, boostMax, boostType);
	elseif nargout == 2
		[subg, count] = mexgspan (G, minsup, d_size_req, 0, ...
			boostY, boostWeights, boostTau, boostN, boostMax, boostType);
	elseif nargout == 1
		[subg] = mexgspan (G, minsup, d_size_req, 0, ...
			boostY, boostWeights, boostTau, boostN, boostMax, boostType);
	end
else
	disp(['Starting normal gspan run...']);
	if nargout == 3
		[subg, count, GY] = mexgspan (G, minsup, d_size_req, 0);
	elseif nargout == 2
		[subg, count] = mexgspan (G, minsup, d_size_req, 0);
	elseif nargout == 1
		[subg] = mexgspan (G, minsup, d_size_req, 0);
	end
end

% Filter out duplicates (this is due to undirected edges)
for i=1:length(subg)
	subg{i}.edges = unique (subg{i}.edges, 'rows');
end


function [valid, reason] = verifygraph (G);

% Verify if the directed graph G is valid, that is:
%
% 1. Does not contain any self-referential or duplicate edges.
% 2. Uses proper indices.

nodes = size(G.nodelabels,1);

% Check for indices
if size(G.edges, 1) > 0
	I = find (G.edges(:,1) <= 0 | G.edges(:,1) > nodes | ...
		G.edges(:,2) <= 0 | G.edges(:,2) > nodes);
	if length(I) ~= 0
		valid = 0;
		if nargout >= 2
			reason = ['Edge indices not within 1-', num2str(nodes), '.'];
		end
		return;
	end

	% Check for self referential edges (self-loops)
	I = find (G.edges(:,1) == G.edges(:,2));
	if length(I) ~= 0
		valid = 0;
		if nargout >= 2
			reason = ['Graph has ', num2str(length(I)), ' self-loops.'];
		end
		return;
	end
end

% Check for duplicate edges
% DISABLED as both our gspan modification as well as mexrelabel support
%    multiple edge attributes
% TODO
if false
	E = G.edges(:,1:2);
	Eu = unique(E,'rows');
	if size(E,1) ~= size(Eu,1)
		valid = 0;
		if nargout >= 2
			reason = ['There are ', num2str(size(E,1)-size(Eu,1)), ...
				' duplicate edges.'];
		end
		return;
	end
end

if nargout >= 2
	reason = 'Graph is ok.';
end
valid = 1;

