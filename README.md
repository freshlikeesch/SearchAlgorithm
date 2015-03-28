# SearchAlgorithm
Mathematical PageRanks for a simple network, expressed as percentages. (Google uses a logarithmic scale.) Page C has a higher PageRank than Page E, even though there are fewer links to C; the one link to C comes from an important page and hence is of high value. If web surfers who start on a random page have an 85% likelihood of choosing a random link from the page they are currently visiting, and a 15% likelihood of jumping to a page chosen at random from the entire web, they will reach Page E 8.1% of the time. (The 15% likelihood of jumping to an arbitrary page corresponds to a damping factor of 85%.) Without damping, all web surfers would eventually end up on Pages A, B, or C, and all other pages would have PageRank zero. In the presence of damping, Page A effectively links to all pages in the web, even though it has no outgoing links of its own.

Expressed in code, this looks like:
% Parameter M adjacency matrix where M_i,j represents the link from 'j' to 'i', such that for all 'j' sum(i, M_i,j) = 1
% Parameter d damping factor
% Parameter v_quadratic_error quadratic error for v
% Return v, a vector of ranks such that v_i is the i-th rank from [0, 1]
 
function [v] = rank(M, d, v_quadratic_error)
 
N = size(M, 2); % N is equal to half the size of M
v = rand(N, 1);
v = v ./ norm(v, 2);
last_v = ones(N, 1) * inf;
M_hat = (d .* M) + (((1 - d) / N) .* ones(N, N));
 
while(norm(v - last_v, 2) > v_quadratic_error)
	last_v = v;
	v = M_hat * v;
	v = v ./ norm(v, 2);
end
 
endfunction
 
function [v] = rank2(M, d, v_quadratic_error)
 
N = size(M, 2); % N is equal to half the size of M
v = rand(N, 1);
v = v ./ norm(v, 1);   % This is now L1, not L2
last_v = ones(N, 1) * inf;
M_hat = (d .* M) + (((1 - d) / N) .* ones(N, N));
 
while(norm(v - last_v, 2) > v_quadratic_error)
	last_v = v;
	v = M_hat * v;  
        % removed the L2 norm of the iterated PR
end
 
endfunction


Example of code calling the rank function defined above:
M = [0 0 0 0 1 ; 0.5 0 0 0 0 ; 0.5 0 0 0 0 ; 0 1 0.5 0 0 ; 0 0 0.5 1 0];
rank(M, 0.80, 0.001)


This example takes 13 iterations to converge.
The following is a proof that rank.m is incorrect. It is based on the first graphic example. My understanding is that rank.m uses the wrong norm on the input, then continues to renormalize L2, which is unnecessary.
% This represents the example graph, correctly normalized and accounting for sinks (Node A) 
% by allowing it to effectively random transition 100% of time, including to itself.
% While RANK.m doesn't actually handle this incorrectly, it does not show exactly how one should
% handle sink nodes (one possible solution would be a SELF-TRANSITION of 1.0), which does not
% give the correct result.
 
test_graph = ...
[  0.09091   0.00000   0.00000   0.50000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000;
   0.09091   0.00000   1.00000   0.50000   0.33333   0.50000   0.50000   0.50000   0.50000   0.00000   0.00000;
   0.09091   1.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000;
   0.09091   0.00000   0.00000   0.00000   0.33333   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000;
   0.09091   0.00000   0.00000   0.00000   0.00000   0.50000   0.50000   0.50000   0.50000   1.00000   1.00000;
   0.09091   0.00000   0.00000   0.00000   0.33333   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000;
   0.09091   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000;
   0.09091   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000;
   0.09091   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000;
   0.09091   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000;
   0.09091   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000 ]
 
pr = rank(test_graph, 0.85, 0.001)   % INCORRECT is not normalized.
 
%  0.062247
%  0.730223
%  0.650829
%  0.074220
%  0.153590
%  0.074220
%  0.030703
%  0.030703
% 0.030703
%  0.030703
%  0.030703
 
pr / norm(pr,1)    % CORRECT once normalized. I still don't know why the L2 normalization happens ( v = v/norm(v, 2))
 
%   0.032781
%   0.384561
%   0.342750
%   0.039087
%   0.080886
%   0.039087
%   0.016170
%   0.016170
%   0.016170
%   0.016170
%   0.016170
 
pr = rank2(test_graph, 0.85, 0.001) % CORRECT, only requires input PR normalization (make sure it sums to 1.0)
 
%   0.032781
%   0.384561
%   0.342750
%   0.039087
%   0.080886
%   0.039087
%   0.016170
%   0.016170
%   0.016170
%   0.016170
%   0.016170


