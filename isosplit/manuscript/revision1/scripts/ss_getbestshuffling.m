function [j ashuf] = getbestshuffling(a,b,opts)
% Snapshot of getbestshuffling.m on 5/22/2015
% GETBESTSHUFFLING - return best permutation between two lists of labels
%
% j = getbestshuffling(a,b) where a and b are row or column vectors of
% integer labels in 1 to K, tries to returns a permutation j taking a's
% labels to b's with fewest errors. That is, b(i) = j(a(i)) holds for the
% most possible choices of i from 1...N, where N is length of a and b.
%
% [j ashuf] = getbestshuffling(a,b,opts) also returns ashuf = j(a) the relabeled
%  first input.
%
%  opts is optional
%  opts.method = 
%     'h' for hungarian (default) -- depends on Hungarian.m
%     'g' for greedy
%
% Notes: crude greedy alg.
% Jeremy mentions Hungarian alg for this.
%
% Barnett 12/9/14. ashuf 1/9/15
% jfm added Hungarian as default option 1/9/15
% jfm fixed a problem in the case of K1>K2 1/10/15

addpath([fileparts(mfilename('fullpath')),'/internal']);

if nargin<1, test_getbestshuffling; return; end
if numel(a)~=numel(b), error('a and b must have same length'); end

if (nargin<3) opts.method='h'; end; %default to hungarian algorithm

K1=max(a);
K2=max(b);
if (K1>K2) %This case will cause problems because 2 labels will get mapped into the same label
	
	% This part is trickier than one would expect
	% If K1>K2, then 2 labels of a will get mapped to the same label of b
	% This is not desirable since we want a true permutation of labels of a
	% So the trick is to map the labels of b to a (which will be
	% one-to-one). Then we invert the map. But of course there are some
	% labels of a that aren't in the range of the inversion. What to do with those labels? At
	% first I thought, let's just map them to themselves. However, that
	% doesn't always work. Instead we need to map them to unique labels
	% that haven't yet been used.
	% 
	% If the reader has a more elegant way of handling this, please let me
	% know. -jfm
	
	
	[jinv,bshuf]=ss_getbestshuffling(b,a,opts);
	% jinv is K2x1
	next_one_not_used=K2+1;
	j=zeros(1,K1);
	for k1=1:K1
		tmp0=find(jinv==k1);
		if (length(tmp0)>0)
			j(k1)=tmp0(1);
		else
			j(k1)=next_one_not_used; %here's the tricky part!
			next_one_not_used=next_one_not_used+1;
		end;
	end;
	ashuf=zeros(size(a));
	ashuf(find(a>0))=j(a(find(a>0)));
	return;
end;

if (strcmp(opts.method,'h'))
	[j,ashuf]=getbestshuffling_hungarian(a,b,opts);
elseif (strcmp(opts.method,'g'))
	[j,ashuf]=getbestshuffling_greedy(a,b,opts);
else
	error('Unrecognized shuffling method');
end;

end

function [j,ashuf]=getbestshuffling_hungarian(a,b,opts)

K1=max(a);
K2=max(b);

if (K1>K2) error('Unexpected problem'); end;

AA=zeros(K1,K2);
for i=1:K1
  for j=1:K2
    AA(i,j)=sum((a==i).*(b==j));
  end;
end;

MM=Hungarian(-AA);
inds=zeros(1,K1);
for k1=1:K1 
	tmp=find(MM(k1,:)==1);
	if (length(tmp)>0)
		inds(k1)=tmp(1);
	else
		inds(k1)=k1;
	end;
end;

j=inds;

ashuf=zeros(size(a));
ashuf(find(a>0))=j(a(find(a>0)));

end

function [j,ashuf]=getbestshuffling_greedy(a,b,opts)

K1=max(a);
K2=max(b);

if (K1>K2) error('Unexpected problem'); end;

sizea = size(a);
a = a(:); b = b(:);
K = max([a;b]);
j = nan(1,K); remain = 1:K;  % remain = which labels in b available each step
for i=1:K               % loop over labels in a
  B = bsxfun(@eq,b,remain); % N-by-(K-i+1) logical array; note @eq not @isequal
  [~,k] = max(sum(bsxfun(@eq, a==i, B),1)); % k is best of remaining in b
  j(i) = remain(k);
  remain = remain([1:k-1, k+1:(K-i+1)]);  % kill off this option in b
end
% shuffle classified labels only...
ashuf = a; clas = find(a>0); ashuf(clas) = j(a(clas));
ashuf = reshape(ashuf,sizea);

end

function test_getbestshuffling  % todo: test ashuf output
disp('testing K1=K2...')
N = 1e3;  % number of cases
K = 5;    % how many label types
a = randi(K,1,N);
je = randperm(K);   % the true shuffling
b = je(a);
j = getbestshuffling(a,b);
fprintf('noise-free case: errors in the shuffling = %d\n',numel(find(j~=je)))

i = rand(1,N)<0.05; % some indices to mess up in b
b(i) = randi(K,1,numel(find(i))); % mess them up
j = getbestshuffling(a,b);
fprintf('noisy case: errors in the shuffling = %d\n',numel(find(j~=je)))

opts.method='h'; disp('testing K1>K2...')
[pp ashuf]=getbestshuffling([1,1,1,2,2,2,3,3,3,4,4,4,5,5,5],[1,1,1,2,2,3,3,3,3,4,4,4,2,2,2])

end