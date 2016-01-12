function [M,map12,map21]=confusion_matrix(T1,L1,T2,L2,opts)
% CONFUSION_MATRIX - confusion matrix between two event lists
%   
% [M,map12,map21]=confusion_matrix(T1,L1,T2,L2,opts)
%   Returns confusion matrix M, and the optimal mappings m12 and m21
%   between the mappings. 

% Inputs:
%   T1,L1 - list of times (real) and labels (in 1..K1) for firing event list 1
%   T2,L2 - list of times (real) and labels (in 1..K2) for firing event list 2
%   opts  - (optional) controls various parameters such as:
%           opts.max_matching_offset : max time difference to be considered a
%                                      match (positive real, default 3.0)
%           opts.map12: give priority to matches that coincide with this
%                       mapping
% Outputs:
%  M - (K1+1)-by-(K2+1) confusion matrix where K1=max(L1) and K2=max(L2).
%  map12 - integer vector of size K1 the optimal mapping of labels from
%          T1/L1 to labels from T2/L2
%  map21 - integer vector of size K2 the optimal mapping of labels from
%          T2/L2 to labels from T1/L1
%
% Run without arguments gives a self-test

% run the test if no arguments
if nargin==0, test_confusion_matrix; return; end
if (nargin<5) opts=struct; end;
if (~isfield(opts,'max_matching_offset')) opts.max_matching_offset=3; end;

K1=max(L1); K2=max(L2);
CM=zeros(K1+1,K2+1);

if (~isfield(opts,'map12'))
    opts.map12=zeros(1,K1);
    [~,opts.map12]=confusion_matrix(T1,L1,T2,L2,opts);
end;

for pass=1:2
    %on the first pass we are giving priority to matches that agree with opts.map12
    for offset=0:opts.max_matching_offset
        inds1_to_remove=zeros(1,length(T1)); %once paired up we remove those events
        inds2_to_remove=zeros(1,length(T2)); %once paired up we remove those events

        if (length(T2)>0)
            ptr2 = 1;  % moving index to sorted T2 list
            for ii1=1:length(T1) %this loop makes things slow (should we fix?)
                %find the indices of T2/L2 which should be paired with the event at T1(ii1)
                [ii2,ptr2] = indexlist(T2,T1(ii1),offset,ptr2);
                %only use those that haven't been removed!
                ii2=ii2(find(inds2_to_remove(ii2)==0));
                if (pass==1)
                    %only use those that agree with the mapping
                    ii2=ii2(find(opts.map12(L1(ii1))==L2(ii2)));
                end;
                if (length(ii2)>0)
                    %let's only use the first
                    ii2=ii2(1);
                    CM(L1(ii1),L2(ii2))=CM(L1(ii1),L2(ii2))+1; %increment the entry in the confusion matrix
                    inds1_to_remove(ii1)=1; %we've handled the event, so let's remove it!
                    inds2_to_remove(ii2)=1; %we've handled the event, so let's remove it!
                end;
            end;
        end;

        %Now let's remove the events that were marked above
        T1=T1(inds1_to_remove==0); 
        L1=L1(inds1_to_remove==0);
        T2=T2(inds2_to_remove==0);
        L2=L2(inds2_to_remove==0);
    end;
end;

%The rest are unclassified
for j=1:K1
    CM(j,K2+1)=length(find(L1==j));
end;
for j=1:K2
    CM(K1+1,j)=length(find(L2==j));
end;

M=CM;

HH=Hungarian(-M(1:K1,1:K2));

map12=zeros(1,K1);
for j=1:K1
    ind0=find(HH(j,1:K2)==1);
    if (length(ind0)==0) map12(j)=0;
    else map12(j)=ind0(1);
    end;
end;
map21=zeros(1,K2);
for j=1:K2
    ind0=find(HH(1:K1,j)==1);
    if (length(ind0)==0) map21(j)=0;
    else map21(j)=ind0(1);
    end;
end;

end

%%%%%%%
function [ii2,ptr2] = indexlist(t2,t1,off,ptr2) % find indices of t2 for which
% abs(t2-t1)<off, where t1 is a scalar. Update ptr2 which gives the rough
% center value of this index. 
N = numel(t2);
i = ptr2; while(t2(i)>=t1-off & i>1), i=i-1; end  % go down until t2's too early
j = ptr2; while(t2(j)<=t1+off & j<N), j=j+1; end  % go up until t2's too late
ii2 = i:j;
ii2=ii2(find(abs(t2(ii2)-t1)<=off));
if (length(ii2)>0) ptr2=ii2(1); end; %update the pointer
end


function test_confusion_matrix

T1=[10,20,   30,40,50,60,80   ];
L1=[1, 1,    1, 2, 2, 3, 4    ];
T2=[11,22,21,30,32,52,60      ];
L2=[2, 2, 1, 4, 2  1, 3       ];
opts=struct;
%opts.map12=[1,2,0];
[M,map12,map21]=confusion_matrix(T1,L1,T2,L2,opts);
M

end
