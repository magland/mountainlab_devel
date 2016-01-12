#include "confusion_matrix.h"

bool confusion_matrix(const char *clusters1_path, const char *clusters2_path, const char *output_path, int max_matching_offset)
{

	for (int offset=0; offset<=max_matching_offset; offset++) {
		QSet<int> inds1_to_remove;
	}

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
}
