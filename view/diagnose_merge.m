function diagnose_merge(info)
% DIAGNOSE_MERGE.  plot figures from merge_across_channels's diagnostic output
%
% Input: info is the struct returned by ms_merge_across_channels
% Barnett 4/22/16

Kpm = numel(info.premergelabelchans);  % # label types per-merge
info.premergelabelchans
S = nan(Kpm); r = S; pr = S; cf = S; sdc = S; % extract merge info struct...
for j=1:Kpm, for k=1:Kpm, s = info.Sinfo{j,k};
    if isfield(s,'s'), S(j,k) = s.s; end
    if isfield(s,'r12'), r(j,k) = s.r12; end
    if isfield(s,'peakratio'), pr(j,k) = s.peakratio; end
    if isfield(s,'coincfrac'), cf(j,k) = s.coincfrac; end
    if isfield(s,'stddevC'), sdc(j,k) = s.stddevC; end
  end, end
figure; imagesc(S); colorbar; title('to merge Boolean');
figure;
subplot(2,2,1);imagesc(r);title('r (corr coeff)');colorbar
subplot(2,2,2); imagesc(pr); title('peakratio');colorbar
subplot(2,2,3);imagesc(cf);title('coincfrac');colorbar
subplot(2,2,4); imagesc(sdc); title('stddevC');colorbar
