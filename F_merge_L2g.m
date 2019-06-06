function outp = F_merge_L2g(collocated_cris_l2,fn)
% merging two L2g mat files together. written by Kang Sun on 2019/02/21 for
% AMoN intercomparison

outp = [];
if length(fn) == 1
    outp.merged.output_subset = collocated_cris_l2.(fn{1}).output_subset;
    return
else
    localfn = fieldnames(collocated_cris_l2.(fn{1}).output_subset);
    for ilocalfield = 1:length(localfn)
        outp.merged.output_subset.(localfn{ilocalfield}) = [];
        for ifield = 1:length(fn)
            outp.merged.output_subset.(localfn{ilocalfield}) = cat(1,...
                outp.merged.output_subset.(localfn{ilocalfield}),...
                collocated_cris_l2.(fn{ifield}).output_subset.(localfn{ilocalfield}));
        end
    end
end