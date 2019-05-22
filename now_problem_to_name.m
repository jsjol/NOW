function fn_l = now_problem_to_name(p)
% function fn_l = now_problem_to_name(p)
% Get a set of reasonable names for the waveform files based on the problem
% definition. This function may be replaced with header information stored
% in the file one such headers are supported by the sequence code.

gampstr = sprintf('%0.0f', p.gMax);
slewstr = sprintf('%0.0f', p.sMax);
normstr = num2str(p.useMaxNorm);
nosstr  = num2str(p.N);
etastr  = sprintf('%0.2f', p.eta);
mxwlstr = num2str(p.doMaxwellComp);

TT = eig(p.targetTensor);
TT = TT / max(TT);

tensstr = [num2str(TT(1),'%0.2f') '_' num2str(TT(2),'%0.2f') '_' num2str(TT(3),'%0.2f')];

fn = ['NOW_gMax-' gampstr '_sMax-' slewstr '_MaxNorm-' normstr '_DoMxwl-' mxwlstr '_N-' nosstr '_eta-' etastr '_T-' tensstr];

if p.durationZeroGradientRequested > 0
    
    t1str   = sprintf('%0.2f', p.durationFirstPartActual);
    t2str   = sprintf('%0.2f', p.durationSecondPartActual);
    tpstr   = sprintf('%0.2f', p.durationZeroGradientActual);
    
    fn = [fn '_dur-' t1str '_' tpstr '_' t2str];
    
    fn_l{1} = [fn '_AB'];
    fn_l{2} = [fn '_A'];
    fn_l{3} = [fn '_B'];
    
else
    ttotstr = num2str(p.totalTimeActual);
    
    fn = [fn '_dur-' ttotstr];
    
    fn_l{1} = [fn '_AB'];
end

