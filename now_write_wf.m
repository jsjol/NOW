function fn = now_write_wf(r, p, o_dir, fn)
% function fn = now_write_wf(r, p, o_dir, fn)
% Write result from NOW_RUN to a .txt file.
% If the file names (fn) are not defined, a name will be created based on
% the problem definition.


if nargin < 3
    o_dir = pwd;
end

if nargin < 4
    fn = now_problem_to_name(p);
end


g = r.g / max(abs(r.g(:)));

mkdir(o_dir);

fn{1} = wf_save_wf_to_txt_file( g(:,1),  g(:,2),  g(:,3), fn{1}, o_dir);


if ~isempty(r.zind)
    end_a = r.zind(1) +1;
    beg_b = r.zind(end) +1;
    
    ga =  g(1:end_a, :);
    gb = -g(beg_b:end, :); % Minus sign to avoid rev gamp on scanner (ie assumes 180 pulse)
    
    fn{2} = wf_save_wf_to_txt_file( ga(:,1),  ga(:,2),  ga(:,3), fn{2}, o_dir);
    fn{3} = wf_save_wf_to_txt_file( gb(:,1),  gb(:,2),  gb(:,3), fn{3}, o_dir);
end


end



function out_fn = wf_save_wf_to_txt_file(gx, gy, gz, out_name, out_dir)

out_fn = [out_dir filesep out_name '.txt'];

f_h = fopen(out_fn, 'w');

formatspec = '%8.6f %8.6f %8.6f\r\n';

out_mat_1st = [gx(:) gy(:) gz(:)]';

fprintf(f_h, '%8.0i\r\n', size(out_mat_1st, 2));
fprintf(f_h, formatspec, out_mat_1st);


fclose(f_h);
end