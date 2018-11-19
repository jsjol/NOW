function h = now_plot_all(r)
% function h = now_plot_all(r)
% Plot an overview and analysis of the gradient waveform from the NOW
% result structure. Note that this function serves as a template for integration with 
% the multidimensional diffusion framework, available here: 
% https://github.com/markus-nilsson/md-dmri

h = gwf_plot_all(r.gwf, r.rf, r.dt);




