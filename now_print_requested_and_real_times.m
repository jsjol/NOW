function now_print_requested_and_real_times(o)

fprintf(1, '------------ Requested timing parameters: ------------ \n');
fprintf(1, 'DurPre = %5.3f  DurPost = %5.3f  DurPi = %5.3f  [ms]\n\n', o.durationFirstPartRequested, o.durationSecondPartRequested, o.durationZeroGradientRequested);
fprintf(1, '------------   Actual timing parameters:  ------------ \n');
fprintf(1, 'DurPre = %5.3f  DurPost = %5.3f  DurPi = %5.3f  [ms]\n\n', o.durationFirstPartActual, o.durationSecondPartActual, o.durationZeroGradientActual);
