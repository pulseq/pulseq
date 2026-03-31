function [grad, times, amplitudes] = makeExtendedTrapezoidArea(channel, grad_start, grad_end, area, sys)
% Make the shortest possible extended trapezoid for a given area and edge values
% This version is the one with the fixed flat top and was derived from the
% corresponding PyPulseq version by by Mehmet Emin Öztürk. This implementation
% is both faster and more accurate than the previous one and it runs in Octave.
% Main methodology in explained in Python version
% Some variable names might also be different

    if nargin < 5 || isempty(sys)
        sys = default_opts(); % Define your default system
    end

    max_slew = sys.maxSlew * 0.99;
    max_grad = sys.maxGrad * 0.99;
    raster_time = sys.gradRasterTime;

    min_duration = max(round(calc_ramp_time(grad_end, grad_start, max_slew, raster_time) / raster_time), 2);

    % Estimate upper bound duration
    max_duration = max( ...
        [round(calc_ramp_time(0, grad_start, max_slew, raster_time) / raster_time), ...
        round(calc_ramp_time(0, grad_end, max_slew, raster_time) / raster_time), ...
        min_duration]);

    % Try to find a solution linearly
    solution = [];
    for duration = min_duration:max_duration
        solution = find_solution(duration, area, grad_start, grad_end, max_slew, max_grad, raster_time);
        if ~isempty(solution)
            break;
        end
    end

    % Binary search if linear search failed
    if isempty(solution)
        duration = max_duration;
        while isempty(solution)
            duration = duration * 2;
            solution = find_solution(duration, area, grad_start, grad_end, max_slew, max_grad, raster_time);
        end

        solution = binary_search(@(d) find_solution(d, area, grad_start, grad_end, max_slew, max_grad, raster_time), ...
                                 floor(duration/2), duration);
    end

    time_ramp_up = solution(1) * raster_time;
    flat_time    = solution(2) * raster_time;
    time_ramp_down = solution(3) * raster_time;
    grad_amp     = solution(4);

    % Generate final time vector and amplitudes
    if flat_time > 0
        times = cumsum([0, time_ramp_up, flat_time, time_ramp_down]);
        amplitudes = [grad_start, grad_amp, grad_amp, grad_end];
    else
        times = cumsum([0, time_ramp_up, time_ramp_down]);
        amplitudes = [grad_start, grad_amp, grad_end];
    end

    grad=mr.makeExtendedTrapezoid(channel,'system',sys,'times',times, 'amplitudes', amplitudes);
    if abs(grad.area - area) >= 1e-3
        error('Could not find a solution for area=%.6f.', area);
    end
end


function time = to_raster(time_val, raster_time)
    time = ceil(time_val / raster_time) * raster_time;
end

function t = calc_ramp_time(g1, g2, max_slew, raster_time)
    t = to_raster(abs(g1 - g2) / max_slew, raster_time);
end

function sol = binary_search(fun, low, high)
    while low < high - 1
        mid = floor((low + high) / 2);
        if ~isempty(fun(mid))
            high = mid;
        else
            low = mid;
        end
    end
    sol = fun(high);
end


function sol = find_solution(duration, area, grad_start, grad_end, max_slew, max_grad, raster_time)
    sign_area = sign(area);
    grad_amp = sign_area * max_grad;

    % Convert to raster steps
    ru_min = abs(grad_amp - grad_start) / max_slew / raster_time;
    rd_min = abs(grad_amp - grad_end) / max_slew / raster_time;
    flat_time = max(duration - ru_min - rd_min, 0);

    % Check if feasible
    approx_area = ru_min * (grad_amp + grad_start) + ...
                  rd_min * (grad_amp + grad_end) + ...
                  2 * flat_time * grad_amp;

    if abs(2 * area / raster_time) > abs(approx_area)
        sol = [];
        return;
    end

    % Easy early solution: max_grad
    ru = (duration * max_slew * raster_time + sign_area * (grad_end - grad_start)) / (2 * max_slew * raster_time);
    if sign_area * grad_start + ru * max_slew * raster_time > max_grad + 1e-5
        ru_steps = round(abs(grad_start - sign_area * max_grad) / max_slew / raster_time);
        rd_steps = round(abs(grad_end - sign_area * max_grad) / max_slew / raster_time);
        flat_steps = duration - ru_steps - rd_steps;
        if flat_steps > 0
            grad_amp = -(ru_steps * raster_time * grad_start + ...
                         rd_steps * raster_time * grad_end - 2 * area) / ...
                        ((ru_steps + 2 * flat_steps + rd_steps) * raster_time);
            amps = [grad_start, grad_amp, grad_amp, grad_end];
            t = cumsum([0, ru_steps, flat_steps, rd_steps]) * raster_time;
            slew = diff(amps) ./ diff(t);
            if max(abs(slew)) < max_slew + 1e-5 && max(abs(amps)) < max_grad
                sol = [ru_steps, flat_steps, rd_steps, grad_amp];
                return;
            end
        end
    end

    % Conservative downscaling
    while abs(2 * area / raster_time) < abs(approx_area)
        grad_amp = grad_amp / 2;
        if abs(grad_amp) < abs(max_grad) / 10
            ru_min = 0; rd_min = 0;
            flat_time = max(duration - ru_min - rd_min, 0);
            break;
        end
        ru_min = abs(grad_amp - grad_start) / max_slew / raster_time;
        rd_min = abs(grad_amp - grad_end) / max_slew / raster_time;
        flat_time = max(duration - ru_min - rd_min, 0);
        approx_area = ru_min * (grad_amp + grad_start) + ...
                      rd_min * (grad_amp + grad_end) + ...
                      2 * flat_time * grad_amp;
    end

    % Convert to integer steps
    ru_min = floor(ru_min);
    rd_min = floor(rd_min);
    ru_limit = ceil(abs(sign_area * max_grad - grad_start) / max_slew / raster_time) + 1;
    rd_limit = ceil(abs(sign_area * max_grad - grad_end) / max_slew / raster_time) + 1;

    % All combinations
    [RU, RD] = meshgrid(ru_min:ru_limit, rd_min:rd_limit);
    RU = RU(:);
    RD = RD(:);
    valid_mask = RD < (duration - RU);
    RU = RU(valid_mask);
    RD = RD(valid_mask);

    % Flat = 0 case
    num = (2 * area - duration * (grad_end + grad_start) * raster_time);
    denom_min = (grad_start - grad_end + sign_area * duration * max_slew * raster_time);
    denom_max = (-grad_start + grad_end + sign_area * duration * max_slew * raster_time);
    ru_flat0_min = round(num / denom_min / raster_time);
    ru_flat0_max = duration - round(num / denom_max / raster_time);

    RU = [RU; (ru_flat0_min:ru_flat0_max)'];
    RD = [RD; (duration - (ru_flat0_min:ru_flat0_max))'];

    % Filter invalid ones
    flat = duration - RU - RD;
    valid = flat >= 0 & RU > 0 & RD > 0;
    RU = RU(valid); RD = RD(valid); flat = flat(valid);

    % Calculate amp
    grad_amp = -(RU * raster_time * grad_start + RD * raster_time * grad_end - 2 * area) ./ ...
               ((RU + 2 * flat + RD) * raster_time);

    % Slew
    slew1 = abs(grad_start - grad_amp) ./ (RU * raster_time);
    slew2 = abs(grad_end - grad_amp) ./ (RD * raster_time);

    % Valid gradient/slew combinations
    valid = abs(grad_amp) <= max_grad + 1e-5 & slew1 <= max_slew + 1e-5 & slew2 <= max_slew + 1e-5;
    if ~any(valid)
        sol = [];
        return;
    end

    % Pick lowest slew
    ind = find(valid);
    [~, min_idx] = min(slew1(ind) + slew2(ind));
    best = ind(min_idx);

    sol = [RU(best), flat(best), RD(best), grad_amp(best)];
end
