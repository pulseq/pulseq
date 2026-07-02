function [grad, times, amplitudes] = makeHexagonGradientArea(channel, grad_start, grad_end, area, sys)
% Make the shortest possible hexagonal gradien (creating an extended
% trapezoid object) for a given area and edge values. In contrast to
% mr.makeExtendedTrapezoidArea(), this function creates a generic
% polynomial gradient object without a plato between the vortex2 and
% vortex3. We expect this function to find better solutions by relaxing the
% constraint of two vertices having the same amplitude, however, as seen in 
% testCase_12, this objective cannot be achieved in all cases by the current 
% implementation. 
% Parameters and methodology are explained in mr.makeExtendedTrapezoidArea(). 
% This function was implemented by Mehmet Emin Öztürk during his visit to
% Freiburg with some input from Maxim Zaitsev.

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
    times = [];
    amplitudes = [];
    for duration = min_duration:max_duration
        [times, amplitudes] = find_solution(duration, area, grad_start, grad_end, max_slew, max_grad, raster_time);
        if ~isempty(times)
            break;
        end
    end

    % Binary search if linear search failed
    if isempty(times)
        duration = max_duration;
        while isempty(times)
            duration = duration * 2;
            [times, ~] = find_solution(duration, area, grad_start, grad_end, max_slew, max_grad, raster_time);
        end

        [times, amplitudes] = binary_search(@(d) find_solution(d, area, grad_start, grad_end, max_slew, max_grad, raster_time), ...
                                 floor(duration/2), duration);
    end
    if any(diff(times) == 0)
        times = [times(1), times(2), times(4)];
        amplitudes = [amplitudes(1), amplitudes(2), amplitudes(4)];
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

function [times, amplitudes] = binary_search(fun, low, high)
    while low < high - 1
        mid = floor((low + high) / 2);
        if ~isempty(fun(mid))
            high = mid;
        else
            low = mid;
        end
    end
    [times, amplitudes] = fun(high);
end


function [times, amplitudes] = find_solution(duration, area, grad_start, grad_end, max_slew, max_grad, raster_time)
% Find extended trapezoid gradient waveform for given duration

    sign_area = sign(area);
    grad_amp = sign_area * max_grad;
    ramp_up_times = [];
    ramp_down_times = [];

    % Early estimation
    ru_min = abs(grad_amp - grad_start) / max_slew / raster_time;
    rd_min = abs(grad_amp - grad_end) / max_slew / raster_time;
    flat_time = max(duration - ru_min - rd_min, 0);

    area_check = ru_min * (grad_amp + grad_start) + rd_min * (grad_amp + grad_end) + 2 * flat_time * grad_amp;

    if abs(2 * area / raster_time) > abs(area_check)
        times = []; amplitudes = [];
        return;
    end

    % Fast-case: max_grad ramp with valid flat
    ru = (duration * max_slew * raster_time + sign_area * (grad_end - grad_start)) / (2 * max_slew * raster_time);
    if sign_area * grad_start + ru * max_slew * raster_time > max_grad + 1e-5
        ru_steps = round(abs(grad_start - sign_area * max_grad) / max_slew / raster_time);
        rd_steps = round(abs(grad_end - sign_area * max_grad) / max_slew / raster_time);
        flat_steps = duration - ru_steps - rd_steps;
        if flat_steps > 0
            grad_amp = -(ru_steps * raster_time * grad_start + rd_steps * raster_time * grad_end - 2 * area) / ...
                      ((ru_steps + 2 * flat_steps + rd_steps) * raster_time);
            amps = [grad_start, grad_amp, grad_amp, grad_end];
            t = cumsum([0, ru_steps, flat_steps, rd_steps]) * raster_time;
            slew = diff(amps) ./ diff(t);
            if max(abs(slew)) < max_slew + 1e-5 && max(abs(amps)) < max_grad
                times = t;
                amplitudes = amps;
                return;
            end
        end
    end

    % Gradually reduce grad_amp if area too large
    while abs(2 * area / raster_time) < abs(area_check)
        grad_amp = grad_amp / 2;
        if abs(grad_amp) < abs(max_grad) / 10
            ru_min = 0; rd_min = 0;
            flat_time = max(duration - ru_min - rd_min, 0);
            break;
        end
        ru_min = abs(grad_amp - grad_start) / max_slew / raster_time;
        rd_min = abs(grad_amp - grad_end) / max_slew / raster_time;
        flat_time = max(duration - ru_min - rd_min, 0);
        area_check = ru_min * (grad_amp + grad_start) + rd_min * (grad_amp + grad_end) + 2 * flat_time * grad_amp;
    end

    % Discrete timing limits
    ru_min = floor(ru_min);
    rd_min = floor(rd_min);
    ru_limit = ceil(abs(sign_area * max_grad - grad_start) / max_slew / raster_time);
    rd_limit = ceil(abs(sign_area * max_grad - grad_end) / max_slew / raster_time);

    flat_time = duration - min(rd_min, rd_limit) - min(ru_min, ru_limit);
    flat_time_min = duration - rd_limit - ru_limit;

    min_dif_area = area;
    i = -1;

    while flat_time > max(flat_time_min, -1)
        i = i + 1;
        ru_max = ru_min + i;
        rd_max = rd_min + i;
        flat_time = duration - min(rd_min + i, rd_limit) - min(ru_min + i, ru_limit);

        if flat_time <= 0
            % Fallback: flat_time = 0 or 1
            ru_0min = (2 * area - duration * (grad_end + grad_start) * raster_time) / ...
                      (grad_start - grad_end + sign_area * duration * max_slew * raster_time) / raster_time;
            ru_0max = duration - (2 * area - duration * (grad_end + grad_start) * raster_time) / ...
                      (-grad_start + grad_end + sign_area * duration * max_slew * raster_time) / raster_time;

            ru_0min = floor(ru_0min); ru_0max = ceil(ru_0max);
            for ru_try = ru_0min:ru_0max
                ramp_up_times = [ramp_up_times, ru_try, ru_try];
                ramp_down_times = [ramp_down_times, duration - ru_try - 1, duration - ru_try];
            end
            break;
        end

        % Gradients at corners
        grad_p0 = sign_area * max_slew * ru_max * raster_time + grad_start;
        grad_p1 = sign_area * max_slew * rd_max * raster_time + grad_end;

        if abs(grad_p0) >= max_grad
            limit_option1 = abs(sign_area * max_slew * (ru_limit - 1) * raster_time + grad_start);
            limit_option2 = abs(((ru_limit + flat_time) * sign_area * max_grad + grad_start - grad_p1) / (ru_limit + flat_time));
            if limit_option1 > limit_option2
                grad_p0 = sign_area * max_slew * (ru_limit - 1) * raster_time + grad_start;
                ru_max = floor(abs(sign_area * max_grad - grad_start) / max_slew / raster_time);
                ru_limit = ru_max;
            else
                grad_p0 = sign_area * max_grad;
                ru_max = ceil(abs(sign_area * max_grad - grad_start) / max_slew / raster_time);
                ru_limit = ru_max;
            end
            flat_time = duration - rd_max - ru_max;
        end

        if abs(grad_p1) >= max_grad
            limit_option1 = abs(sign_area * max_slew * (rd_limit - 1) * raster_time + grad_end);
            limit_option2 = abs(((rd_limit + flat_time) * sign_area * max_grad + grad_end - grad_p0) / (rd_limit + flat_time));
            if limit_option1 > limit_option2
                grad_p1 = sign_area * max_slew * (rd_limit - 1) * raster_time + grad_end;
                rd_max = floor(abs(sign_area * max_grad - grad_end) / max_slew / raster_time);
                rd_limit = rd_max;
            else
                grad_p1 = sign_area * max_grad;
                rd_max = ceil(abs(sign_area * max_grad - grad_end) / max_slew / raster_time);
                rd_limit = rd_max;
            end
            flat_time = duration - rd_max - ru_max;
        end

        % Slope too high from p0 to p1
        if abs(grad_p0 - grad_p1) / (flat_time * raster_time) > max_slew
            if abs(grad_p0) < abs(grad_p1)
                grad_p1 = grad_p0 + flat_time * raster_time * max_slew * sign_area;
                rd_max = ceil(abs(grad_end - grad_p1) / max_slew / raster_time);
            else
                grad_p0 = grad_p1 + flat_time * raster_time * max_slew * sign_area;
                ru_max = ceil(abs(grad_start - grad_p0) / max_slew / raster_time);
            end
            flat_time = duration - rd_max - ru_max;
        end

        % Compute area
        area_current = raster_time/2 * ...
            (ru_max * (grad_p0 + grad_start) + ...
             flat_time * (grad_p0 + grad_p1) + ...
             rd_max * (grad_p1 + grad_end));

        if sign(min_dif_area) ~= sign(area - area_current)
            t = round(cumsum([0, ru_max, flat_time, rd_max]) * raster_time, 5);
            if abs(grad_p0) < abs(grad_p1)
                Gtest = grad_p1;
                corner_grad = -(grad_start * ru_max + grad_end * rd_max + Gtest * (flat_time + rd_max) - 2 * area / raster_time) / (flat_time + ru_max);
                amps = [grad_start, corner_grad, Gtest, grad_end];
            else
                Gtest = grad_p0;
                corner_grad = -(grad_start * ru_max + grad_end * rd_max + Gtest * (flat_time + ru_max) - 2 * area / raster_time) / (flat_time + rd_max);
                amps = [grad_start, Gtest, corner_grad, grad_end];
            end

            if any(round(diff(t)/raster_time) < 1)
                continue;
            end

            slew = diff(amps) ./ diff(t);
            if max(abs(slew)) <= max_slew + 1e-5
                times = t;
                amplitudes = amps;
                return;
            end

            % fallback search if slope still too large
            for ru_try = ru_max-1 : duration - rd_max
                for rd_try = rd_max-1 : duration - ru_try
                    flat = duration - ru_try - rd_try;
                    grad_amp = -(ru_try * raster_time * grad_start + rd_try * raster_time * grad_end - 2 * area) / ...
                               ((ru_try + 2 * flat + rd_try) * raster_time);
                    amps = [grad_start, grad_amp, grad_amp, grad_end];
                    t = cumsum([0, ru_try, flat, rd_try]) * raster_time;
                    slew = diff(amps) ./ diff(t);
                    if max(abs(slew)) < max_slew + 1e-5 && max(abs(amps)) < max_grad
                        times = t;
                        amplitudes = amps;
                        return;
                    end
                end
            end
        end

        if abs(area_current - area) < abs(min_dif_area)
            min_dif_area = area - area_current;
        end
    end

    % Fallback to triangle search
    ru_vec = ramp_up_times(:);
    rd_vec = ramp_down_times(:);
    valid = ru_vec .* rd_vec > 0;
    ru_vec = ru_vec(valid); rd_vec = rd_vec(valid);
    flat = duration - ru_vec - rd_vec;
    valid = flat >= 0;
    ru_vec = ru_vec(valid); rd_vec = rd_vec(valid); flat = flat(valid);

    grad_amp = -(ru_vec * raster_time * grad_start + rd_vec * raster_time * grad_end - 2 * area) ./ ...
               ((ru_vec + 2 * flat + rd_vec) * raster_time);

    slew1 = abs(grad_start - grad_amp) ./ (ru_vec * raster_time);
    slew2 = abs(grad_end - grad_amp) ./ (rd_vec * raster_time);

    valid = abs(grad_amp) <= max_grad + 1e-5 & slew1 <= max_slew + 1e-5 & slew2 <= max_slew + 1e-5;
    idx = find(valid, 1);
    if isempty(idx)
        times = []; amplitudes = [];
        return;
    end

    t = cumsum([0, ru_vec(idx), flat(idx), rd_vec(idx)]) * raster_time;
    amps = [grad_start, grad_amp(idx), grad_amp(idx), grad_end];

    times = t;
    amplitudes = amps;
end