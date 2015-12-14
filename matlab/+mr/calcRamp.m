function [kout, success] = calcRamp(k0,kend,varargin)
%
% the aim of joink is to join the points k0 and kend in three-dimensional 
% k-space in minimal time, observing the gradient and slew limits, and the
% gradient strength G0 before k0(:,2) and Gend after kend(:,1)
%
% In the context of a fixed gradient dwell time this is a discrete problem
% with an a priori unknown number of discretization steps. Therefore joink
% tries out the optimization with 0 steps, then 1 step, and so on, until 
% all conditions can be fulfilled, thus yielding a short connection
%
% N.B. The connection found this way is not necessarily always the shortest 
% (there are some counterexamples) but still quite short. Improvements
% possible.
%
% Usage: [kout success] = joink(k0,kend,MaxGrad,MaxSlew,GradDwell,MaxPoints)
%
% [kout]      connecting k-space points without k0 and kend, size = [3,Nt], 
%             where Nt = number of steps between k0 and kend.
%             k-space units: 1/m
%
% [success]   a flag indicating if a solution was found with up to
%             MaxPoints k-space points: success (1), no solution (0)
%
% [k0]        Two preceding points in k-space, size = [3,2]. From these
%             points, the starting gradient will be calculated.
%
% [kend]      Two following points in k-space, size = [3,2]. From these
%             points, the target gradient will be calculated.
%
% [MaxGrad]   maximum total vector gradient strength, size = [1,1]
%    or
% [MaxGrad]   maximum gradient strength per coordinate, size = [3,1]
%
% [MaxSlew]   maximum total vector slew rate, size = [1,1], 
%    or       
% [MaxSlew]   maximum slew rate per coordinate, size = [3,1]
%             all slew units: T/(m*s)
%
% [GradDwell] time between two k-space points, size = [1,1], unit: s
%
% [MaxPoints] maximum number of k-space points to be used in connecting k0
%             with kend. Keep at a reasonable order of magnitude!
%

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'calcRamp';
    parser.addRequired('k0',@isnumeric);
    parser.addRequired('kend',@isnumeric);
    parser.addOptional('system',mr.opts(),@isstruct);
    parser.addParamValue('MaxPoints',500,@isnumeric);
    parser.addParamValue('maxGrad',0,@isnumeric);
    parser.addParamValue('maxSlew',0,@isnumeric);
    
end
parse(parser,k0,kend,varargin{:});
opt = parser.Results;

maxSlew=opt.system.maxSlew;
maxGrad=opt.system.maxGrad;
if opt.maxGrad>0
    maxGrad=opt.maxGrad;
end
if opt.maxSlew>0
    maxSlew=opt.maxSlew;
end

GradRaster = opt.system.gradRasterTime;
MaxPoints = opt.MaxPoints;

SaveRecLimit = get(0,'RecursionLimit'); 
set(0,'RecursionLimit',MaxPoints+10);

% Determine whether we are in componentwise limited mode or in total vector
% limited mode. 
% mode = 0 --> total vector limited
% mode = 1 --> componentwise limited

if     isequal(size(maxGrad),[1 1]) && isequal(size(maxSlew),[1 1])
    mode = 0;
elseif isequal(size(maxGrad),[3 1]) && isequal(size(maxSlew),[3 1])
    mode = 1;
else
    error('Input value MaxGrad or MaxSlew in invalid format.');
end

G0   = (k0  (:,2)-k0  (:,1))/GradRaster;
Gend = (kend(:,2)-kend(:,1))/GradRaster;
k0   = k0  (:,2);
kend = kend(:,1);

success = 0;
kout = zeros(3,0);
UsePoints = 0;                             % first try: connecting directly

while (success == 0) && (UsePoints <= MaxPoints)
    if mode == 0
        if (norm(G0)>maxGrad) || (norm(Gend)>maxGrad) 
            break; 
        end;
        kout = joinleft0(k0,kend,G0,Gend,UsePoints);
    else
        if (abs(G0)>abs(maxGrad)) || (abs(Gend)>abs(maxGrad))
            break; 
        end;        
        kout = joinleft1(k0,kend,G0,Gend,UsePoints);
    end
    UsePoints = UsePoints + 1;
end

set(0,'RecursionLimit',SaveRecLimit);      % set previous recursion limit

% -------------------------------------------------------------------------
function ok = InsideLimits(Grad,Slew)
%
% check if both gradient and slew rates are inside the respective limits
%
if mode == 0
    Grad2 = sum(Grad.^2,1);                % gradient vector norm squared
    Slew2 = sum(Slew.^2,1);                % slew vector norm squared
    ok = (max(Grad2) <= maxGrad^2) && (max(Slew2) <= maxSlew^2);
else
    ok = (sum(max(abs(Grad),[],2) <= maxGrad) == 3) && ...
         (sum(max(abs(Slew),[],2) <= maxSlew) == 3);
end

end % function InsideLimits

% -------------------------------------------------------------------------
function koutleft = joinleft0(k0,kend,G0,Gend,UsePoints)
%
% Add one k-space point close to k0. Gradient and slew limits apply in 
% total vector limited mode.
% 
% Rationale: 
%
%  0. If UsePoints == 0 the recursion stops. If k0 and kend can be joined
%     in one GradDwell time, return success, else return "no success".
%
%  1. Calculate optimal k-space point kopt that would lie on a straight 
%     line of N=UsePoints evenly spaced points to kend. If this kopt can be 
%     reached within gradient and slew limts, kopt is the solution of this
%     function call. 
%
%  2. If kopt cannot be reached, calculate the gradient limited point kgl 
%     closest to kopt. If this point can be reached in one GradDwell time
%     without violating the slew limit, kgl is the solution of this 
%     function call.
%
%  3. If not kgl is not inside the slew limit, the slew limited point  
%     closest to kopt, ksl, is calculated. If ksl is inside the gradient 
%     limit, ksl is the solution.
% 
%  4. If neither kgl nor ksl are possible find the point kglsl closest to
%     kopt that satisfies both limits at the same time. See
%     illustration.fig / illustration.png 
%
%  5. Call joinright0 to obtain the other points starting with a point 
%     next to kend.
%
if UsePoints == 0
    G = [G0 (kend-k0)/GradRaster Gend];
    S = (G(:,2:end)-G(:,1:end-1))/GradRaster;

    koutleft = zeros(3,0);                 % no additional k-space point
    success = InsideLimits(G,S);
   
    return;    
end

dk = (kend-k0)/(UsePoints+1);
kopt = k0+dk;                              % this would be on the direct 
Gopt = (kopt-k0)/GradRaster;                % line
Sopt = (Gopt-G0)/GradRaster;

okGopt = (sum(Gopt.^2,1) <= maxGrad^2);
okSopt = (sum(Sopt.^2,1) <= maxSlew^2);

if okGopt && okSopt
    kLeft = kopt;
else
    a = maxGrad*GradRaster;
    b = maxSlew*GradRaster^2;
    
    dkprol = G0*GradRaster;                 % prolonged point with no change
    dkconn = dk-dkprol;                    % in gradient 
                                           % slew limited closest to kopt
    ksl = k0 + dkprol + dkconn/norm(dkconn)*b;
    Gsl = (ksl-k0)/GradRaster;
    okGsl = (sum(Gsl.^2,1) <= maxGrad^2);
                                           
    kgl = k0 + dk/norm(dk)*a;              % gradient limited closest to
    Ggl = (kgl-k0)/GradRaster;              % kopt
    Sgl = (Ggl-G0)/GradRaster;
    okSgl = (sum(Sgl.^2,1) <= maxSlew^2);

    if okGsl
        kLeft = ksl;
    elseif okSgl
        kLeft = kgl;
    else
        c  = norm(dkprol);
        c1 = (a^2-b^2+c^2)/(2*c);          % if e.g. |Gend|<Gmax, then
        h  = sqrt(a^2-c1^2);               % a^2-c1^2 is always positive
        kglsl = k0 + c1*dkprol/norm(dkprol);
        projondkprol = (kgl*dkprol.') * dkprol/norm(dkprol);
        hdirection = kgl - projondkprol;
        kglsl = kglsl + h*hdirection/norm(hdirection);
        kLeft = kglsl;
    end
end

k = joinright0(kLeft,kend,(kLeft-k0)/GradRaster,Gend,UsePoints-1);

koutleft = [kLeft k];                      % pass along result

end % function joinleft0

% -------------------------------------------------------------------------
function koutright = joinright0(k0,kend,G0,Gend,UsePoints)
%
% Add one k-space point close to kend. Gradient and slew limits apply in 
% total vector limited mode. Rationale see joinleft0.
%
if UsePoints == 0
    G = [G0 (kend-k0)/GradRaster Gend];
    S = (G(:,2:end)-G(:,1:end-1))/GradRaster;

    koutright = zeros(3,0);                % no additional k-space point    
    success = InsideLimits(G,S);
            
    return;    
end

dk = (k0-kend)/(UsePoints+1);
kopt = kend+dk;                            % this would be on the direct 
Gopt = (kend-kopt)/GradRaster;    % line
Sopt = (Gend-Gopt)/GradRaster;

okGopt = (sum(Gopt.^2,1) <= maxGrad^2);
okSopt = (sum(Sopt.^2,1) <= maxSlew^2);

if okGopt && okSopt
    kRight = kopt;
else
    a = maxGrad*GradRaster;
    b = maxSlew*GradRaster^2;
    
    dkprol = -Gend*GradRaster;              % prolonged point with no change
    dkconn = dk-dkprol;                    % in gradient
                                           % slew limited closest to kopt
    ksl = kend + dkprol + dkconn/norm(dkconn)*b;
    Gsl = (kend-ksl)/GradRaster;
    okGsl = (sum(Gsl.^2,1) <= maxGrad^2);
                                           
    kgl = kend + dk/norm(dk)*a;            % gradient limited closest to
    Ggl = (kend-kgl)/GradRaster;  % kopt
    Sgl = (Gend-Ggl)/GradRaster;
    okSgl = (sum(Sgl.^2,1) <= maxSlew^2);

    if okGsl
        kRight = ksl;
    elseif okSgl
        kRight = kgl;
    else
        c  = norm(dkprol);
        c1 = (a^2-b^2+c^2)/(2*c);          % if e.g. |Gend|<Gmax, then
        h  = sqrt(a^2-c1^2);               % a^2-c1^2 is always positive
        kglsl = kend + c1*dkprol/norm(dkprol);
        projondkprol = (kgl*dkprol.') * dkprol/norm(dkprol);
        hdirection = kgl - projondkprol;
        kglsl = kglsl + h*hdirection/norm(hdirection);
        kRight = kglsl;
    end
end

k = joinleft0(k0,kRight,G0,(kend-kRight)/GradRaster,UsePoints-1);

koutright = [k kRight];                    % pass along result

end % function joinright0

% -------------------------------------------------------------------------
function koutleft = joinleft1(k0,kend,G0,Gend,UsePoints)
%
% Add one k-space point close to k0. Gradient and slew limits apply in 
% componentwise limited mode. Rationale is the same as in joinleft0 but
% it's much easier to find the point kglsl in step 4.
% 
if UsePoints == 0
    G = [G0 (kend-k0)/GradRaster Gend];
    S = (G(:,2:end)-G(:,1:end-1))/GradRaster;

    koutleft = zeros(3,0);                 % no additional k-space point
    success = InsideLimits(G,S);
   
    return;    
end

kLeft = zeros(3,1);

dk = (kend-k0)/(UsePoints+1);
kopt = k0+dk;                              % this would be on the direct 
Gopt = (kopt-k0)/GradRaster;               % line
Sopt = (Gopt-G0)/GradRaster;
okGopt = (abs(Gopt) <= maxGrad);
okSopt = (abs(Sopt) <= maxSlew);

dkprol = G0*GradRaster;
dkconn = dk-dkprol;
                                           % slew limited
ksl = k0 + dkprol + sign(dkconn).*maxSlew*GradRaster^2;
Gsl = (ksl-k0)/GradRaster;
okGsl = (abs(Gsl) <= maxGrad);
                                           % gradient limited
kgl = k0 + sign(dk).*maxGrad*GradRaster;
Ggl = (kgl-k0)/GradRaster;
Sgl = (Ggl-G0)/GradRaster;
okSgl = (abs(Sgl) <= maxSlew);

for ii=1:3
    if     (okGopt(ii)==1) && (okSopt(ii)==1)
        kLeft(ii) = kopt(ii);
    elseif (okGsl(ii)==1)
        kLeft(ii) = ksl(ii);
    elseif (okSgl(ii)==1)
        kLeft(ii) = kgl(ii);
    else
        display('Moment mal - hier dürfte ich niemals hinkommen!');        
    end
end

k = joinright1(kLeft,kend,(kLeft-k0)/GradRaster,Gend,UsePoints-1);

koutleft = [kLeft k];                      % pass along result

end % function joinleft1

% -------------------------------------------------------------------------
function koutright = joinright1(k0,kend,G0,Gend,UsePoints)
%
% Add one k-space point close to kend. Gradient and slew limits apply in 
% componentwise limited mode. Rationale is the same as in joinright1
% 
if UsePoints == 0
    G = [G0 (kend-k0)/GradRaster Gend];
    S = (G(:,2:end)-G(:,1:end-1))/GradRaster;

    koutright = zeros(3,0);                % no additional k-space point    
    success = InsideLimits(G,S);
            
    return;    
end

kRight = zeros(3,1);

dk = (k0-kend)/(UsePoints+1);
kopt = kend+dk;                            % this would be on the direct 
Gopt = (kend-kopt)/GradRaster;              % line
Sopt = (Gend-Gopt)/GradRaster;
okGopt = (abs(Gopt) <= maxGrad);
okSopt = (abs(Sopt) <= maxSlew);

dkprol = -Gend*GradRaster;
dkconn = dk-dkprol;
                                           % slew limited 
ksl = kend + dkprol + sign(dkconn).*maxSlew*GradRaster^2;
Gsl = (kend-ksl)/GradRaster;
okGsl = (abs(Gsl) <= maxGrad);
                                           % gradient limited
kgl = kend + sign(dk).*maxGrad*GradRaster;
Ggl = (kend-kgl)/GradRaster;    
Sgl = (Gend-Ggl)/GradRaster;
okSgl = (abs(Sgl) <= maxSlew);

for ii = 1:3
    if     (okGopt(ii)==1) && (okSopt(ii)==1)
        kRight(ii) = kopt(ii);
    elseif (okGsl(ii)==1)
        kRight(ii) = ksl(ii);
    elseif (okSgl(ii)==1)
        kRight(ii) = kgl(ii);
    else
        error('Unknown error. Code should not execute');
    end
end

k = joinleft1(k0,kRight,G0,(kend-kRight)/GradRaster,UsePoints-1);

koutright = [k kRight];                    % pass along result

end % function joinright1

end % main
% =========================================================================
