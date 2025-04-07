function [B, m1, m2, m3] = calcMomentsBtensor(obj, varargin)
%give input arguments of calcB,calcm1,calcm2,calcm3 with true or false and
%Ndummy with 0,1,... as the Number of Dummy scans (scans that should be
%skipped for the calculation)

%m1/m1 are of the structure array(R,3), where R denotes the number of the
%repetition readout (for single shot sequences) with the three orientations
%x,y,z respectively
%B is of a similar structure with array(R,3,3)
%CONSIDER: This so far only works for one slice or rather one measurement PER TR
defaultB = true;
defaultM = false;
defaultD = 0;
%temporalTolerance=1e-9; % 1ns

p=inputParser;
addRequired(p, 'obj');

% addParameter(p, 'calcB', defaultB);
% addParameter(p, 'calcm1', defaultM);
% addParameter(p, 'calcm2', defaultM);
% addParameter(p, 'calcm3', defaultM);
% addParameter(p, 'Ndummy', defaultD);

% Add optional arguments with default values
addOptional(p, 'calcB', defaultB);
addOptional(p, 'calcm1', defaultM);
addOptional(p, 'calcm2', defaultM);
addOptional(p, 'calcm3', defaultM);
addOptional(p, 'Ndummy', defaultD);

% Parse the input
parse(p, obj, varargin{:});

% Assign the results to variables
calcB = p.Results.calcB;
calcm1 = p.Results.calcm1;
calcm2 = p.Results.calcm2;
calcm3 = p.Results.calcm3;
Ndummy = p.Results.Ndummy;

%%
% calcB = true;
% calcm1 = true;
% calcm2 = true;
% calcm3 = true;
% Ndummy = 0;
%% first get gradient shape and timing properties out of the sequence
warning('OFF', 'mr:restoreShape')
[~, ~, ~, ~, t_excitation, t_refocusing, ~, ~, gw_pp] = obj.calculateKspacePP();

R = length(t_excitation); %repetions defined over excitations

gx_pp_l = gw_pp{1};
gy_pp_l = gw_pp{2};
gz_pp_l = gw_pp{3};
%% splitting pp-functions into equal parts
tSeq = [];
t_echo = [];
for i=1:R
    t_echo(end+1)=(2*t_refocusing(i) - t_excitation(i)); % TODO: fixme for double-refocused sequences
    tSeq = [tSeq [t_excitation(i), t_refocusing(i), t_echo(i)]];
end

t1 = gx_pp_l.breaks;
t2 = gy_pp_l.breaks;
t3 = gz_pp_l.breaks;
tn  = unique([t1,t2,t3,tSeq]);

% % debuging / visualization
% tnew = linspace(0,tn(end),10000);
% figure; hold on;
% plot(tnew,ppval(gx_pp_l,tnew));
% plot(tnew,ppval(gy_pp_l,tnew));
% plot(tnew,ppval(gz_pp_l,tnew));

gx_pp_coefs=fillPpCoefs(gx_pp_l,tn);
gx_pp_l = mkpp(tn,gx_pp_coefs);
gy_pp_coefs=fillPpCoefs(gy_pp_l,tn);
gy_pp_l = mkpp(tn,gy_pp_coefs);
gz_pp_coefs=fillPpCoefs(gz_pp_l,tn);
gz_pp_l = mkpp(tn,gz_pp_coefs);

%% splitting pp at multiples of TR
n = (1+Ndummy);
gz_pp = cell(1,R);
gy_pp = cell(1,R);
gx_pp = cell(1,R);
for i = 1:size(gz_pp_l.breaks,2)
    if n == R
        break
    elseif abs(gz_pp_l.breaks(i) - t_excitation(n+1)) == 0 %skipping first excitation
        gz_pp{n-Ndummy} = fnbrk(gz_pp_l,[t_excitation(n) t_excitation(n+1)]);
        gy_pp{n-Ndummy} = fnbrk(gy_pp_l,[t_excitation(n) t_excitation(n+1)]);
        gx_pp{n-Ndummy} = fnbrk(gx_pp_l,[t_excitation(n) t_excitation(n+1)]);
        n = n +1;
    end
end
%add the last block as well
gx_pp{R}= fnbrk(gx_pp_l,[t_excitation(end) gx_pp_l.breaks(end)]);
gy_pp{R}= fnbrk(gy_pp_l,[t_excitation(end) gy_pp_l.breaks(end)]);
gz_pp{R}= fnbrk(gz_pp_l,[t_excitation(end) gz_pp_l.breaks(end)]);


%% considering effects of the rf pulse -> effective gradients
%i.e gradient becomes negative after
%the refocusing pulse
%as gradients are split at excitation, nothing necessary here
for j=1:R
    for i = 1:size(gx_pp{j}.coefs,1)
        if gx_pp{j}.breaks(:,i) == t_refocusing(j+Ndummy)
            gx_pp{j}.coefs(i:end,:) = -gx_pp{j}.coefs(i:end,:);
            gy_pp{j}.coefs(i:end,:) = -gy_pp{j}.coefs(i:end,:);
            gz_pp{j}.coefs(i:end,:) = -gz_pp{j}.coefs(i:end,:);
        end
    end
end
%% B-Tensor calculation
B = zeros(R,3,3);
if calcB
    
    qz_pp = cell(1,R);
    qy_pp = cell(1,R);
    qx_pp = cell(1,R);
   
    for i = 1:R
        qx_pp{i}=fnint(gx_pp{i});
        qy_pp{i}=fnint(gy_pp{i});
        qz_pp{i}=fnint(gz_pp{i});

        qx_pp{i}.coefs = qx_pp{i}.coefs*2*pi;
        qy_pp{i}.coefs = qy_pp{i}.coefs*2*pi;
        qz_pp{i}.coefs = qz_pp{i}.coefs*2*pi;
    end

    % tnew = linspace(0,gx_pp_l.breaks(end),1000);
    % figure;
    % plot(tnew,ppval(qz_pp{2},tnew),'LineWidth',2,'Color','r');
    % % hold on;
    % % plot(tnew,ppval(qy_pp{1},tnew),'LineWidth',2,'Color','b');
    % %plot(tnew,ppval(qz_pp{2},tnew),'LineWidth',2,'Color','y');
    % hold off
    % title('wave vector')
    % legend('kx','ky','kz')
    % xlabel('Time [ms]')
    % ylabel('Wave Vector [1/m]')
    % grid;

    for m = 1:R
        q_pp{1} = qx_pp{m};
        q_pp{2} = qy_pp{m};
        q_pp{3} = qz_pp{m};
        for i = 1:3
            for j = 1:3
                coefs=zeros(size(qx_pp{m}.coefs,1),5);
                for k=1:size(qx_pp{m}.coefs,1)
                    coefs(k,:)=conv(q_pp{i}.coefs(k,:),q_pp{j}.coefs(k,:));
                end
                Bpp_div = mkpp(qx_pp{m}.breaks,coefs);
                Bpp = fnint(Bpp_div);
                %evaluate b-value at TE
                B(m,i,j) = ppval(Bpp,tSeq(3 + 3*(m-1+Ndummy)));
            end
        end
    end
end

%% calculate first moment m1
m1 = zeros(R,3);

if calcm1

    m1z_pp = cell(1,R);
    m1y_pp = cell(1,R);
    m1x_pp = cell(1,R);
    for i = 1:R
        %defining a new piecewise polynomial aquivalent to g(t)=t
        t_pp_coefs = zeros(length(gx_pp{i}.breaks)-1,2);
        t_pp_coefs(:,1) = 1;
        t_pp_coefs(:,2) = gx_pp{i}.breaks(1:(end-1))-gx_pp{i}.breaks(1);

        new_cx = zeros(size(gx_pp{i}.coefs,1),size(gx_pp{i}.coefs,2)+1);
        new_cy = zeros(size(gy_pp{i}.coefs,1),size(gy_pp{i}.coefs,2)+1);
        new_cz = zeros(size(gz_pp{i}.coefs,1),size(gz_pp{i}.coefs,2)+1);

        for k=1:size(t_pp_coefs,1)
            new_cx(k,:)=conv(gx_pp{i}.coefs(k,:),t_pp_coefs(k,:));
            new_cy(k,:)=conv(gy_pp{i}.coefs(k,:),t_pp_coefs(k,:));
            new_cz(k,:)=conv(gz_pp{i}.coefs(k,:),t_pp_coefs(k,:));
        end
        
        tgx_pp = mkpp(gx_pp{i}.breaks, new_cx);
        tgy_pp = mkpp(gy_pp{i}.breaks, new_cy);
        tgz_pp = mkpp(gz_pp{i}.breaks, new_cz);

        m1x_pp{i}=fnint(tgx_pp);
        m1y_pp{i}=fnint(tgy_pp);
        m1z_pp{i}=fnint(tgz_pp);
        m1x_pp{i}.coefs = m1x_pp{i}.coefs*2*pi;
        m1y_pp{i}.coefs = m1y_pp{i}.coefs*2*pi;
        m1z_pp{i}.coefs = m1z_pp{i}.coefs*2*pi;
    end

    %as an integration from 0 to TE, I will evaluate the function at TE to
    %gain the first order vector
    
    for m=1:R
        m1(m,1) = ppval(m1x_pp{m},tSeq(3 + 3*(m-1+Ndummy)));%-ppval(m1x_pp{m},m1x_pp{m}.breaks(1));
        m1(m,2) = ppval(m1y_pp{m},tSeq(3 + 3*(m-1+Ndummy)));%-ppval(m1y_pp{m},m1y_pp{m}.breaks(1));
        m1(m,3) = ppval(m1z_pp{m},tSeq(3 + 3*(m-1+Ndummy)));%-ppval(m1z_pp{m},m1z_pp{m}.breaks(1));
    end

end
%% calculate second moment m2
m2 = zeros(R,3);

if calcm2

    m2z_pp = cell(1,R);
    m2y_pp = cell(1,R);
    m2x_pp = cell(1,R);
    for i = 1:R
        t2_pp_coefs = zeros(length(gx_pp{i}.breaks)-1,3);
        t2_pp_coefs(:,1) = 1;
        t2_pp_coefs(:,2) = 2*(gx_pp{i}.breaks(1:(end-1))-gx_pp{i}.breaks(1));
        t2_pp_coefs(:,3) = (gx_pp{i}.breaks(1:(end-1))-gx_pp{i}.breaks(1)).^2;
        
        new_cx = zeros(size(gx_pp{i}.coefs,1),size(gx_pp{i}.coefs,2)+2);
        new_cy = zeros(size(gy_pp{i}.coefs,1),size(gy_pp{i}.coefs,2)+2);
        new_cz = zeros(size(gz_pp{i}.coefs,1),size(gz_pp{i}.coefs,2)+2);

        for k=1:size(t2_pp_coefs,1)
            new_cx(k,:)=conv(gx_pp{i}.coefs(k,:),t2_pp_coefs(k,:));
            new_cy(k,:)=conv(gy_pp{i}.coefs(k,:),t2_pp_coefs(k,:));
            new_cz(k,:)=conv(gz_pp{i}.coefs(k,:),t2_pp_coefs(k,:));
        end
        
        tgx_pp = mkpp(gx_pp{i}.breaks, new_cx);
        tgy_pp = mkpp(gy_pp{i}.breaks, new_cy);
        tgz_pp = mkpp(gz_pp{i}.breaks, new_cz);
        
        m2x_pp{i}=fnint(tgx_pp);
        m2y_pp{i}=fnint(tgy_pp);
        m2z_pp{i}=fnint(tgz_pp);
        m2x_pp{i}.coefs = m2x_pp{i}.coefs*2*pi;
        m2y_pp{i}.coefs = m2y_pp{i}.coefs*2*pi;
        m2z_pp{i}.coefs = m2z_pp{i}.coefs*2*pi;
    end

    %as an integration from 0 to TE, I will evaluate the function at TE to
    %gain the first order vector
    
    for m=1:R
        m2(m,1) = ppval(m2x_pp{m},tSeq(3 + 3*(m-1+Ndummy)))-ppval(m2x_pp{m},m2x_pp{m}.breaks(1));
        m2(m,2) = ppval(m2y_pp{m},tSeq(3 + 3*(m-1+Ndummy)))-ppval(m2y_pp{m},m2y_pp{m}.breaks(1));
        m2(m,3) = ppval(m2z_pp{m},tSeq(3 + 3*(m-1+Ndummy)))-ppval(m2z_pp{m},m2z_pp{m}.breaks(1));
    end

end


%% calc third moment m3:
m3 = zeros(R,3);

if calcm3

    m3z_pp = cell(1,R);
    m3y_pp = cell(1,R);
    m3x_pp = cell(1,R);
    for i = 1:R
        t3_pp_coefs = zeros(length(gx_pp{i}.breaks)-1,4);
        t3_pp_coefs(:,1) = 1;
        t3_pp_coefs(:,2) = 3*(gx_pp{i}.breaks(1:(end-1))-gx_pp{i}.breaks(1));
        t3_pp_coefs(:,3) = 3*(gx_pp{i}.breaks(1:(end-1))-gx_pp{i}.breaks(1)).^2;
        t3_pp_coefs(:,4) = (gx_pp{i}.breaks(1:(end-1))-gx_pp{i}.breaks(1)).^3;
        % t3_pp_breaks = gx_pp{i}.breaks;
        % t3_pp = mkpp(t3_pp_breaks,t3_pp_coefs);
        % tnew = linspace(0,gx_pp_l.breaks(end),1000);
        % figure;
        % plot(tnew,ppval(t3_pp,tnew),'LineWidth',2,'Color','r');
        % hold on
        % plot(tnew,(tnew-gx_pp{i}.breaks(1)).^3,'Color','b');
        % hold off

        new_cx = zeros(size(gx_pp{i}.coefs,1),size(gx_pp{i}.coefs,2)+3);
        new_cy = zeros(size(gy_pp{i}.coefs,1),size(gy_pp{i}.coefs,2)+3);
        new_cz = zeros(size(gz_pp{i}.coefs,1),size(gz_pp{i}.coefs,2)+3);

        for k=1:size(t3_pp_coefs,1)
            new_cx(k,:)=conv(gx_pp{i}.coefs(k,:),t3_pp_coefs(k,:));
            new_cy(k,:)=conv(gy_pp{i}.coefs(k,:),t3_pp_coefs(k,:));
            new_cz(k,:)=conv(gz_pp{i}.coefs(k,:),t3_pp_coefs(k,:));
        end
        
        tgx_pp = mkpp(gx_pp{i}.breaks, new_cx);
        tgy_pp = mkpp(gy_pp{i}.breaks, new_cy);
        tgz_pp = mkpp(gz_pp{i}.breaks, new_cz);
        
        m3x_pp{i}=fnint(tgx_pp);
        m3y_pp{i}=fnint(tgy_pp);
        m3z_pp{i}=fnint(tgz_pp);
        m3x_pp{i}.coefs = m3x_pp{i}.coefs*2*pi;
        m3y_pp{i}.coefs = m3y_pp{i}.coefs*2*pi;
        m3z_pp{i}.coefs = m3z_pp{i}.coefs*2*pi;
    end

    %as an integration from 0 to TE, I will evaluate the function at TE to
    %gain the first order vector
    
    for m=1:R
        m3(m,1) = ppval(m3x_pp{m},tSeq(3 + 3*(m-1+Ndummy)))-ppval(m3x_pp{m},m3x_pp{m}.breaks(1));
        m3(m,2) = ppval(m3y_pp{m},tSeq(3 + 3*(m-1+Ndummy)))-ppval(m3y_pp{m},m3y_pp{m}.breaks(1));
        m3(m,3) = ppval(m3z_pp{m},tSeq(3 + 3*(m-1+Ndummy)))-ppval(m3z_pp{m},m3z_pp{m}.breaks(1));
    end

end
end
%% functions 
function pp1_coefs=fillPpCoefs(pp1,xn)
    idx1 = slookup(xn(1:end-1),pp1.breaks(1:end-1));
    pp1_coefs = zeros(length(xn)-1,pp1.order);
    for i=1:size(pp1_coefs,1)
        if idx1(i)>0
            % simple copy
            pp1_coefs(i,:)=pp1.coefs(idx1(i),:);
        elseif i > 1
            % copy from the left with a reference shift (else leave at 0)
            for k=0:(pp1.order-1)
                for l=0:k
                    pp1_coefs(i,end-l) = pp1_coefs(i,end-l) + pp1_coefs(i-1,end-k)*nchoosek(k,l)*(xn(i)-xn(i-1))^(k-l);
                end
            end
        end
    end
end

function idx=slookup(what,where)
% finds indices of values given by the sorted vector 'what' in the sorted vector 'where'  
% for failing searches indices of 0 are returned
    idx=zeros(size(what));
    wb=1; % where bound
    for c=1:length(what)
        i=find(what(c)==where(wb:end),1);
        if isempty(i), continue; end
        idx(c)=wb+i-1;
        wb=wb+i;
    end
end

function idx=sintlookup(what,where)
% finds indices of intervals to which the sorted vector 'what' belongs in the sorted vector 'where'  
% for failing searches indices of 0 are returned
    idx=zeros(size(what));
    wb=1; % where bound
    for c=1:length(what)
        i=find(what(c)>=where(wb:end));
        if isempty(i), continue; end
        idx(c)=wb+i(end)-1;
        wb=idx(c);
    end
end