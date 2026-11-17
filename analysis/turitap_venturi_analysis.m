function out = turitap_venturi_analysis()
% TURITAP_VENTURI_ANALYSIS  First-pass 1-D hydraulic model of the TuriTap
% Venturi on a handpump (vertical, free discharge to atmosphere).
%
%   Run:  out = turitap_venturi_analysis;
%
%   Works in MATLAB (R2016b+) and GNU Octave. No toolboxes needed.
%
%   What it does
%     Fig 1  2-D cross-sections: current part (from "1-2 Venturi 03.step")
%            and one parametric handpump candidate.
%     Fig 2  Required pump-head water level and throat pressure vs flow.
%     Fig 3  Pressure along the axis at low flow and peak flow.
%     Fig 4  Design map over throat diameter and exit diameter.
%     Fig 5  Pump-head water level through piston strokes (tank buffering).
%
%   Model: steady, single-phase, pipe-full Bernoulli + Crane TP-410 cone
%   losses + Darcy friction. Exit kinetic energy is lost (free jet).
%   Expect roughly +/-25 % on losses. It CANNOT see air ingestion,
%   diffuser separation or jet detachment - it flags when the pipe-full
%   assumption is violated (inlet gauge head < 0) instead.
%
%   Sign conventions: z measured DOWN from the inlet plane [m];
%   pressures are gauge; heads are metres of water.

% ---------------------------------------------------------------- inputs
fl.rho = 997;  fl.nu = 8.9e-7;  fl.g = 9.81;  fl.pv = 3170;  fl.patm = 101325;

op.Qlow  = 12;      % L/min  dosing check point
op.Qpeak = 50;      % L/min  peak piston-stroke flow
op.Qplot = linspace(4, 60, 57);

% --- ASSUMPTIONS TO CONFIRM (see questions) ---------------------------
site.hFree = 0.15;  % m  allowable water rise above the Venturi inlet plane
site.Kent  = 1.1;   % -  upstream loss on inlet-bore velocity: reducer + threaded 3/4" 90 deg elbow (0.4 for a 1" sweep elbow)
site.Atank = 0.015; % m2 free-surface area in the pump head that stores water
site.Lmax  = 0.1165;% m  overall length of the candidate (0.1165 = same as current part)

pump.Qavg  = 25;    % L/min  mean delivery (bucket test: 20-25, up to 30)
pump.Qpk   = 50;    % L/min  instantaneous peak
pump.spm   = 40;    % strokes per minute

mdl.eps       = 5e-6;  % m  wall roughness (SPI B-3 moulded HDPE)
mdl.kMinor    = 1.0;   % multiplier on cone losses (sensitivity knob)
mdl.alphaExit = 1.05;  % kinetic-energy factor of the exit jet

% ------------------------------------------------------------ geometries
cur = geomCurrent();

p.Din   = 21.0;   % mm  inlet bore
p.Lin   = 15.0;   % mm  straight inlet socket
p.thc   = 12.0;   % deg converging half-angle
p.Dt    = 14.0;   % mm  throat diameter
p.Lt    = 14.0;   % mm  cylindrical throat length
p.thd   = 5.0;    % deg diffuser half-angle
p.Dex   = 26.0;   % mm  diffuser exit / tailpipe bore
p.Ltail = [];     % mm  [] = fill up to site.Lmax
p.R1    = 1.0;    % throat entry fillet radius / Dt
p.R2    = 1.0;    % throat exit  fillet radius / Dt
p.Lmax  = site.Lmax*1e3;
cand = geomParam(p);

% ------------------------------------------------------- flow-rate sweep
Sc = sweepQ(cur,  op.Qplot, fl, mdl, site);
Sn = sweepQ(cand, op.Qplot, fl, mdl, site);

fprintf('\n=== Steady pipe-full results ===\n');
reportPoint('CURRENT  ', cur,  [10 12 20 30 50 60], fl, mdl, site);
reportPoint('CANDIDATE', cand, [10 12 20 30 50 60], fl, mdl, site);

Qfc = qAtLevel(cur,  0,          fl, mdl, site);
Qoc = qAtLevel(cur,  site.hFree, fl, mdl, site);
Qfn = qAtLevel(cand, 0,          fl, mdl, site);
Qon = qAtLevel(cand, site.hFree, fl, mdl, site);
fprintf('\nPipe-full operating window (steady):\n');
fprintf('  CURRENT  : inlet starves below %4.1f L/min, overflows above %4.1f L/min\n', Qfc, Qoc);
fprintf('  CANDIDATE: inlet starves below %4.1f L/min, overflows above %4.1f L/min\n', Qfn, Qon);

% --------------------------------------------------- pump-head transient
Tc = pumpHead(cur,  pump, fl, mdl, site);
Tn = pumpHead(cand, pump, fl, mdl, site);
fprintf('\nPump-head level with stroke buffering (Atank = %.0f cm2):\n', site.Atank*1e4);
fprintf('  CURRENT  : peak level %5.1f cm above inlet (limit %.0f cm)\n', Tc.hmax*100, site.hFree*100);
fprintf('  CANDIDATE: peak level %5.1f cm above inlet (limit %.0f cm)\n', Tn.hmax*100, site.hFree*100);

% -------------------------------------------------------- design options
designOptions(p, fl, mdl, site);

% ------------------------------------------------------------ design map
Dt  = 11:0.5:17;   Dex = 17:1:32;
[HH, PP, QF] = designMap(Dt, Dex, p, op, fl, mdl, site);

% ----------------------------------------------------------------- plots
figure(1); clf;
subplot(1,2,1); drawSection(cur,  'Current TuriTap (from STEP)');
subplot(1,2,2); drawSection(cand, 'Handpump candidate (parametric)');

figure(2); clf;
subplot(2,1,1);
plot(op.Qplot, Sc.hTank*100, 'r-', op.Qplot, Sn.hTank*100, 'b-', 'LineWidth', 1.6); hold on;
plot(op.Qplot([1 end]), site.hFree*100*[1 1], 'k--');
plot(op.Qplot([1 end]), [0 0], 'k:');
grid on; xlabel('Flow  [L/min]'); ylabel('Water level above inlet  [cm]');
title('Pump-head level needed to push the flow through');
legend('current', 'candidate', 'overflow limit', 'inlet starved below this', 'Location', 'northwest');
subplot(2,1,2);
plot(op.Qplot, Sc.pThroat/1e3, 'r-', op.Qplot, Sn.pThroat/1e3, 'b-', 'LineWidth', 1.6);
grid on; xlabel('Flow  [L/min]'); ylabel('Throat gauge pressure  [kPa]');
title('Throat pressure (pipe-full assumption)');
legend('current', 'candidate', 'Location', 'southwest');

figure(3); clf;
Q2 = [op.Qlow op.Qpeak];
for k = 1:2
  subplot(1,2,k);
  a = hydraulics(cur,  Q2(k), fl, mdl, site);
  b = hydraulics(cand, Q2(k), fl, mdl, site);
  plot(a.pProf/1e3, cur.z*1e3, 'r-', b.pProf/1e3, cand.z*1e3, 'b-', 'LineWidth', 1.6);
  set(gca, 'YDir', 'reverse'); grid on;
  xlabel('Gauge pressure  [kPa]'); ylabel('Distance below inlet  [mm]');
  title(sprintf('Axial pressure at %d L/min', Q2(k)));
  legend('current', 'candidate', 'Location', 'southeast');
end

figure(4); clf;
subplot(1,3,1); [c,h] = contour(Dt, Dex, HH'*100, [0 5 10 15 20 30 50 80]); clabel(c,h);
xlabel('Throat D_t  [mm]'); ylabel('Exit D_{ex}  [mm]'); grid on;
title(sprintf('Level above inlet at %d L/min  [cm]', op.Qpeak));
subplot(1,3,2); [c,h] = contour(Dt, Dex, PP'/1e3, -6:0.5:0); clabel(c,h);
xlabel('Throat D_t  [mm]'); ylabel('Exit D_{ex}  [mm]'); grid on;
title(sprintf('Throat pressure at %d L/min  [kPa]', op.Qlow));
subplot(1,3,3); [c,h] = contour(Dt, Dex, QF', 10:5:60); clabel(c,h);
xlabel('Throat D_t  [mm]'); ylabel('Exit D_{ex}  [mm]'); grid on;
title('Lowest pipe-full flow  [L/min]');

figure(5); clf;
plot(Tc.t, Tc.Qp, 'k:', Tc.t, Tc.Qv, 'r-', Tn.t, Tn.Qv, 'b-', 'LineWidth', 1.3); hold on;
ylabel('Flow  [L/min]'); xlabel('Time  [s]'); grid on;
title('Piston pulse in, Venturi flow out (last 3 strokes)');
legend('pump delivery', 'through current', 'through candidate', 'Location', 'northeast');
xlim([Tc.t(end)-3*60/pump.spm, Tc.t(end)]);

out.cur = cur; out.cand = cand; out.Sc = Sc; out.Sn = Sn;
out.Tc = Tc; out.Tn = Tn; out.map.Dt = Dt; out.map.Dex = Dex;
out.map.hTank = HH; out.map.pThroat = PP; out.map.Qfull = QF;
end

% ======================================================================
function g = geomCurrent()
% Bore of "1-2 Venturi 03" read from the STEP surfaces. [z r] in mm,
% z down from the inlet face. Drafted (1 deg) sections included.
v = [  0.00  10.772     % inlet face, dia 21.54
      15.60  10.500     % end of inlet socket, dia 21.00
      34.95   6.350     % end of 12.1 deg cone, dia 12.70
      47.50   6.131     % minimum bore, dia 12.26 (mould parting line)
      60.35   6.350     % start of diffuser, dia 12.70
      79.40   7.938     % end of 4.76 deg diffuser, dia 15.88
     116.50   8.585 ];  % outlet face, dia 17.17
g.vert  = v*1e-3;
g.zPort = 47.7e-3;      % side port axis (1/4-28 UNF, 1.67 mm orifice, 5 deg tilt)
g.dPort = 1.67e-3;
g = finishGeom(g, [], []);
end

function g = geomParam(p)
% Parametric handpump Venturi. Lengths in mm, angles in deg.
Lc = (p.Din - p.Dt)/2/tand(p.thc);
Ld = (p.Dex - p.Dt)/2/tand(p.thd);
Lfix = p.Lin + Lc + p.Lt + Ld;
if isempty(p.Ltail), Ltail = max(p.Lmax - Lfix, 0); else, Ltail = p.Ltail; end
z = cumsum([0 p.Lin Lc p.Lt Ld]);
r = [p.Din p.Din p.Dt p.Dt p.Dex]/2;
if Ltail > 0, z(end+1) = z(end) + Ltail; r(end+1) = p.Dex/2; end
g.vert  = [z(:) r(:)]*1e-3;
g.zPort = (z(3) + p.Lt/2)*1e-3;
g.dPort = 1.67e-3;
g.Ltot  = z(end)*1e-3;
g = finishGeom(g, [3 4], [p.R1 p.R2]*p.Dt*1e-3);
end

function g = finishGeom(g, idx, R)
% Build the fine wall profile (with optional fillets at vertices idx).
P = g.vert;
Q = P(1,:);
for i = 2:size(P,1)-1
  k = find(idx == i, 1);
  if isempty(k)
    Q(end+1,:) = P(i,:); %#ok<AGROW>
  else
    Q = [Q; filletArc(P(i-1,:), P(i,:), P(i+1,:), R(k))]; %#ok<AGROW>
  end
end
Q(end+1,:) = P(end,:);
zf = linspace(P(1,1), P(end,1), 600)';
[zu, iu] = unique(Q(:,1));
g.z = zf;  g.r = interp1(zu, Q(iu,2), zf);
g.Ltot = P(end,1);
end

function A = filletArc(Pa, Pb, Pc, R)
% Round the corner at Pb with radius R (tangent to both legs).
u1 = (Pa - Pb)/norm(Pa - Pb);  u2 = (Pc - Pb)/norm(Pc - Pb);
phi = acos(max(-1, min(1, dot(u1, u2))));
if abs(phi - pi) < 1e-6, A = Pb; return; end
t = min(R/tan(phi/2), 0.45*min(norm(Pa - Pb), norm(Pc - Pb)));
R = t*tan(phi/2);
T1 = Pb + u1*t;  T2 = Pb + u2*t;
b  = (u1 + u2)/norm(u1 + u2);
C  = Pb + b*R/sin(phi/2);
a1 = atan2(T1(2) - C(2), T1(1) - C(1));
a2 = atan2(T2(2) - C(2), T2(1) - C(1));
da = mod(a2 - a1 + pi, 2*pi) - pi;
a  = a1 + da*linspace(0, 1, 15)';
A  = [C(1) + R*cos(a), C(2) + R*sin(a)];
end

% ======================================================================
function R = hydraulics(g, Qlpm, fl, mdl, site)
% Pipe-full energy balance, marched upward from the atmospheric exit.
Q  = Qlpm/60000;  g0 = fl.g;
z  = g.vert(:,1);  r = g.vert(:,2);  n = numel(z);
V  = Q./(pi*r.^2);
hs = zeros(n-1, 1);
for i = 1:n-1
  L  = z(i+1) - z(i);
  rs = min(r(i), r(i+1));  rl = max(r(i), r(i+1));
  b2 = (rs/rl)^2;  th = atan((rl - rs)/L);  Vs = Q/(pi*rs^2);
  if r(i+1) < r(i)                       % contraction (Crane TP-410)
    if th <= pi/8, K = 0.8*sin(th)*(1 - b2);
    else,          K = 0.5*sqrt(sin(th))*(1 - b2); end
  else                                   % expansion / diffuser
    if th <= pi/8, K = 2.6*sin(th)*(1 - b2)^2;
    else,          K = (1 - b2)^2; end
  end
  rr = linspace(r(i), r(i+1), 21);  rm = 0.5*(rr(1:end-1) + rr(2:end));
  Vm = Q./(pi*rm.^2);  Dm = 2*rm;
  f  = fric(Vm.*Dm/fl.nu, mdl.eps./Dm);
  hf = sum(f.*(L/20)./Dm.*Vm.^2/(2*g0));
  hs(i) = mdl.kMinor*K*Vs^2/(2*g0) + hf;
end
E = zeros(n, 1);                         % total head, datum = exit plane
E(n) = mdl.alphaExit*V(n)^2/(2*g0);
for i = n-1:-1:1, E(i) = E(i+1) + hs(i); end
elev = z(end) - z;
pv   = fl.rho*g0*(E - elev - V.^2/(2*g0));

Ef = interp1(z, E, g.z);  Vf = Q./(pi*g.r.^2);
R.pProf   = fl.rho*g0*(Ef - (z(end) - g.z) - Vf.^2/(2*g0));
R.pInlet  = pv(1);                       % static gauge at inlet plane [Pa]
R.hInlet  = pv(1)/(fl.rho*g0);           % same, m of water
R.hTank   = E(1) - elev(1) + site.Kent*V(1)^2/(2*g0);  % still-water level above inlet [m]
[R.pThroat, it] = min(R.pProf);
R.zThroat = g.z(it);
R.pPort   = interp1(g.z, R.pProf, g.zPort);
R.hLoss   = sum(hs);                     % internal losses [m]
R.hExitKE = E(n);                        % exit jet energy thrown away [m]
R.Vthroat = max(Vf);
R.cavMargin = (fl.patm + R.pThroat - fl.pv)/(fl.rho*g0);  % m, >0 = OK
ie = find(diff(g.r) > 1e-9, 1);          % diffuser recovery efficiency
if isempty(ie), R.etaDiff = NaN; else
  Vt = max(Vf);  Ve = Vf(end);
  R.etaDiff = 1 - sum(hs(z(1:end-1) >= g.z(it) - 1e-9))/((Vt^2 - Ve^2)/(2*g0));
end
end

function f = fric(Re, rel)
% Darcy friction factor: laminar / Haaland, blended.
f = 64./max(Re, 1);
t = Re > 2300;
f(t) = (-1.8*log10((rel(t)/3.7).^1.11 + 6.9./Re(t))).^-2;
end

function S = sweepQ(g, Qv, fl, mdl, site)
for k = numel(Qv):-1:1
  R = hydraulics(g, Qv(k), fl, mdl, site);
  S.hTank(k) = R.hTank;  S.hInlet(k) = R.hInlet;  S.pThroat(k) = R.pThroat;
  S.pPort(k) = R.pPort;  S.etaDiff(k) = R.etaDiff;
end
end

function Q = qAtLevel(g, h, fl, mdl, site)
% Flow [L/min] at which the still-water level sits h above the inlet.
Q = fzero(@(q) getfield(hydraulics(g, q, fl, mdl, site), 'hTank') - h, [0.5 400]); %#ok<GFLD>
end

function reportPoint(name, g, Qs, fl, mdl, site)
fprintf('%s  L = %.1f mm\n', name, g.Ltot*1e3);
fprintf('   Q[L/min]  level[cm]  p_in[cm]  p_throat[kPa]  V_t[m/s]  loss[cm]  exitKE[cm]  eta_diff  note\n');
for q = Qs
  R = hydraulics(g, q, fl, mdl, site);  s = '';
  if R.hTank < 0,          s = 'inlet starved (air ingestion)'; end
  if R.hTank > site.hFree, s = 'OVERFLOW'; end
  fprintf('   %6.0f   %8.1f  %8.1f  %12.2f  %8.2f  %8.1f  %9.1f  %8.2f  %s\n', q, ...
    R.hTank*100, R.hInlet*100, R.pThroat/1e3, R.Vthroat, R.hLoss*100, R.hExitKE*100, R.etaDiff, s);
end
end

% ======================================================================
function T = pumpHead(g, pump, fl, mdl, site)
% Water level in the pump head: Atank*dh/dt = Qpump(t) - Qventuri(h).
% Pump pulse = half-sine of peak Qpk, widened so the mean equals Qavg.
Tper = 60/pump.spm;
duty = min(pump.Qavg*pi/(2*pump.Qpk), 1);
qg = linspace(0.5, 300, 400);
hg = arrayfun(@(q) getfield(hydraulics(g, q, fl, mdl, site), 'hTank'), qg); %#ok<GFLD>
Q0 = interp1(hg, qg, 0);                 % Venturi capacity with level at inlet
dt = 2e-3;  t = (0:dt:120*Tper)';  n = numel(t);   % long enough to reach a repeating cycle
h = zeros(n, 1);  Qv = zeros(n, 1);
ph = mod(t, Tper)/(duty*Tper);
Qp = pump.Qpk*sin(pi*ph).*(ph < 1);
for i = 1:n-1
  if h(i) <= 0 && Qp(i) <= Q0
    Qv(i) = Qp(i);  h(i+1) = 0;          % everything passes, pipe not full
  else
    Qv(i) = interp1(hg, qg, h(i), 'linear', 'extrap');
    h(i+1) = max(h(i) + dt*(Qp(i) - Qv(i))/60000/site.Atank, 0);
  end
end
Qv(n) = Qv(n-1);
last = t > t(end) - Tper;
T.t = t;  T.h = h;  T.Qp = Qp;  T.Qv = Qv;  T.hmax = max(h(last));
T.fracFull = mean(h(last) > 0);          % share of the stroke running pipe-full
end

function [HH, PP, QF] = designMap(Dt, Dex, p, op, fl, mdl, site)
HH = nan(numel(Dt), numel(Dex));  PP = HH;  QF = HH;
for i = 1:numel(Dt)
  for j = 1:numel(Dex)
    if Dex(j) <= Dt(i) + 1, continue; end
    q = p;  q.Dt = Dt(i);  q.Lt = Dt(i);  q.Dex = Dex(j);  q.Ltail = [];
    Lfix = q.Lin + (q.Din - q.Dt)/2/tand(q.thc) + q.Lt + (q.Dex - q.Dt)/2/tand(q.thd);
    if Lfix > q.Lmax, continue; end      % does not fit in the length limit
    g = geomParam(q);
    a = hydraulics(g, op.Qpeak, fl, mdl, site);
    b = hydraulics(g, op.Qlow,  fl, mdl, site);
    HH(i,j) = a.hTank;  PP(i,j) = b.pThroat;
    QF(i,j) = qAtLevel(g, 0, fl, mdl, site);
  end
end
end

function designOptions(p, fl, mdl, site)
% For each overall length and throat diameter, find the exit diameter,
% diffuser angle and cone angle that give the lowest water level at
% 30 L/min, then report the level at 20-50 L/min.
%   Q_on  = flow where throat velocity head equals that of the current
%           part at 6 L/min (its observed dosing threshold) -> scales Dt^2.
%   port  = side-port orifice that keeps the same chlorine:water ratio.
Dt0 = 12.26;  Qon0 = 6;  dPort0 = 1.67;  Qs = [20 25 30 40 50];
fprintf('\n=== Design options (inlet bore %.0f mm, entry K = %.2f) ===\n', p.Din, site.Kent);
for Ltot = [116.5 150 200 250]
  fprintf('Overall length %.1f mm\n', Ltot);
  fprintf('   Dt   Dex  th_d  th_c | level above inlet [cm] at Q =  20     25     30     40     50 | Q_on  port\n');
  for Dt = 12:1:18
    best = inf;  gb = [];  qb = [];
    for Dex = Dt+2:1:32
      for thd = [3.5 5 6.5]
        for thc = [12 20]
          q = p;  q.Dt = Dt;  q.Lt = Dt;  q.Dex = Dex;  q.thd = thd;  q.thc = thc;  q.Lmax = Ltot;  q.Ltail = [];
          Lfix = q.Lin + (q.Din - Dt)/2/tand(thc) + q.Lt + (Dex - Dt)/2/tand(thd);
          if Lfix > Ltot, continue; end
          g = geomParam(q);
          R = hydraulics(g, 30, fl, mdl, site);
          if R.hTank < best, best = R.hTank;  gb = g;  qb = q; end
        end
      end
    end
    if isempty(gb), continue; end
    lv = arrayfun(@(x) getfield(hydraulics(gb, x, fl, mdl, site), 'hTank'), Qs)*100; %#ok<GFLD>
    fprintf('  %4.1f  %3.0f  %4.1f  %4.0f |                               %6.1f %6.1f %6.1f %6.1f %6.1f | %4.1f  %4.2f\n', ...
      Dt, qb.Dex, qb.thd, qb.thc, lv, Qon0*(Dt/Dt0)^2, dPort0*Dt/Dt0);
  end
end
end

% ======================================================================
function drawSection(g, ttl)
% 2-D cross-section: bore in white, 2.5 mm nominal wall in grey.
z = g.z*1e3;  r = g.r*1e3;  w = 2.5;
fill([r + w; flipud(r)], [z; flipud(z)], [0.75 0.75 0.75]); hold on;
fill(-[r + w; flipud(r)], [z; flipud(z)], [0.75 0.75 0.75]);
plot([0 0], [z(1) - 5, z(end) + 5], 'k-.');
rp = interp1(z, r, g.zPort*1e3);
plot(-[rp, rp + w + 6], g.zPort*1e3*[1 1], 'b-', 'LineWidth', 2);
text(-(rp + w + 7), g.zPort*1e3, 'Cl port', 'HorizontalAlignment', 'right', 'Color', 'b');
v = g.vert*1e3;
for i = 1:size(v,1)
  plot(v(i,2) + w + [1 5], v(i,1)*[1 1], 'k-');
  text(v(i,2) + w + 6, v(i,1), sprintf('dia %.2f @ z=%.1f', 2*v(i,2), v(i,1)), 'FontSize', 8);
end
axis equal; set(gca, 'YDir', 'reverse'); xlim([-32 48]);
xlabel('r  [mm]'); ylabel('z below inlet  [mm]'); title(ttl);
end
