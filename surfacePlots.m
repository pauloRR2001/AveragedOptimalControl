%% Figure 1: full trajectory as thrust-colored surface + target/initial curves

R = recon.r_cart(:,1:3);
Tsurf = recon.T(:);

% Make sure thrust vector matches reconstructed Cartesian samples
if length(Tsurf) ~= size(R,1)
    tau_R = recon.tau(:);
    tau_T = linspace(recon.tau(1), recon.tau(end), length(Tsurf)).';
    Tsurf = interp1(tau_T, Tsurf, tau_R, 'linear', 'extrap');
end

% Tube settings
tube_radius = 0.03;     % [DU], tune for visual thickness
nCirc = 24;

[Xs,Ys,Zs,Cs] = curve_tube_surface_colored(R, Tsurf, tube_radius, nCirc);

figure;
hold on;

% Controlled transfer as surface colored by reconstructed thrust magnitude
surf(Xs, Ys, Zs, Cs, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.92, ...
    'DisplayName', 'Controlled transfer');

% Target orbit
plot3(out_f.cart(:,1), out_f.cart(:,2), out_f.cart(:,3), ...
    'r', ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Target orbit');

% Initial orbit
plot3(out_0.cart(:,1), out_0.cart(:,2), out_0.cart(:,3)+0.01, ...
    'm', ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Initial Orbit');

earthy(1, 'Earth', 1,[0,0,0])

legend('Controlled Transfer','Target Orbit','Initial Orbit','','Location','best')

grid on; box on; axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
title('Transfer Trajectory Colored by Reconstructed Thrust Magnitude');

colormap(turbo);
cb = colorbar;
ylabel(cb,'T');

caxis([params.T_min params.T_max]);

view(3);
camlight headlight;
lighting gouraud;

function [Xs,Ys,Zs,Cs] = curve_tube_surface_colored(R, Cline, radius, nCirc)
% Build a tube surface around a 3D curve R [N x 3].
% Color field Cline [N x 1] is replicated around each tube cross-section.

    R = R(:,1:3);
    Cline = Cline(:);

    N = size(R,1);

    if length(Cline) ~= N
        error('Cline must have the same number of samples as R.');
    end

    % Remove any duplicate or invalid trajectory points
    good = all(isfinite(R),2) & isfinite(Cline);
    R = R(good,:);
    Cline = Cline(good);

    N = size(R,1);

    % Tangent vectors
    That = zeros(N,3);
    That(2:N-1,:) = R(3:N,:) - R(1:N-2,:);
    That(1,:)     = R(2,:)   - R(1,:);
    That(N,:)     = R(N,:)   - R(N-1,:);

    tnorm = vecnorm(That,2,2);
    tnorm(tnorm < 1e-14) = 1;
    That = That ./ tnorm;

    % Allocate normal/binormal frame
    Nvec = zeros(N,3);
    Bvec = zeros(N,3);

    % Initial reference vector
    ref = [0 0 1];

    if abs(dot(ref,That(1,:))) > 0.9
        ref = [0 1 0];
    end

    Nvec(1,:) = cross(That(1,:),ref);

    if norm(Nvec(1,:)) < 1e-14
        ref = [1 0 0];
        Nvec(1,:) = cross(That(1,:),ref);
    end

    Nvec(1,:) = Nvec(1,:)/norm(Nvec(1,:));
    Bvec(1,:) = cross(That(1,:),Nvec(1,:));
    Bvec(1,:) = Bvec(1,:)/norm(Bvec(1,:));

    % Propagate frame along curve
    for i = 2:N

        Nvec(i,:) = cross(Bvec(i-1,:),That(i,:));

        if norm(Nvec(i,:)) < 1e-14
            Nvec(i,:) = Nvec(i-1,:);
        else
            Nvec(i,:) = Nvec(i,:)/norm(Nvec(i,:));
        end

        Bvec(i,:) = cross(That(i,:),Nvec(i,:));

        if norm(Bvec(i,:)) < 1e-14
            Bvec(i,:) = Bvec(i-1,:);
        else
            Bvec(i,:) = Bvec(i,:)/norm(Bvec(i,:));
        end

    end

    theta = linspace(0,2*pi,nCirc);

    Xs = zeros(N,nCirc);
    Ys = zeros(N,nCirc);
    Zs = zeros(N,nCirc);
    Cs = zeros(N,nCirc);

    for i = 1:N
        for j = 1:nCirc

            offset = radius*cos(theta(j))*Nvec(i,:) + ...
                     radius*sin(theta(j))*Bvec(i,:);

            P = R(i,:) + offset;

            Xs(i,j) = P(1);
            Ys(i,j) = P(2);
            Zs(i,j) = P(3);

            Cs(i,j) = Cline(i);

        end
    end

end

% Figure 4: costates vs full transfer time with terminal target values
% Use reconstructed costate history if available
if isfield(recon,'lambda')
    lambda_plot = recon.lambda;
    t_costate_days = recon.tau(:) * alpha0 * TU / 86400;
else
    lambda_plot = lambda;
    t_costate_days = tau(:) * alpha0 * TU / 86400;
end

figure;

costate_idx = [1 2 3 4 5 6 9];
costate_labels = { ...
    '$\lambda_p$', ...
    '$\lambda_f$', ...
    '$\lambda_g$', ...
    '$\lambda_h$', ...
    '$\lambda_k$', ...
    '$\lambda_L$', ...
    '$\lambda_m$'};

for j = 1:numel(costate_idx)

    subplot(4,2,j);

    idx = costate_idx(j);
    plot(t_costate_days, lambda_plot(:,idx), 'b', 'LineWidth', 1.5); 
    hold on;

    % Terminal transversality targets
    if idx == 6 || idx == 9
        yline(0, 'r--', 'LineWidth', 1.5);
        legend(costate_labels{j}, 'target', ...
               'Location', 'best', 'Interpreter', 'latex');
    else
        legend(costate_labels{j}, ...
               'Location', 'best', 'Interpreter', 'latex');
    end

    grid on; box on;
    xlabel('time [days]', 'Interpreter', 'latex');
    ylabel(costate_labels{j}, 'Interpreter', 'latex');
    title([costate_labels{j} ' vs time'], 'Interpreter', 'latex');

end

sgtitle('Costate History with Terminal Transversality Targets', ...
        'Interpreter', 'latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot 2: full reconstructed controlled trajectories for all satellites
%         as colored tube surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;

tube_radius = 0.02;   % [DU], adjust for appearance
nCirc = 20;           % points around tube cross-section

clr = lines(N_s);     % one color per satellite

for sat = 1:N_s

    recon = recon_all{sat};
    R = recon.r_cart(:,1:3);

    [Xs,Ys,Zs] = curve_tube_surface(R, tube_radius, nCirc);

    surf(Xs, Ys, Zs, ...
        'FaceColor', clr(sat,:), ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.95, ...
        'DisplayName', sprintf('Sat %02d, RAAN %.1f deg', ...
            sat, target_all(sat).raan_deg));

end

earthy(1, 'Earth', 1,[0,0,0]);

% Prevent Earth object from appearing as an empty legend entry
hEarth = findobj(gcf,'Type','patch','-or','Type','surface');
set(hEarth(1),'HandleVisibility','off');

grid on; box on; axis equal;
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
zlabel('z', 'Interpreter', 'latex');
title('Full Reconstructed Controlled Trajectory Surfaces for All Satellites', 'Interpreter', 'latex');
legend('Location','bestoutside', 'Interpreter', 'latex');

view(3);
camlight headlight;
lighting gouraud;

function [Xs,Ys,Zs] = curve_tube_surface(R, radius, nCirc)
% Build a tube surface around a 3D curve R [N x 3]

    R = R(:,1:3);
    N = size(R,1);

    % Remove invalid points
    good = all(isfinite(R),2);
    R = R(good,:);
    N = size(R,1);

    if N < 2
        error('Trajectory must contain at least 2 valid points.');
    end

    % Tangent vectors
    That = zeros(N,3);
    That(2:N-1,:) = R(3:N,:) - R(1:N-2,:);
    That(1,:)     = R(2,:)   - R(1,:);
    That(N,:)     = R(N,:)   - R(N-1,:);

    tnorm = vecnorm(That,2,2);
    tnorm(tnorm < 1e-14) = 1;
    That = That ./ tnorm;

    % Normal/binormal frames
    Nvec = zeros(N,3);
    Bvec = zeros(N,3);

    % Initial reference direction
    ref = [0 0 1];
    if abs(dot(ref,That(1,:))) > 0.9
        ref = [0 1 0];
    end

    Nvec(1,:) = cross(That(1,:),ref);
    if norm(Nvec(1,:)) < 1e-14
        ref = [1 0 0];
        Nvec(1,:) = cross(That(1,:),ref);
    end
    Nvec(1,:) = Nvec(1,:)/norm(Nvec(1,:));

    Bvec(1,:) = cross(That(1,:),Nvec(1,:));
    Bvec(1,:) = Bvec(1,:)/norm(Bvec(1,:));

    % Propagate frame
    for i = 2:N

        Nvec(i,:) = cross(Bvec(i-1,:),That(i,:));

        if norm(Nvec(i,:)) < 1e-14
            Nvec(i,:) = Nvec(i-1,:);
        else
            Nvec(i,:) = Nvec(i,:)/norm(Nvec(i,:));
        end

        Bvec(i,:) = cross(That(i,:),Nvec(i,:));

        if norm(Bvec(i,:)) < 1e-14
            Bvec(i,:) = Bvec(i-1,:);
        else
            Bvec(i,:) = Bvec(i,:)/norm(Bvec(i,:));
        end

    end

    theta = linspace(0,2*pi,nCirc);

    Xs = zeros(N,nCirc);
    Ys = zeros(N,nCirc);
    Zs = zeros(N,nCirc);

    for i = 1:N
        for j = 1:nCirc

            offset = radius*cos(theta(j))*Nvec(i,:) + ...
                     radius*sin(theta(j))*Bvec(i,:);

            P = R(i,:) + offset;

            Xs(i,j) = P(1);
            Ys(i,j) = P(2);
            Zs(i,j) = P(3);

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FULL TRAJECTORY SURFACES OF ALL SATELLITES
% Break trajectory into continuous chunks before tubing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;

tube_radius = 0.02;      % [DU]
nCirc = 20;
clr = lines(N_s);

jump_factor = 8;         % increase if too many breaks, decrease if bridges remain

for sat = 1:N_s

    r_transfer = recon_transfer_all{sat}.r_cart(:,1:3);
    r_circ     = recon_circ_all{sat}.r_cart(:,1:3);

    % Keep legs separate to avoid artificial bridge between transfer/circ
    legs = {r_transfer, r_circ};

    for ileg = 1:numel(legs)

        R = legs{ileg};
        R = R(all(isfinite(R),2),:);

        if size(R,1) < 3
            continue;
        end

        % Distance between consecutive reconstructed points
        dR = vecnorm(diff(R,1,1),2,2);

        % Robust threshold for detecting artificial jumps
        dR_med = median(dR(dR > 0));

        if isempty(dR_med) || ~isfinite(dR_med)
            continue;
        end

        jump_idx = find(dR > jump_factor*dR_med);

        % Segment boundaries
        i_start = [1; jump_idx + 1];
        i_end   = [jump_idx; size(R,1)];

        for iseg = 1:length(i_start)

            Rseg = R(i_start(iseg):i_end(iseg),:);

            if size(Rseg,1) < 5
                continue;
            end

            [Xs,Ys,Zs] = curve_tube_surface(Rseg, tube_radius, nCirc);

            if ileg == 1 && iseg == 1
                name_i = sprintf('Sat %02d', sat);
                hv = 'on';
            else
                name_i = '';
                hv = 'off';
            end

            surf(Xs,Ys,Zs, ...
                'FaceColor', clr(sat,:), ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.90, ...
                'DisplayName', name_i, ...
                'HandleVisibility', hv);

        end

    end

end

earthy(1, 'Earth', 1, [0,0,0]);

hEarth = findobj(gcf,'Type','patch','-or','Type','surface');
set(hEarth(1),'HandleVisibility','off');

grid on; box on; axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
title('Full Trajectory Surfaces of All Satellites');
legend('Location','bestoutside');
view(3);
camlight headlight;
lighting gouraud;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ground tracks for all satellites from rotating-frame final orbits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;

vernalEq = 30;
clr = lines(N_s);

for sat = 1:N_s

    r_rot = orbit_rot_all{sat}.r_rot;

    plotGroundTrackFromRotatingState( ...
        r_rot, ...
        vernalEq, ...
        'NewFigure', false, ...
        'Color', clr(sat,:), ...
        'LineWidth', 1.2, ...
        'DisplayName', sprintf('Sat %02d, RAAN %.1f deg', ...
            sat, target_circ_all(sat).raan_deg), ...
        'Title', 'Ground Tracks of Final Circularized Orbits');

end

legend('Location','bestoutside');