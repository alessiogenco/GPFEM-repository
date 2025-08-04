close all
clear all

%% ...........................Determination wedge failure geometry .......................................
%% Properties anchor taper 
D_outer = 0.30; % [m]
D_inner = 0.25; % [m]
H_taper = 0.35; % [m]

%% Properties rock (mudstone from Dorth)

gamma_dry = 25.7; % [kN/m3]
E_rock = 32.5; % [GPa]
nu = 0.3; % [-]

% HB parameters

GSI = 100; % [-] Intact rock
% GSI = 50; % [-] Reduced properties

mi = 9; % [-]
D = 0; %[-]

Sigma_ci = 6300; % [kPa]

UCS = Sigma_ci/1000; % [MPa]

% Parameters calculated

mb = mi*exp((GSI-100)/(28-14*D)); % [-]

s = exp((GSI-100)/(9-3*D)); % [-]

a = 1/2+1/6*(exp(-GSI/15)-exp(-20/3)); % [-]

Sigmac = abs(Sigma_ci)*s^(a); % [kPa]

Sigmat = (s*abs(Sigma_ci))/mb; % [kPa]

%% Parameters for the taper surface
r1 = D_inner/2;      % Bottom radius
r2 = D_outer/2;      % Top radius
h = H_taper;         % Height
n = 50;              % Resolution (number of points around circumference)

% Create theta and height vectors
theta = linspace(0, pi/2, n);
theta_tot = linspace(0, pi, n);
z = linspace(0, h, 2); % Two layers: bottom and top

% Create meshgrid
[Theta, Z] = meshgrid(theta, z);
[Theta2, Z] = meshgrid(theta_tot, z);

% Radius varies linearly with height
R = r1 + (r2 - r1) * (Z / h);

% Parametric equations for the taper
X = R .* cos(Theta);
Y = R .* sin(Theta);

X_tot = R .* cos(Theta2);
Y_tot = R .* sin(Theta2);

% Plot the taper (cyan surface)
figure
surf(X, Y, Z, 'FaceAlpha', 0.5, 'FaceColor', 'cyan', 'EdgeColor', 'none')
hold on
axis equal
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
title('Shallow wedge of failure from field tests')
grid on
view(3)

%% Plot side triangle (red patch)

% Y_max = 0.388; % [m]
% 
% x_tri = [0, 0, 0];                          % Constant X
% y_tri = [D_inner/2, D_outer/2, Y_max];      % Y coordinates
% z_tri = [0, H_taper, H_taper];               % Z coordinates
% fill3(x_tri, y_tri, z_tri, 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k')
% hold on

%% Build the top surface perimeter

Y_max = 0.358; % [m]

% ① Line along top of red triangle
x_top_red = [0, 0];
y_top_red = [D_outer/2, Y_max];
z_top_red = [H_taper, H_taper];

% ② Top surface: Circular arc (connect red to cyan smoothly)
circle_diameter = 1.55 * D_outer;
circle_radius = circle_diameter / 2;
theta_circle = linspace(-pi/2 + 0.21*pi, pi/2, 200);  % Arc angles

% Full circle coordinates
x_circle = circle_radius * cos(theta_circle);
y_circle = circle_radius * sin(theta_circle) + ((D_outer/2) + 0.21215*D_outer);
z_circle = ones(size(x_circle)) * H_taper;  % Constant Z

% Find indexes on circular arc
start_circle = [0, 0.4551];              % From top of red patch
end_circle = [0.148017, 0.0228226];       % To top of taper patch

[~, idx_start_circle] = min(vecnorm([x_circle; y_circle] - start_circle', 2, 1));
[~, idx_end_circle] = min(vecnorm([x_circle; y_circle] - end_circle', 2, 1));

% Extract arc portion
if idx_start_circle < idx_end_circle
    x_arc = x_circle(idx_start_circle:idx_end_circle);
    y_arc = y_circle(idx_start_circle:idx_end_circle);
else
    x_arc = x_circle(idx_end_circle:idx_start_circle);
    y_arc = y_circle(idx_end_circle:idx_start_circle);
end
z_arc = H_taper * ones(size(x_arc));

% ③ Line along top of cyan taper surface
X_top = X(2,:);   % Top row at Z = H_taper
Y_top = Y(2,:);   % Top row at Z = H_taper

% Find matching index on (X_top,Y_top) from end of arc
start_taper = [0.148017, 0.0228226];
[~, idx_start_taper] = min(vecnorm([X_top; Y_top] - start_taper', 2, 1));

% Extract taper part from arc end to taper end
x_taper = X_top(idx_start_taper:end);
y_taper = Y_top(idx_start_taper:end);
z_taper = H_taper * ones(size(x_taper));

x_top_surf = [x_arc 0 flip(x_taper)];
y_top_surf = [y_arc 0.15 flip(y_taper)];
z_top_surf = H_taper * ones(size(y_top_surf));

% %% ④ Assemble the full perimeter: red line -> arc -> taper
% x_fill = [x_top_red(1), x_top_red(2), x_arc, x_taper];
% y_fill = [y_top_red(1), y_top_red(2), y_arc, y_taper];
% z_fill = H_taper * ones(size(x_fill));

% Plot the top surface (blue patch)
% patch(x_fill, y_fill, z_fill, 'b', 'FaceAlpha', 0.6)
patch(x_top_surf, y_top_surf, z_top_surf, 'b', 'FaceAlpha', 0.6)


%%

x_back_surf = [x_circle X(1,end) flip(X(1,6:end)) x_circle(1)];
y_back_surf = [y_circle Y(1,end) flip(Y(1,6:end)) y_circle(1)];
z_back_surf = [H_taper * ones(size(x_circle)) 0 zeros(size(X(1,6:end))) 0.35];

patch(x_back_surf, y_back_surf, z_back_surf, 'g', 'FaceAlpha', 0.6,'EdgeColor', 'k')

% figure
% plot3(x_back_surf,y_back_surf,z_back_surf);

%% 

%% ----------------------- Simpified wedge geometry ----------------------------------

figure(2)

% Plot the taper (cyan surface)
figure
surf(X_tot, Y_tot, Z, 'FaceAlpha', 0.5, 'FaceColor', 'red', 'EdgeColor', 'none')
hold on
axis equal
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
% title('Shallow wedge simplified')
grid on
view(3)
% fill3(x_tri, y_tri, z_tri, 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k')
% hold on

% Plot top surface 

% x_top_simpl = [0 0 circle_diameter/2 (X(2,1:end))];
x_top_simpl = [min(X_tot(2,:)) -circle_diameter/2 circle_diameter/2 X_tot(2,:)];
y_top_simpl = [0 Y_max Y_max Y_tot(2,:)];
% y_top_simpl = [0.15 0.4551 0.4551 (Y(2,1:end))];
z_top_simp = H_taper * ones(size(x_top_simpl)); 
% 
fill3(x_top_simpl, y_top_simpl, z_top_simp, 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k')
hold on

% Plot back surfaces

x_sid1_simpl = [min(X_tot(2,:)) -circle_diameter/2  min(X_tot(1,:))];
y_sid1_simpl = [0 Y_max 0];
z_sid1_simpl = [H_taper H_taper 0];

x_sid2_simpl = [max(X_tot(2,:)) circle_diameter/2 max(X_tot(1,:))];
y_sid2_simpl = [0 Y_max 0];
z_sid2_simpl = [H_taper H_taper 0];

fill3(x_sid1_simpl, y_sid1_simpl, z_sid1_simpl, 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k')

hold on

fill3(x_sid2_simpl, y_sid2_simpl, z_sid2_simpl, 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k')

x_bac_simp = [circle_diameter/2 X_tot(1,:) -circle_diameter/2 circle_diameter/2];
y_bac_simp = [Y_max Y_tot(1,:) Y_max Y_max];
z_bac_simp = [H_taper 0 * ones(size(X_tot(1,:))) H_taper H_taper];
% 
fill3(x_bac_simp, y_bac_simp, z_bac_simp, 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k')

%% Calculation areas and volume of the wedge simplified 

V_taper = (1/3)*pi*H_taper*((D_outer/2)^2 + (D_inner/2)^2 + ((D_outer/2)*(D_inner/2))); % [m3]

theta = atand(((circle_diameter/2)-(max(X_tot(2,:))))/(Y_max)); % [deg]

beta = atand(Y_max/H_taper);  % [deg]

area_section_triangle = (H_taper * Y_max)/2; % [m2]

area_sides = 0.5*(H_taper^2)*tand(beta)*secd(theta); % [m2]

% area_back = (D_inner*H_taper*secd(beta)) + ((H_taper*tand(beta)*tand(theta)) + 0.5*D_outer - 0.5*D_inner)*H_taper*secd(beta);

% area_back = H_taper*secd(beta) * (D_inner + ((H_taper*tand(beta)*tand(theta)) + 0.5*D_outer - 0.5*D_inner)); % [m2] - Calculated as sum of rectangle are + 2 times triangle ares
% area_baack = 0.5 * H_taper*secd(beta) * (D_outer + 2* H_taper*tand(beta)*tan(theta) + D_inner); 
area_back = (((D_outer + D_inner)/2) + H_taper*tand(beta)*tand(theta))*H_taper*secd(beta);

%V_wedge = (area_section_triangle*(circle_diameter+max(Y(2,:))*2)/2) - (V_taper/2); % [m3]

Volume_wedge = area_section_triangle * 0.5*((2*(H_taper*tand(beta)*tand(theta))+D_outer) + D_outer) - (V_taper/2); % [m3]

fprintf('Volume wedge simplified = %.5f m^3\n', Volume_wedge)

W = gamma_dry*Volume_wedge; % [kN]

%% Constitutive modelling formulations

% Liang et al. (2009)

% MC parameters derivation from Yang (2006) -> Hoek et al. (2002)

Sigma_3 = gamma_dry * H_taper/3;

Sigma_1 = Sigma_3 + Sigma_ci*(mb*(Sigma_3/Sigma_ci)+s)^a; % [kPa]

Sigma_n = Sigma_3 + (((Sigma_1-Sigma_3)^2)/(2*(Sigma_1-Sigma_3)+0.5*mb*Sigma_ci)); % [kPa]
tau = (Sigma_n-Sigma_3)*sqrt(1 + (mb*Sigma_ci)/(2*(Sigma_1-Sigma_3))); % [kPa]

phi = 90 - asind((2*tau)/(Sigma_1-Sigma_3)); % [deg]
c = tau - Sigma_n * tand(phi); % [kPa]

% Phi and c derivation from my HB fitting 

phi_fittinghb = 38; % [deg]
c_fittinghb = 150; % [kPa]

% phi = phi_fittinghb;
% c = c_fittinghb;

% Carranza-Torres (2004) MC parameters derivations from Elsokary et al.
% (2023)

Sigma_v = gamma_dry*H_taper; % [kPa] Minor principal stress

Sigma_3rect = Sigma_v + (gamma_dry*H_taper/3); % [kPa]
Sigma_3tria = Sigma_v + (gamma_dry*H_taper/2); % [kPa]


Sigma_3c = Sigma_3rect; % Sigma_3c2 > Sigma_V -> Major principal stress (Wedge under Passive stress)
%Sigma_3c = Sigma_3tria; % Sigma_3c2 > Sigma_V -> Major principal stress (Wedge under Passive stress)

phic = asind(1-(2/(2+a*mb*(mb*(Sigma_3c/Sigma_ci)+s)^(a-1)))); % [deg]

K_0 = (1-sind(phi)); % [-]

Sigma_3c2 = K_0*(Sigma_v + (gamma_dry*H_taper/3)); % [kPa] Sigma_3c2 < Sigma_V

Sigma_nc = Sigma_3c2 + (Sigma_ci/2)*((mb*(Sigma_3c2/Sigma_ci)+s)^a)*(1 - ((a*mb*((mb*(Sigma_3c2/Sigma_ci)+s)^(a-1)))/(2+a*mb*(mb*(Sigma_3c2/Sigma_ci)+s)^(a-1)))); % [kPa]
tau_c = Sigma_ci*((mb*(Sigma_3c2/Sigma_ci)+s)^a)*((sqrt(1+a*mb*(mb*(Sigma_3c2/Sigma_ci)+s)^(a-1)))/(2+a*mb*(mb*(Sigma_3c2/Sigma_ci)+s)^(a-1))); % [kPa]


%% ------------------------- Determination of the forces acting on the system --------------------------------------

% Side surfaces forces 

F_ns = area_sides*Sigma_n; % [kN]
F_ss = area_sides*tau; % [kN]

% Back surface forces 

F_nb = area_back*Sigma_n; % [kN]
F_sb = area_back*tau; % [kN]

% Total resistance in horizontal lateral direction

P_t = 2*F_ss*cosd(theta)*sind(beta) + F_sb*sind(beta) + F_nb * cosd(beta) - 2*F_ns*sind(theta); % [kN]


fprintf('P_t = %.5f kN', P_t)
%% Derivation for determination of P_ult as ultimate rock resistance [kN/m]

% Known/constant value substitution BEFORE declaring symbolic H_taper:
H_taper_val = 0.3;

% Start symbolic work
syms H_taper 

% Define expressions
N_1 = (mb * ((gamma_dry * H_taper)/(3 * Sigma_ci)) + s);
N_2 = 0.5 * mb * Sigma_ci;

% τ(H_taper)
tau_expr = ((Sigma_ci^2) * N_1^(2*a)) / (2 * Sigma_ci * N_1^a + N_2) * sqrt(1 + (mb * Sigma_ci) / (2 * Sigma_ci * N_1^a));
d_F_ss = (((Sigma_ci^2) * N_1^(2*a)) / (2 * Sigma_ci * N_1^a + N_2) * sqrt(1 + (mb * Sigma_ci) / (2 * Sigma_ci * N_1^a)))*(0.5*(H_taper^2)*tand(beta)*secd(theta));


% Derivative of τ w.r.t. H_taper

d_F_ss_dH = diff(d_F_ss, H_taper);
d_F_ss_dH = double(subs(d_F_ss_dH, H_taper, H_taper_val));

%% Normal stress σ_n derivation and force
% σ_n(H_taper)
sigma_n_expr = (gamma_dry * H_taper / 3) + (Sigma_ci * N_1^(2*a)) / (2 * Sigma_ci * N_1^a + N_2);
d_F_ns = ((gamma_dry * H_taper / 3) + (Sigma_ci * N_1^(2*a)) / (2 * Sigma_ci * N_1^a + N_2))*(0.5*(H_taper^2)*tand(beta)*secd(theta));

% Derivative

d_F_ns_dH = diff(d_F_ns, H_taper);
d_F_ns_dH = double(subs(d_F_ns_dH, H_taper, H_taper_val));

%% Calculation derivatives F_sb and F_nb related to the back surface

d_F_sb = ((Sigma_ci^2) * N_1^(2*a)) / (2 * Sigma_ci * N_1^a + N_2) * sqrt(1 + (mb * Sigma_ci) / (2 * Sigma_ci * N_1^a))*(H_taper*secd(beta) * (D_inner + ((H_taper*tand(beta)*tand(theta)) + 0.5*D_outer - 0.5*D_inner)));

% Derivative
d_F_sb_dH = diff(d_F_sb, H_taper);
d_F_sb_dH = double(subs(d_F_sb_dH, H_taper, H_taper_val));

%

d_F_nb = ((gamma_dry * H_taper / 3) + (Sigma_ci * N_1^(2*a)) / (2 * Sigma_ci * N_1^a + N_2))*(H_taper*secd(beta) * (D_inner + ((H_taper*tand(beta)*tand(theta)) + 0.5*D_outer - 0.5*D_inner)));

% Derivative
d_F_nb_dH = diff(d_F_nb, H_taper);
d_F_nb_dH = double(subs(d_F_nb_dH, H_taper, H_taper_val));


%% ----------------------- P_ult determination  ---------------------

P_ult = 2*d_F_ss_dH*cosd(theta)*sind(beta) + d_F_sb_dH*sind(beta) + d_F_nb_dH*cosd(beta) - 2*d_F_ns_dH*sind(theta); 


fprintf('P_ult = %.5f kN/m', P_ult)

%% ------------------------- p-y curve deriation -----------------------------

% K_i determination 

K_i = 48.719*(UCS^(-0.5088))*(E_rock^(0.5378)); 

% Shape of the curve 

y_D = linspace(0.00001,0.1,100); % [m]

p_p_ult = (0.7787*exp(2.501.*y_D)) - (0.8258*exp(-49.53.*y_D));

p = p_p_ult*P_ult; % [kN/m]




%% Plotting the p-y curve normalized to import into Benjamin s code from numerical data GPFEM 2D

figure12=figure

fig_settings_print4

subplot1 = subplot(1,1,1,'Parent',figure12);

box(subplot1,'on');
hold(subplot1,'all');

set(subplot1,'XMinorTick','on');

subplot1 (1,1,1)



y_d_vec=y_D(1:round(length(y_D)/5*4/20):round(length(y_D)/5*4)); % x vector discretized
y_d_vec(1)=0;
y_d_vec(end)=y_D(end);

p_2D=p(1:round(length(y_D)/5*4/20):round(length(y_D)/5*4)); % y vector discretized
p_2D(end)= p(end);
p_2D(1)=0;
 
y = y_d_vec*0.3;

plot(y_D(1:end)*0.3,p ...
    (1:end),'Parent',subplot1,'MarkerFaceColor','k','Marker','none','MarkerSize',3,'LineStyle','-','Color','k','LineWidth',2,...
        'DisplayName',('$$P-y$$ $$LEM$$' ));
hold on
s=plot(y_d_vec*0.3,p_2D ...
    (1:end),'Parent',subplot1,'MarkerFaceColor','none','Marker','o','MarkerSize',7,'LineStyle','-','Color','b','LineWidth',2,...
        'DisplayName',('$$P-y$$ $$interp$$' ));

% xVectors = get(s,'XData');
% yVectors = get(s,'YData');
% 
% xVector_s = xVectors.';
% yVector_s = yVectors.';
% 
% Py_Vectors = [xVector_s(1:end) yVector_s(1:end)];

 set(gca, 'FontSize', fsz)
   set(gcf,'color','w');
   axis square
   box off
   grid off
   
ylabel('$$p$$ $$[kN/m]$$','Interpreter','latex');
xlabel('$$Y$$ $$[m]$$','Interpreter','latex');

box off
grid off
axis square
% axis equal
h=legend(subplot1,'show','Location', 'northwest');

 set(h,'FontSize', fsz,'Interpreter','latex','Location', 'northwest')
 legend('Boxoff')
%    ylim([0 20]);
%xlim([0 0.5]);
% if print_yn==1
% print('p-y depth GPFEM','-dpng','-r300')
% end

%% Plotting parametric investigation UCS, GSI, and H_tap

figure4=figure

subplot1 = subplot(1,1,1,'Parent',figure4);

box(subplot1,'on');
hold(subplot1,'all');

set(subplot1,'XMinorTick','on');

subplot1 (1,1,1)


plot([3 6 9 12 15],[281.9 561.85 841.78 1121.72 1401.65],'Parent',subplot1,'MarkerFaceColor','none','Marker','square','MarkerSize',6,'LineStyle','-','Color','k','LineWidth',1,...
        'DisplayName',['$$LEM P_u$']);

 set(gca, 'FontSize', fsz)
   set(gcf,'color','w');
   axis square
   box on
   grid off
   
ylabel('$$P_u$$ $$[kN]$$','Interpreter','latex');
xlabel('$$UCS$$ $$[MPa]$$','Interpreter','latex');

box on
grid off
axis square
h=legend(subplot1,'show','Location', 'northwest');

 set(h,'FontSize', fsz,'Interpreter','latex','Location', 'northwest')
 legend('Boxoff')
  ylim([0 1500]);
  xlim([0 18]);
% if print_yn==1
% print('Pulloutcomp','-dpng','-r300')
% end

figure5=figure

subplot1 = subplot(1,1,1,'Parent',figure5);

box(subplot1,'on');
hold(subplot1,'all');

set(subplot1,'XMinorTick','on');

subplot1 (1,1,1)


plot([100 80 60 40 20],[589.84 169.17 48.92 15.04 4.78],'Parent',subplot1,'MarkerFaceColor','none','Marker','square','MarkerSize',6,'LineStyle','-','Color','r','LineWidth',1,...
        'DisplayName',['$$LEM P_u$']);

 set(gca, 'FontSize', fsz)
   set(gcf,'color','w');
   axis square
   box on
   grid off
   
ylabel('$$P_u$$ $$[kN]$$','Interpreter','latex');
xlabel('$$GSI$$ $$[-]$$','Interpreter','latex');

box on
grid off
axis square
h=legend(subplot1,'show','Location', 'northwest');

 set(h,'FontSize', fsz,'Interpreter','latex','Location', 'northwest')
 legend('Boxoff')
  ylim([0 800]);
  xlim([0 100]);
% if print_yn==1
% print('Pulloutcomp','-dpng','-r300')
% end

figure6=figure

subplot1 = subplot(1,1,1,'Parent',figure6);

box(subplot1,'on');
hold(subplot1,'all');

set(subplot1,'XMinorTick','on');

subplot1 (1,1,1)


plot([0.15 0.25 0.35 0.45 0.55],[444.06 525.56 589.84 641.70 685.52],'Parent',subplot1,'MarkerFaceColor','none','Marker','square','MarkerSize',6,'LineStyle','-','Color','b','LineWidth',1,...
        'DisplayName',['$$LEM P_u$']);

 set(gca, 'FontSize', fsz)
   set(gcf,'color','w');
   axis square
   box on
   grid off
   
ylabel('$$P_u$$ $$[kN]$$','Interpreter','latex');
xlabel('$$H_{tap}$$ $$[m]$$','Interpreter','latex');

box on
grid off
axis square
h=legend(subplot1,'show','Location', 'northwest');

 set(h,'FontSize', fsz,'Interpreter','latex','Location', 'northwest')
 legend('Boxoff')
  ylim([0 800]);
  xlim([0 0.6]);
% if print_yn==1
% print('Pulloutcomp','-dpng','-r300')
% end

