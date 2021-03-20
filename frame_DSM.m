clear; clc; close all;
% unit: N, m

nodes = [0 0;
        0 6;
        0 12;
        18 12;
        18 6;
        18 0];

eles = [1 2;
       2 3;
       3 4;
       4 5;
       5 6;
       2 5
       ];


resolution = 3;  % larger, finer

% Geometry properties
mag_factor = 5000;  % magnification factor for visualization
B = .50;
H_c = .50;
H_b = 1.;
A_c = B*H_c;
I_c = 1./12.*B*H_c^3;
A_b = B*H_b;
I_b = 1./12.*B*H_b^3;

elenum = size(eles, 1);
E = 30e9 * ones(elenum, 1);
A = [A_c A_c A_b A_c A_c A_b]';
I = [I_c I_c I_b I_c I_c A_b]';

% 3 means xyz constrained, 2 means xy constrained
% 13 means z constrained, 12 means y constrained, 11 means x constrained,
% hinge not contained
Supporting = [3 0 0 0 0 3 ]';

% nodal forces, moments (global), [Fx Fy node#]
ExF = [3e3 0 2;
       -3e3 0 3]; 

% element forces
% element dist load: [mag(include local dir) eletag]
w = [-1e3 3;
    -3e3 6;
    ]; 
% element force: [mag(include local dir) eletag DistFromStartNode]
eleF = [];  % waiting for update
%%

% length container
L = [];
for i = 1:size(eles, 1)
    start_node = eles(i, 1);
    end_node = eles(i, 2);
    vec = nodes(end_node, :) - nodes(start_node, :);
    L(i) = norm(vec);
end
L = L';


%% assemble new "nodes"
num_nodes_inserted_set = floor(L * resolution);  % dim: num_ele by 1

nodes_old = nodes;
num_nodes_inserted_cm_set = [0 sum(num_nodes_inserted_set*...
                    ones(size(num_nodes_inserted_set')).*...
                    triu(ones(numel(num_nodes_inserted_set)), 0), 1)]';  % culumative, [1,2,3]->[1,3,6] 

nodes = zeros(size(nodes_old, 1) + sum(num_nodes_inserted_set), 2);  % nodes(0, :) = [0 0] should hold
%nodes(oldnodes_ind, :) = nodes_old;

oldnodes_ind = [1];
nodes_eles_table = zeros(size(nodes, 1), 1);  % only for inserted nodes
for i=1:size(eles, 1)
  start_node = nodes_old(eles(i, 1), :);  %coordinate
  end_node = nodes_old(eles(i, 2), :);  %coordinate

  x_set = linspace( start_node(1), end_node(1), num_nodes_inserted_set(i)+2 )';
  x_set = x_set(2:size(x_set, 1)-1);
  y_set = linspace( start_node(2), end_node(2), num_nodes_inserted_set(i)+2 )';
  y_set = y_set(2:size(y_set, 1)-1);

  nodes(i+1+num_nodes_inserted_cm_set(i) : i+num_nodes_inserted_cm_set(i+1), 1) = x_set;
  nodes(i+1+num_nodes_inserted_cm_set(i) : i+num_nodes_inserted_cm_set(i+1), 2) = y_set;
  
  nodes_eles_table(i+1+num_nodes_inserted_cm_set(i) : i+num_nodes_inserted_cm_set(i+1)) = i;

  if ismember(end_node, nodes, 'rows')
    % do nothing
  else
    nodes(i+1+num_nodes_inserted_cm_set(i+1), 1) = end_node(1);
    nodes(i+1+num_nodes_inserted_cm_set(i+1), 2) = end_node(2);

    oldnodes_ind = [oldnodes_ind; i+1+num_nodes_inserted_cm_set(i+1)];
  end

end

%% assemble new "eles"
eles_old = eles;

for i=1:numel(eles)
  eles(i) = oldnodes_ind(eles(i));
end

eles_org = eles;  % eles tags using new nodes coordinates

eles = [];
A_old = A;
I_old = I;
E_old = E;
A = []; E = []; I = [];

orgeletag_set = [];
for i=1:numel(num_nodes_inserted_set)
    orgeletag_set = [orgeletag_set; i*ones(num_nodes_inserted_set(i)+1, 1)];
end
orgeletag_set = [orgeletag_set (1:1:numel(orgeletag_set))'];

for i=1:size(eles_org, 1)
    %eles = [eles; eles_org(i, :)];
    if i==6
        a=1;
    end

    num_eles = sum(orgeletag_set(:, 1)==i);
    start_node = eles_org(i, 1);
    end_node = eles_org(i, 2);
    inserted_nodes = find(nodes_eles_table==i)';
    ele_table = [start_node inserted_nodes end_node];

    from_nodes = ele_table(1: numel(ele_table)-1);
    to_nodes = ele_table(2: numel(ele_table));

    eles = [eles; [from_nodes' to_nodes']];

    A = [A; A_old(i)*ones(num_eles, 1)];
    E = [E; E_old(i)*ones(num_eles, 1)];
    I = [I; I_old(i)*ones(num_eles, 1)];
end


%% assemble the rest stuffs

elenum = size(eles, 1);

Supporting_old = Supporting;
Supporting_ind = oldnodes_ind(Supporting_old~=0);
Supporting_type = Supporting_old(Supporting_old~=0);

Supporting = zeros(size(nodes, 1), 1);
Supporting(Supporting_ind) = Supporting_type;

ExF_old = ExF;
ExF = zeros(size(ExF_old));
ExF(:, 1:2) = ExF_old(:, 1:2);
ExF(:, 3) = oldnodes_ind(ExF_old(:, 3));

w_old = w;
w = [];  % [mag eletag]

for i=1:size(w_old, 1)
    orgele_under_w = w_old(i, 2);
    eles_under_w = orgeletag_set(orgeletag_set(:, 1)==orgele_under_w, 2);
    w_mag = w_old(i, 1) * ones(numel(eles_under_w), 1);
    w = [w; [w_mag eles_under_w]];
end

L = [];
for i = 1:size(eles, 1)
    start_node = eles(i, 1);
    end_node = eles(i, 2);
    vec = nodes(end_node, :) - nodes(start_node, :);
    L(i) = norm(vec);
end
L = L';


%% equivalant external element element end forces (LOCAL)
% ExEF = [Fyi Mzi Fyj Mzj nodetag]
w_mag = w(:, 1);
w_tag = w(:, 2);
L_w = L(w(:, 2));
ExEF = [w_mag.*L_w/2 w_mag.*L_w.^2/12 w_mag.*L_w/2 -w_mag.*L_w.^2/12 w_tag];

% element freedom table, each node has 3 DoFs
EFT = 1:1:3*size(nodes, 1);

%% Construct K_ff
K = zeros(size(EFT, 2));

for i = 1:size(eles, 1)
    start_node = eles(i, 1);
    end_node = eles(i, 2);
    vec = nodes(end_node, :) - nodes(start_node, :);
    th = asin(vec(2)/norm(vec));
    if vec(1) < 0
        th = th + pi;
    end
    c = cos(th);
    s = sin(th);
    
    % transformation matrix
    T_e = [c s 0 0 0 0;
           -s c 0 0 0 0;
           0 0 1 0 0 0;
           0 0 0 c s 0;
           0 0 0 -s c 0;
           0 0 0 0 0 1];
       
    A_i = A(i);
    L_i = L(i);
    E_i = E(i);
    I_i = I(i);
    % local stifness matrix for element i
    k_e_b = [A_i*E_i/L_i    0                   0                   -A_i*E_i/L_i    0                   0;
             0              12*E_i*I_i/L_i^3    6*E_i*I_i/L_i^2     0               -12*E_i*I_i/L_i^3   6*E_i*I_i/L_i^2;
             0              6*E_i*I_i/L_i^2     4*E_i*I_i/L_i       0               -6*E_i*I_i/L_i^2    2*E_i*I_i/L_i;
             -A_i*E_i/L_i   0                   0                   A_i*E_i/L_i     0                   0;
             0              -12*E_i*I_i/L_i^3   -6*E_i*I_i/L_i^2    0               12*E_i*I_i/L_i^3    -6*E_i*I_i/L_i^2;
             0              6*E_i*I_i/L_i^2     2*E_i*I_i/L_i       0               -6*E_i*I_i/L_i^2    4*E_i*I_i/L_i;
            ];
    
    % global stifness matrix for element i
    k_e = T_e'*k_e_b*T_e;
    
    % element free table for element i
    EFTe = [EFT(3*start_node - 2:3*start_node), EFT(3*end_node - 2:3*end_node)];
    
    % modify the global stifness matrix
    K(EFTe, EFTe) = K(EFTe, EFTe) + k_e;
end

EFTf = [];
for i = 1:size(Supporting, 1)
    if Supporting(i, 1) == 11
        EFTf = [EFTf EFT(3*i-1) EFT(3*i)];
    elseif Supporting(i, 1) == 12
        EFTf = [EFTf EFT(3*i-2) EFT(3*i)];
    elseif Supporting(i, 1) == 13
        EFTf = [EFTf EFT(3*i-2) EFT(3*i-1)];
    elseif Supporting(i, 1) == 2
        EFTf = [EFTf EFT(3*i)];
    elseif Supporting(i, 1) == 0
        EFTf = [EFTf EFT(3*i-2:3*i)];
    end
end

Kff = K(EFTf, EFTf);

%% Construct f_ff (contained moments)
f = zeros(size(EFT, 2), 1);
for i = 1:size(ExF, 1)
    f_EFT = [3*ExF(i, 3)-2, 3*ExF(i, 3)-1];
    f(f_EFT) = f(f_EFT) + ExF(i, 1:2)';
end

for i = 1:size(ExEF, 1)
    ele = ExEF(i, 5);
    start_node = eles(ele, 1);
    end_node = eles(ele, 2);
    vec = nodes(end_node, :) - nodes(start_node, :);
    th = asin(vec(2)/norm(vec));
    if vec(1) < 0
        th = th + pi;
    end
    c = cos(th);
    s = sin(th);
    L_e = L(ele);
    
    % transformation matrix
    T_e = [c s 0 0 0 0;
           -s c 0 0 0 0;
           0 0 1 0 0 0;
           0 0 0 c s 0;
           0 0 0 -s c 0;
           0 0 0 0 0 1];
    % local equivalant force
    f_EF_i_b = [0 ExEF(i, 1:2) 0 ExEF(i, 3:4)]';
    
    % global equivalant force
    f_EF_i = T_e'*f_EF_i_b;
    
    f_EFT = [(3*start_node-2:3*start_node) (3*end_node-2:3*end_node)];
    
    f(f_EFT) = f(f_EFT) + f_EF_i;
end

f_ff = f(EFTf);

%% Caculate displacement for free degrees
u_ff = Kff^-1*f_ff;
u = zeros(size(EFT, 2), 1);
u(EFTf) = u_ff;
disp('u(m)=');
disp(reshape(u, [3, numel(u)/3]));

%% Calculate element end forces
f_ele_end = zeros(6, size(eles, 1));


for i = 1:size(eles, 1)
    start_node = eles(i, 1);
    end_node = eles(i, 2);
    vec = nodes(end_node, :) - nodes(start_node, :);
    th = asin(vec(2)/norm(vec));
    if vec(1) < 0
        th = th + pi;
    end
    c = cos(th);
    s = sin(th);
    
    % transformation matrix
    T_e = [c s 0 0 0 0;
           -s c 0 0 0 0;
           0 0 1 0 0 0;
           0 0 0 c s 0;
           0 0 0 -s c 0;
           0 0 0 0 0 1];
       
    A_i = A(i);
    L_i = L(i);
    E_i = E(i);
    I_i = I(i);
    % local stifness matrix for element i
    k_e_b = [A_i*E_i/L_i    0                   0                   -A_i*E_i/L_i    0                   0;
             0              12*E_i*I_i/L_i^3    6*E_i*I_i/L_i^2     0               -12*E_i*I_i/L_i^3   6*E_i*I_i/L_i^2;
             0              6*E_i*I_i/L_i^2     4*E_i*I_i/L_i       0               -6*E_i*I_i/L_i^2    2*E_i*I_i/L_i;
             -A_i*E_i/L_i   0                   0                   A_i*E_i/L_i     0                   0;
             0              -12*E_i*I_i/L_i^3   -6*E_i*I_i/L_i^2    0               12*E_i*I_i/L_i^3    -6*E_i*I_i/L_i^2;
             0              6*E_i*I_i/L_i^2     2*E_i*I_i/L_i       0               -6*E_i*I_i/L_i^2    4*E_i*I_i/L_i;
            ];
    
    % global stifness matrix for element i
    k_e = T_e'*k_e_b*T_e;
    
    % element free table for element i
    EFTe = [EFT(3*start_node - 2:3*start_node), EFT(3*end_node - 2:3*end_node)];
    
    % f = k_e*u_e - f_EF_e
    if ismember(i, ExEF(:, 5)) == 0
        f_EF_e = zeros(6, 1);
    else
        f_EF_e_b = ExEF(ismember(ExEF(:, 5), i), :);
        f_EF_e_b = [0 f_EF_e_b(1:2) 0 f_EF_e_b(3:4)]';
        f_EF_e = T_e' * f_EF_e_b;
    end
    u_e = u(EFTe);
    f = k_e * u_e - f_EF_e;
    
    f_ele_end(:, i) = f;
end

disp('Element end forces (global) (N)=');
disp(f_ele_end);

f_ele_end_b = zeros(size(f_ele_end));

for i = 1:size(eles, 1)
    start_node = eles(i, 1);
    end_node = eles(i, 2);
    vec = nodes(end_node, :) - nodes(start_node, :);
    th = asin(vec(2)/norm(vec));
    if vec(1) < 0
        th = th + pi;
    end
    c = cos(th);
    s = sin(th);
    
    % transformation matrix
    T_e = [c s 0 0 0 0;
           -s c 0 0 0 0;
           0 0 1 0 0 0;
           0 0 0 c s 0;
           0 0 0 -s c 0;
           0 0 0 0 0 1];
    
    
    f_ele_end_b(:, i) = T_e * f_ele_end(:, i);
end

disp('Element end forces (local) (N)=');
disp(f_ele_end_b);

%% Calculate reactions
support_nodes = find(Supporting~=0);
reactions = zeros(3 * size(support_nodes, 1), 1);
for i = 1 : size(support_nodes, 1)
    % find which element this support node is
    e = find(eles == support_nodes(i));
    e_ind = ceil(e/2);
    n_ind = 2 - ceil(e_ind - e/2);
    
    % reaction = f_ele_end
    reactions(3*i-2:3*i) = f_ele_end(3*n_ind-2:3*n_ind, e_ind);
    if Supporting(support_nodes(i)) == 2
        reactions(3*i) = NaN;
    elseif Supporting(support_nodes(i)) == 11
        reactions(3*i-1) = NaN;
        reactions(3*i) = NaN;
    elseif Supporting(support_nodes(i)) == 12
        reactions(3*i-2) = NaN;
        reactions(3*i) = NaN;
    elseif Supporting(support_nodes(i)) == 13
        reactions(3*i-2) = NaN;
        reactions(3*i-1) = NaN;
    end
end
disp('Supporting nodes are')
disp(support_nodes)
disp('Reactions(N)=')
disp(reactions)


  
%% visualization
figure()
subplot(221)
title('deformed shape')
for i = 1:elenum
    hold on
    start_node = eles(i, 1);
    end_node = eles(i, 2);
    plot(nodes([start_node, end_node], 1), nodes([start_node, end_node], 2)...
        , 'LineWidth', 1, 'Color', 'red');
end
grid on
xlim([min(min(nodes))-5, max(max(nodes))+5]);
ylim([min(min(nodes))-5, max(max(nodes))+5]);

for i = 1:size(eles, 1)
    hold on
    start_node = eles(i, 1);
    end_node = eles(i, 2);
    
    
    x_disp = u([start_node*3-2, end_node*3-2])*mag_factor;
    y_disp = u([start_node*3 - 1, end_node*3 - 1])*mag_factor;
    
    plot(nodes([start_node, end_node], 1) + x_disp,...
         nodes([start_node, end_node], 2) + y_disp...
        , 'LineWidth', 5, 'Color', 'blue');
end

subplot(222)
for i = 1:elenum
    hold on
    start_node = eles(i, 1);
    end_node = eles(i, 2);
    plot(nodes([start_node, end_node], 1), nodes([start_node, end_node], 2)...
        , 'LineWidth', 1, 'Color', 'red');
end
grid on

subplot(223)
for i = 1:elenum
    hold on
    start_node = eles(i, 1);
    end_node = eles(i, 2);
    plot(nodes([start_node, end_node], 1), nodes([start_node, end_node], 2)...
        , 'LineWidth', 1, 'Color', 'red');
end
grid on

subplot(224)
for i = 1:elenum
    hold on
    start_node = eles(i, 1);
    end_node = eles(i, 2);
    plot(nodes([start_node, end_node], 1), nodes([start_node, end_node], 2)...
        , 'LineWidth', 1, 'Color', 'red');
end
grid on

orgele_num = numel(unique(orgeletag_set(:, 1)));

for i=1:orgele_num
    start_node = eles_org(i, 1);
    end_node = eles_org(i, 2);
    vec = nodes(end_node, :) - nodes(start_node, :);
    th = asin(vec(2)/norm(vec));
    if vec(1) < 0
        th = th + pi;
    end
    c = cos(th);
    s = sin(th);
    R = [c -s; s c];  % rotation matrix
    
    V_set = [];
    
    ele_temp = orgeletag_set(orgeletag_set(:, 1)==i, :);  % orgeletag_set for this element
    V_set = [f_ele_end_b(2, ele_temp(1, 2)) -f_ele_end_b(5, orgeletag_set(:, 1)==i)];
    start_node = eles_org(i, 1);
    end_node = eles_org(i, 2);
    eles_table_plot = nodes([start_node find(nodes_eles_table==i)' end_node], :)';
    V_set = R*[zeros(1, numel(V_set)); V_set]/2/mag_factor + eles_table_plot;
    
    subplot(222)
    title('shear force diagram')
    plot(V_set(1, :), V_set(2, :), 'LineWidth', 2, 'Color', 'blue')
    xlim([min(min(V_set))-5, max(max(V_set))+5]);
    ylim([min(min(V_set))-5, max(max(V_set))+5]);
    
    M_set = [];
    
    M_set = [-f_ele_end_b(3, ele_temp(1, 2)) f_ele_end_b(6, orgeletag_set(:, 1)==i)];
    M_set = R*[zeros(1, numel(M_set)); M_set]/2/mag_factor + eles_table_plot;
    
    subplot(223)
    title('moment diagram')
    plot(M_set(1, :), M_set(2, :), 'LineWidth', 2, 'Color', 'blue')
    xlim([min(min(M_set))-5, max(max(M_set))+5]);
    ylim([min(min(M_set))-5, max(max(M_set))+5]);
    
    
    P_set = [];
    
    P_set = [-f_ele_end_b(1, ele_temp(1, 2)) f_ele_end_b(4, orgeletag_set(:, 1)==i)];
    P_set = R*[zeros(1, numel(P_set)); P_set]/2/mag_factor + eles_table_plot;
    
    subplot(224)
    title('normal force diagram')
    plot(P_set(1, :), P_set(2, :), 'LineWidth', 2, 'Color', 'blue')
    xlim([min(min(P_set))-5, max(max(P_set))+5]);
    ylim([min(min(P_set))-5, max(max(P_set))+5]);
    
    
end

%xlim([min(min(nodes))-5, max(max(nodes))+5]);
%ylim([min(min(nodes))-5, max(max(nodes))+5]);


    














































