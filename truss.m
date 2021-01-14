clear; clc; close all;

nodes = [0 0;
         0 0.5;
         0.5 0.5;
         1 0.5;
         1 0;
         0.75 0;
         0.5  0;
         0.25 0;
         0.25 0.25;
         0.75 0.25];
ele = [1 2;
       2 3;
       3 4;
       4 5;
       5 6;
       6 7;
       7 8;
       8 9;
       1 8;
       1 9;
       3 7;
       3 9;
       3 10;
       5 10;
       6 10;
       7 9;
       7 10;
       ];

% Geometry properties
E = 200e9*ones(size(ele, 1));
A = 0.002*ones(size(ele, 1));


% 2 mead xy constrain; 12 means y constrain, 11 means x constrain
Supporting = [2 0 0 0 12 0 0 0 0 0]'; 

ExF = [5000 0 2;
       0 -10000 3;
       0 -5000 4;
       0 -5000 6;
       0 -5000 8]; % each row: [Fx,Fy,node]
ExM = []; % [mag, node]
% [L_function(l), [start_node, end_node], [direction_x, direction_y]]
syms l
ExL = [];

% element freedom table, each node has two DoFs
EFT = 1:1:2*size(ele, 1);


%% Construct K_ff
K = zeros(size(EFT, 2));

for i = 1:size(ele, 1)
    start_node = ele(i, 1);
    end_node = ele(i, 2);
    vec = nodes(end_node, :) - nodes(start_node, :);
    L = norm(vec);
    th = atan(vec(2)/vec(1));
    if vec(1) < 0
        th = th + pi;
    end
    c = cos(th);
    s = sin(th);
    
    % element free table for element i
    EFTe = [EFT(2*start_node - 1:2*start_node), EFT(2*end_node - 1:2*end_node)];
    
    % element stifness matrix (global coordinate)
    ke = E(i)*A(i)/L*[-c -s c s]'*[-c -s c s];
    
    % modify the global stifness matrix
    K(EFTe, EFTe) = K(EFTe, EFTe) + ke;
end

EFTf = [];
for i = 1:size(Supporting, 1)
    if Supporting(i, 1) == 11
        EFTf = [EFTf EFT(2*i)];
    elseif Supporting(i, 1) == 12
        EFTf = [EFTf EFT(2*i-1)];
    elseif Supporting(i, 1) == 0
        EFTf = [EFTf EFT(2*i-1:2*i)];
    end
end

Kff = K(EFTf, EFTf);

%% Construct f_ff
f = zeros(size(EFT, 2), 1);
for i = 1:size(ExF, 1)
    f_EFT = [2*ExF(i, 3)-1, 2*ExF(i, 3)];
    f(f_EFT) = f(f_EFT) + ExF(i, 1:2)';
end

f_ff = f(EFTf);

%% Caculate displacement for free degrees
u_ff = Kff^-1*f_ff;
disp('u_ff(N/m)=');
disp(u_ff);

%% Calculate internal forces for elements
nodes_f = (EFTf + mod(EFTf, 2))/2;
f_internal = zeros(size(ele, 1), 1);
u = zeros(size(EFT, 2), 1);
u(EFTf) = u_ff;

for i = 1:size(ele, 1)
    start_node = ele(i, 1);
    end_node = ele(i, 2);
    vec = nodes(end_node, :) - nodes(start_node, :);
    L = norm(vec);
    th = atan(vec(2)/vec(1));
    if vec(1) < 0
        th = th + pi;
    end
    c = cos(th);
    s = sin(th);

    % element free table for element i
    EFTe = [EFT(2*start_node - 1:2*start_node), EFT(2*end_node - 1:2*end_node)];

    % element stifness matrix (global coordinate)
    ke = E(i)*A(i)/L*[-c -s c s]'*[-c -s c s];

    % f = ke*u
    ue = u([2*start_node-1: 2*start_node, 2*end_node-1: 2*end_node]);
    ue_bar = [c s 0 0; 0 0 c s]*ue;
    epsilon = (ue_bar(2) - ue_bar(1)) / L;
    f = E(i)*A(i)*epsilon;
    f_internal(i) = f;
    
end
disp('f_internal(N)=');
disp(f_internal);

%% Calculate the reactions at each supports
support_nodes = find(Supporting~=0);
ele_end_force = zeros(size(support_nodes, 1), 2);
for i = 1:size(support_nodes, 1)
    related_ele = sum(ele == support_nodes(i), 2);
    related_ele_num = find(related_ele == 1);
    f_internal_related = f_internal(related_ele~=0);
    
    for j = 1:size(related_ele_num)
        start_node = ele(related_ele_num(j), 1);
        end_node = ele(related_ele_num(j), 2);
        vec = nodes(end_node, :) - nodes(start_node, :);
        L = norm(vec);
        th = atan(vec(2)/vec(1));
        if vec(1) < 0
            th = th + pi;
        end
        c = cos(th);
        s = sin(th);
        if ele(related_ele_num(j), 1) == support_nodes(i)
            vec = vec;
        elseif ele(related_ele_num(j), 2) == support_nodes(i)
            vec = -vec;
        end
        
        ele_end_force(i, :) = ele_end_force(i, :) + f_internal_related(j)*vec/L; 
    end
end

f_reaction = -ele_end_force;
disp('Supporting nodes:');
disp(support_nodes);
disp('Reaction(N)=');
disp(f_reaction);





























