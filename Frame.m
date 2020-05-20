clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% below is the overall structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes = [
    0 0;
    6 0;
    16 0;
    26 0;
    32 0
    ]; % input nodes here, first column of matrix 'nodes' is x cord, second colomn is y

elements = [
            1 2;
            2 3;
            3 4;
            4 5
            ]'; % each ![row]! represents an element, from node a to node b

SupportNodes = [2; 4]; % nodes that have constrains
        
SupportTypesOnNodes = [0 2 0 1 0]'; % a vector that indicates how many unknown forces on each node respectively

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% External loads
ExF = [
       0 -30 32 0
       ]; % each row: [Fx,Fy,x,y]

ExM = [0, 0, 0]; % [mag, x, y]

% [L_function(l), [start_node, end_node], [direction_x, direction_y]]
syms l
ExL = [
       0.8, [4, 5], [0, -1];
       0.8, [1, 2], [0, -1]
       ]; 



    
%% create unknown variables of reaction forces on each nodes
syms F_react [length(SupportTypesOnNodes) 3] % each row cooresponds a node

% release free degrees
for j=1:length(SupportTypesOnNodes)
    if SupportTypesOnNodes(j)==1
        F_react(j,1)=0;
        F_react(j,3)=0;
    elseif SupportTypesOnNodes(j)==2
        F_react(j,3)=0;
    elseif SupportTypesOnNodes(j)==0
        F_react(j,1)=0;
        F_react(j,2)=0;
        F_react(j,3)=0;
    end
end


%% process ExL, obtain x,y components and its equivalent bearing point
ExL_f = zeros(size(ExL, 1), 1);
ExL_x = zeros(size(ExL, 1), 1);
ExL_y = zeros(size(ExL, 1), 1);
ExL_point_l = zeros(size(ExL, 1), 1);
for i=1:size(ExL, 1)
    node_from_cor = nodes(ExL(i, 2), :);
    node_to_cor = nodes(ExL(i, 3), :);
    ExL_f(i) = int(ExL(i, 1), l, 0, sum((node_to_cor - node_from_cor).^2)^0.5);
    ExL_x(i) = ExL_f(i) * ExL(i, 4)/sum(ExL(i, 4:5).^2)^0.5;
    ExL_y(i) = ExL_f(i) * ExL(i, 5)/sum(ExL(i, 4:5).^2)^0.5;
    
    % the equivalent force position, respect to l, from node_start
    if ExL(i,1) ~= 0
        ExL_point_l(i) = int(ExL(i, 1)*l, l, 0, sum((node_to_cor - node_from_cor).^2)^0.5)/...
                       int(ExL(i, 1), l, 0, sum((node_to_cor - node_from_cor).^2)^0.5);
    else
        ExL_point_l(i) = 0;
    end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Calculate reaction forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fx = sum(F_react(:,1)) + sum(ExF(:,1));
Fy = sum(F_react(:,2)) + sum(ExF(:,2));

Mo = sum(ExM(:, 1),1);
for i=1:size(ExF, 1)
    Mo = Mo + cross([ExF(i, 3:4) 0], [ExF(i, 1:2) 0]);
end
for k=1:size(F_react,1)
    Mo = Mo + cross([nodes(k,:) 0], [F_react(k,1:2) 0]) + [0 0 F_react(k, 3)];
end
for i=1:size(ExL, 1)
    node_from_cor = nodes(ExL(i, 2), :);
    node_to_cor = nodes(ExL(i, 3), :);
    dir = (node_to_cor - node_from_cor)/sum((node_to_cor - node_from_cor).^2)^0.5;
    r = dir * ExL_point_l(i) + node_from_cor;
    
    Mo = Mo + cross([r 0], [ExL_f(i) * ExL(i, 4:5)/sum(ExL(i, 4:5).^2)^0.5 0]);
    Fx = Fx + ExL_x(i);
    Fy = Fy + ExL_y(i);
end

EqnSet = [ Fx == 0; Fy == 0; Mo(3) == 0];


s = solve(EqnSet);
F_name = fieldnames(s);
sol = struct2cell(s);
sol = double(cat(1,sol{:}));


F_react_f = F_react(:, 1:2); % matrix of reaction forces
react_M = F_react(:, 3); % matrix of reaction moments
NumOfReactF = length(F_react_f(F_react_f~=0)) + length(react_M(react_M~=0));% number of reaction forces

F_react_f_res = zeros(size(F_react_f));
react_M_res = zeros(size(react_M));
for i=1:NumOfReactF
    F_name_this = F_name{i};
    a_ind = strfind(F_name_this, 't'); % index of 't'
    b_ind = strfind(F_name_this, '_'); % index of '_'
    b_ind = b_ind(2);
    node = str2double(F_name_this(a_ind+1:b_ind-1));
    dir = str2double(F_name_this(b_ind+1:length(F_name_this))); % reaction force direction
    if dir ~= 3
        F_react_f_res(node, dir) = sol(i);
    elseif dir == 3
        react_M_res(node) = sol(i);
    end
end


%% generate figure to illustrate the problem
F_react = [F_react_f_res, nodes];
react_M = [react_M_res, nodes];
total_F = [ExF; F_react];
total_M = [ExM; react_M];

gen_fig(nodes, elements, SupportTypesOnNodes)
draw_loads(nodes, total_F, total_M, ExL)



%% display forces in command line
for i=1:size(sol,1)
    F_name_this = F_name{i};
    if_m = '';
    if F_name_this(length(F_name_this)) == '3'
        F_name_this = [F_name_this(1:length(F_name_this)-1) 'M'];
        if_m = 'Moment';
    end
    disp([F_name_this, string(sol(i)), if_m]);
end














%% functions need to be revised
function gen_fig(nodes, elements, SupportTypesOnNodes)
    hold on
    xlim([-2+min(nodes(:,1)), 2 + max(nodes(:,1))])
    ylim([-2 + min(nodes(:,2)), 2 + max(nodes(:,2))])
    grid on

    
    for i=1:size(elements, 2)
        plot( nodes(elements(:, i), 1), nodes(elements(:, i), 2), 'color', 'b', 'LineWidth', 5) 
    end
    
    for i=1:size(nodes, 1)
        if SupportTypesOnNodes(i)==2
            plot(nodes(i,1), nodes(i,2), 'Marker', '^', 'MarkerFaceColor', 'black', 'MarkerSize', 10)
        elseif SupportTypesOnNodes(i)==1
            plot(nodes(i,1), nodes(i,2), 'Marker', 'o', 'MarkerFaceColor', 'black', 'MarkerSize', 10)
        end
    end
end

function draw_loads(nodes, total_F, total_M, ExL)
    syms l
    L_funcs = ExL(:, 1);
    L_max_mag_for_each_L = zeros(size(ExL, 1), 1);
    for i=1:size(ExL, 1)
        node_from_cor = nodes(ExL(i, 2), :);
        node_to_cor = nodes(ExL(i, 3), :);
        l_len = sum((node_to_cor - node_from_cor).^2)^0.5;
        
        L_max_mag_for_each_L(i) = max(abs(subs(L_funcs(i), l_len)), abs(subs(L_funcs(i), 0)));
    end
    
    if sum(L_max_mag_for_each_L) > 1e-5
        max_L = max(L_max_mag_for_each_L);
        for i=1:size(ExL, 1)
            node_from_cor = nodes(ExL(i, 2), :);
            node_to_cor = nodes(ExL(i, 3), :);
            l_len = sum((node_to_cor - node_from_cor).^2)^0.5;
            
            offset_start = subs(L_funcs(i), 0);
            offset_end = subs(L_funcs(i), l_len);
            
            l_dir = ExL(i, 4:5);
            l_dir = l_dir/norm(l_dir);
            
            node_from_cor_offset = node_from_cor - l_dir*offset_start/max_L/2;
            node_to_cor_offset = node_to_cor - l_dir*offset_end/max_L/2;
            x = [node_from_cor(1), node_from_cor_offset(1), node_to_cor_offset(1), node_to_cor(1)];
            y = [node_from_cor(2), node_from_cor_offset(2), node_to_cor_offset(2), node_to_cor(2)];
            
            
            plot(x, y, 'Color', 'black', 'LineWidth', 1)
        end
    end

    M_mag = abs(total_M(:,1));
    if sum(M_mag) > 1e-5
        max_M = max(M_mag);
        for i=1:size(total_M,1)
            if total_M(i,1)<0
                plot(total_M(i,2),total_M(i,3), 'Marker', 'x', 'Color', 'r', 'MarkerSize', 10+20*M_mag(i)/max_M, 'LineWidth', 4*M_mag(i)/max_M)
            elseif total_M(i,1)>0
                plot(total_M(i,2),total_M(i,3), 'Marker', 'o', 'Color', 'b', 'MarkerSize', 10+20*M_mag(i)/max_M, 'LineWidth', 4*M_mag(i)/max_M)
            end
        end
    end
    
    F_mag = sqrt(sum((total_F(:,1:2).^2)')');
    if sum(F_mag) > 1e-5
        max_F = max(F_mag)*0.8;
        for i=1:size(total_F, 1)
            if total_F(i,1)==0&&total_F(i,2)==0
                continue
            end
            Mo = cross([total_F(i, 3:4) 0], [total_F(i, 1:2) 0]);
            
            if Mo(3)<0
                plot(linspace(total_F(i,3) - total_F(i,1)/max_F,...
                    total_F(i,3),50),...
                    linspace(total_F(i,4)-total_F(i,2)/max_F,...
                    total_F(i,4),50), 'r.')
                plot(total_F(i,3), total_F(i,4), 'Marker', 'd', 'Color', 'r', 'MarkerSize', 7, 'MarkerFaceColor', 'r')
            elseif Mo(3)>0
                plot(linspace(total_F(i,3) - total_F(i,1)/max_F,...
                    total_F(i,3),50),...
                    linspace(total_F(i,4)-total_F(i,2)/max_F,...
                    total_F(i,4),50), 'b.')
                plot(total_F(i,3), total_F(i,4), 'Marker', 'd', 'Color', 'b', 'MarkerSize', 7, 'MarkerFaceColor', 'b')
            else
                plot(linspace(total_F(i,3) - total_F(i,1)/max_F,...
                    total_F(i,3),50),...
                    linspace(total_F(i,4)-total_F(i,2)/max_F,...
                    total_F(i,4),50), 'g.')
                plot(total_F(i,3), total_F(i,4), 'Marker', 'd', 'Color', 'green', 'MarkerSize', 7, 'MarkerFaceColor', 'green')
            end
        end
    end
    
end









