clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% below is the overall structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes = [
    0 0;
    5 5;
    10 5
    ]; % input nodes here, first column of matrix 'nodes' is x cord, second colomn is y

elements = [
            1 2;
            2 3
            ]'; % each ![row]! represents an element, from node a to node b

SupportTypesOnNodes = [2;3;1]; % a vector that indicates how many unknown forces on each node respectively

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% External loads
ExF = [
       1 0 10 5
       ]; % each row: [Fx,Fy,x,y]

ExM = [0, 0, 0]; % [mag, x, y]

% [L_function(l), [start_node, end_node], [direction_x, direction_y]]
syms l
ExL = [
       1, [1, 2], [1, -1];
       2*l, [2, 3], [0, -1]
       ]; 

%% Calculate reaction forces

    
% create unknown variables of reaction forces on each nodes
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


% process ExL, obtain x,y components and its equivalent bearing point
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
    ExL_point_l(i) = int(ExL(i, 1)*l, l, 0, sum((node_to_cor - node_from_cor).^2)^0.5)/...
                  int(ExL(i, 1), l, 0, sum((node_to_cor - node_from_cor).^2)^0.5);
end
    

EqnSet = [];
nodes_if_calculated = ones(size(nodes, 1), 1); % 1 means haven't been calculated, -1 means already
for k = 1:size(elements, 2)
    from_node = elements(1, k);
    to_node = elements(2, k);
    node_func = [(nodes(to_node, 2)-nodes(from_node, 2))/(nodes(to_node, 1)-nodes(from_node, 1))...
                , (nodes(to_node, 2)-nodes(from_node, 2))/(nodes(to_node, 1)-nodes(from_node, 1))*...
                nodes(from_node, 1)+nodes(from_node, 2)]; % [k, b], represent a linear function
    
    Fx = F_react(from_node, 1)*nodes_if_calculated(from_node)...
            + F_react(to_node, 1)*nodes_if_calculated(to_node);
    Fy = F_react(from_node, 2)*nodes_if_calculated(from_node)...
            + F_react(to_node, 2)*nodes_if_calculated(to_node);
    
    %Fx = sum(F_react(:,1)) + sum(ExF(:,1)) + sum(ExL_x);
    %Fy = sum(F_react(:,2)) + sum(ExF(:,2)) + sum(ExL_y);
    Mo = 0;
    for i=1:size(ExM, 1)
        if abs(ExM(i, 2)*node_func(1) + node_func(2) - ExM(i, 3)) < 1e-10
            Mo = Mo + ExM(i);
        end
    end
    
    for i=1:size(ExF, 1)
        if abs(ExF(i, 3) * node_func(1) + node_func(2) - ExF(i, 4)) < 1e-10
            Mo = Mo + cross([ExF(i, 3:4)-nodes(from_node, :) 0], [ExF(i, 1:2) 0]);
            Fx = Fx + ExF(i, 1);
            Fy = Fy + ExF(i, 2);
        end
    end
    for i=1:size(F_react, 1)
        if abs(nodes(i, 1) * node_func(1) + node_func(2) - nodes(i, 2)) < 1e-10
            
            Mo = Mo + cross([nodes(i,:)-nodes(from_node, :) 0], ...
                [F_react(i,1:2)*nodes_if_calculated(i) 0]) + [0 0 F_react(i, 3)*nodes_if_calculated(i)];
        end
    end
    for i=1:size(ExL, 1)
        node_from_cor = nodes(ExL(i, 2), :);
        node_to_cor = nodes(ExL(i, 3), :);
        dir = (node_to_cor - node_from_cor)/sum((node_to_cor - node_from_cor).^2)^0.5;
        r = dir * ExL_point_l(i);
        if abs(min(ExL(i, 2:3)) - min(from_node,to_node)) < 1e-10 && abs(max(ExL(i, 2:3)) - max(from_node,to_node)) < 1e-10
            Mo = Mo + cross([r 0], [ExL_f(i) * ExL(i, 4:5)/sum(ExL(i, 4:5).^2)^0.5 0]);
            Fx = Fx + ExL_x(i);
            Fy = Fy + ExL_y(i);
        end
    end

    EqnSet = [EqnSet; Fx == 0; Fy == 0; Mo(3) == 0];
    nodes_if_calculated(from_node) = -1;
    nodes_if_calculated(to_node) = -1;
end

sol = struct2cell(solve(EqnSet));
F_sol = double(cat(1,sol{:}));































