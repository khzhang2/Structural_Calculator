clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% below is the overall structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes = [
    0 0;
    15 0;
    15+15 0
    ]; % input nodes here, first column of matrix 'nodes' is x cord, second colomn is y

elements = [
            1 2;
            2 3
            ]'; % each ![row]! represents an element, from node a to node b

SupportNodes = [1;2]; % nodes that have constrains
        
SupportTypesOnNodes = [2;0;1]; % a vector that indicates how many unknown forces on each node respectively

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% External loads
ExF = [
       0 9375 15 0
       ]; % each row: [Fx,Fy,x,y]

ExM = [0, 0, 0]; % [mag, x, y]

% [L_function(l), [start_node, end_node], [direction_x, direction_y]]
syms l
ExL = [
       500, [1, 3], [0, -1]
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
    ExL_point_l(i) = int(ExL(i, 1)*l, l, 0, sum((node_to_cor - node_from_cor).^2)^0.5)/...
                  int(ExL(i, 1), l, 0, sum((node_to_cor - node_from_cor).^2)^0.5);
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




% %% Calculate shear and moment diagrams
% 
% EqnSet = [];
% F_react_if_calculated = ones(size(F_react, 1), 1); % 1 means haven't been calculated, -1 means already
% ExM_if_calculated = ones(size(ExM, 1), 1);
% ExF_if_calculated = ones(size(ExF, 1), 1);
% ExL_if_calculated = ones(size(ExL, 1), 1);
% V = {}; % shear forces at each element
% M = {}; % internal moments at each element
% 
% 
% 
% 
% 
% for k = 1:size(elements, 2)
%     from_node = elements(1, k);
%     to_node = elements(2, k);
%     l_dir = nodes(to_node, :) - nodes(from_node, :);
%     
%     node_func = [(nodes(to_node, 2)-nodes(from_node, 2))/(nodes(to_node, 1)-nodes(from_node, 1))...
%                 , -(nodes(to_node, 2)-nodes(from_node, 2))/(nodes(to_node, 1)-nodes(from_node, 1))*...
%                 nodes(from_node, 1)+nodes(from_node, 2)]; % [k, b], represent a linear function
%     
%     Fx = F_react(from_node, 1)*F_react_if_calculated(from_node)...
%             + F_react(to_node, 1)*F_react_if_calculated(to_node);
%     Fy = F_react(from_node, 2)*F_react_if_calculated(from_node)...
%             + F_react(to_node, 2)*F_react_if_calculated(to_node);
%     
%     
%     % calculate reaction forces at each nodes
%     Mo = 0;
%     for i=1:size(ExM, 1)
%         if abs(ExM(i, 2)*node_func(1) + node_func(2) - ExM(i, 3)) < 1e-10
%             Mo = Mo + ExM(i);
%         end
%     end
%     
%     for i=1:size(ExF, 1)
%         if abs(ExF(i, 3) * node_func(1) + node_func(2) - ExF(i, 4)) < 1e-10
%             Mo = Mo + cross([ExF(i, 3:4)-nodes(from_node, :) 0], [ExF(i, 1:2) 0]);
%             Fx = Fx + ExF(i, 1);
%             Fy = Fy + ExF(i, 2);
%         end
%     end
%     for i=1:size(F_react, 1)
%         if abs(nodes(i, 1) * node_func(1) + node_func(2) - nodes(i, 2)) < 1e-10
%             Mo = Mo + cross([nodes(i,:)-nodes(from_node, :) 0], ...
%                 [F_react(i,1:2)*F_react_if_calculated(i) 0]) + [0 0 F_react(i, 3)*F_react_if_calculated(i)];
%         end
%     end
%     for i=1:size(ExL, 1)
%         node_from_cor = nodes(ExL(i, 2), :);
%         node_to_cor = nodes(ExL(i, 3), :);
%         dir = (node_to_cor - node_from_cor)/sum((node_to_cor - node_from_cor).^2)^0.5;
%         r = dir * ExL_point_l(i);
%         if abs(min(ExL(i, 2:3)) - min(from_node,to_node)) < 1e-10 &&...
%                 abs(max(ExL(i, 2:3)) - max(from_node,to_node)) < 1e-10
%             Mo = Mo + cross([r 0], [ExL_f(i) * ExL(i, 4:5)/sum(ExL(i, 4:5).^2)^0.5 0]);
%             Fx = Fx + ExL_x(i);
%             Fy = Fy + ExL_y(i);
%         end
%     end
% 
%     EqnSet = [EqnSet; Fx == 0; Fy == 0; Mo(3) == 0];
%     F_react_if_calculated(from_node) = -1;
%     F_react_if_calculated(to_node) = -1;
%     
%     
%     
%     % calculate internal moment and shear forces
%     % initialize variables
%     syms ml vl
%     Meq = 0;
%     Veq = 0;
%     
%     for i = 1:size(ExM, 1)
%         if abs(ExM(i, 2)*node_func(1) + node_func(2) - ExM(i, 3)) < 1e-10 && ExM_if_calculated(i) == 1
%             v_temp = ExM(i, 2:3) - nodes(from_node, :); % direction from start to ExM
%             M_exert_point = sum(v_temp.^2)^0.5;
%             %Meq = Meq + ExM(i, 1) * (l>=0 & l<M_exert_point);
%             Veq = Veq + 0;
%             ExM_if_calculated(i) = -1;
%         end
%         
%     end
%     
%     for i = 1:size(ExF, 1)
%         if abs(ExF(i, 3) * node_func(1) + node_func(2) - ExF(i, 4)) < 1e-10 && ExF_if_calculated(i) == 1
%             v_temp = ExF(i, 3:4) - nodes(from_node, :); % direction from start to ExF
%             F_exert_point = sum(v_temp.^2)^0.5;
%             ExF_dir = ExF(i, 2:3);
%             ExF_dir = ExF_dir/norm(ExF_dir);
%             %Meq = Meq + cross([ExF(i, 3:4)-nodes(from_node, :) 0], [ExF(i, 1:2) 0]) * (l>=0 & l<=F_exert_point);
%             ExF_v = ExF(i,1)*ExF_dir - dot(ExF(i,1)*ExF_dir, v_temp/norm(v_temp))*v_temp/norm(v_temp)*(l>=0 & l<=F_exert_point);
%             Veq = Veq + ExF_v(2);
%             ExF_if_calculated(i) = -1;
%         end
%     end
%     
%     for i = size(ExL, 1)
%         if abs(min(ExL(i, 2:3)) - min(from_node,to_node)) < 1e-10 &&...
%                 abs(max(ExL(i, 2:3)) - max(from_node,to_node)) < 1e-10 &&...
%                 ExL_if_calculated(i) == 1            
%             ExL_dir = ExL(i, 2:3);
%             ExL_dir = ExL_dir/norm(ExL_dir);
%             ExL_mag = ExL(i, 1); % related to l
%             ExL_perp = ExL_mag*ExL_dir - dot(ExL_mag*ExL_dir, l_dir/norm(l_dir))*l_dir/norm(l_dir);
%             ExL_perp = ExL_perp(2);
%             %Meq = Meq + cross([ExF(i, 3:4)-nodes(from_node, :) 0], [ExF(i, 1:2) 0]) * (l>=0 & l<=F_exert_point);
%             Veq = Veq + int(ExL_perp, l, 0, sum(l_dir.^2)^0.5);
%             ExL_if_calculated(i) = -1;
%         end
%     end
%     Veq = Veq - vl*l;
%     V{i} = solve(Veq==0);
%     M{i} = int(V{i}, l);
% end

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









