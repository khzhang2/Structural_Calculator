clear;clc;hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% below is the overall structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes = [
    4, 0;
    4, 3;
    0, 3;
    0, 0
    ]; % input nodes here, first column of matrix 'nodes' is x cord, second colomn is y

elements = [1 1 2 2 3;
            2 4 4 3 4]; % each column represents an element, from node a to node b

SupportTypesOnNodes = [0; 0; 2; 2]; % a vector that indicates how many unknown forces on each node respectively

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExF = [
       0, -6, 4, 0;
       ]; % each row: [Fx,Fy,x,y]

ExM = [0, 0, 0]; % [mag, x, y]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(elements,2)*3 < sum(SupportTypesOnNodes)
    disp('Structure indeterminante!')
    return
end

F_sol = zeros(sum(SupportTypesOnNodes),1);

% reaction due to concentrate forces
for i=1:length(ExF(:,1))
    
    % create unknown variables of reaction forces on each nodes
    syms F [length(SupportTypesOnNodes) 3] % each row cooresponds a node
    
    % release free degrees
    for j=1:length(SupportTypesOnNodes)
        if SupportTypesOnNodes(j)==1
            F(j,1)=0;
            F(j,3)=0;
        elseif SupportTypesOnNodes(j)==2
            F(j,3)=0;
        elseif SupportTypesOnNodes(j)==0
            F(j,1)=0;
            F(j,2)=0;
            F(j,3)=0;
        end
    end
    
    
    Fx = sum(F(:,1)) + ExF(i,1);
    Fy = sum(F(:,2)) + ExF(i,2);
    
    Mo = cross([ExF(i, 3:4) 0], [ExF(i, 1:2) 0]);
    for k=1:size(F,1)
        Mo = Mo + cross([nodes(k,:) 0], [F(k,1:2) 0]);
    end
    sol = struct2cell(solve([ Fx == 0; Fy == 0; Mo(3) == 0]));
    F_sol = F_sol + double(cat(1,sol{:}));
end


% reaction due to moments
for i=1:length(ExM(:,1))
    
    % create unknown variables of reaction forces on each nodes
    syms F [length(SupportTypesOnNodes) 3] 
    
    % release free degrees
    for j=1:length(SupportTypesOnNodes)
        if SupportTypesOnNodes(j)==1
            F(j,1)=0;
            F(j,3)=0;
        elseif SupportTypesOnNodes(j)==2
            F(j,3)=0;
        elseif SupportTypesOnNodes(j)==0
            F(j,1)=0;
            F(j,2)=0;
            F(j,3)=0;
        end
    end
    
    
    Fx = sum(F(:,1));
    Fy = sum(F(:,2));
    Mo = [0, 0, ExM(i,1)];
    for j = 1:length(size(ExF,1))
        for k=1:size(F,1)
            Mo = Mo + cross([nodes(k,:) 0], [F(k,1:2) 0]);
        end
    end
    sol = struct2cell(solve([ Fx == 0; Fy == 0; Mo(3) == 0]));
    F_sol = F_sol + double(cat(1,sol{:}));
end

FF=F';
n_0_ind = find(FF~=0);
FF(n_0_ind) = F_sol;
F_res = FF';
c=1;
for i=1:size(F, 1)
    for j=1:size(F, 2)
        if F(i,j) ~=0
            disp([char(F(i,j)), string(F_sol(c))]);
            c = c+1;
        end
    end
end
    


gen_fig(nodes, elements, SupportTypesOnNodes)
total_F = [ExF; [F_res(:,1:2) nodes]]; % this is the final thing that is interested
draw_loads(total_F, ExM)


syms F_internal [size(elements, 2), size(elements, 2)] % #internal forces equals # elements

EqnSet = [];
ele_mat = zeros(size(elements, 2), size(elements, 2));
for i=1:size(elements, 2)
    ele_mat(elements(1,i), elements(2,i)) = 1;
end

F_internal = F_internal.*ele_mat;

for i=1:size(nodes,1)
    node_x = nodes(i,1);
    node_y = nodes(i,2);
    
    % find the external forces that exert on this node, return the indes of
    % node(s)
    ExF_ind = find(abs(total_F(:,3)-node_x)<1e-10 & abs(total_F(:,4)-node_y)<1e-10);
            
    ExF = total_F(ExF_ind, :);
    ExF = sum(ExF(:, 1:2), 1);
    
    % find the elements on this node (OTN), return the column index of matrix 
    % "elements"
    out = find(elements(1, :)==i);
    out = [out; ones(1, length(out))];
    elements_ONT_out = [elements(:, out(1, :)); ones(1, size(out, 2))]; % the third row=1 means it points out
    
    in = find(elements(2, :)==i);
    in = [in; -ones(1, length(in))];
    elements_ONT_in = [elements(:, in(1, :)); -ones(1, size(in, 2))]; % the third row=-1 means it points in
    
    elements_ONT = [elements_ONT_out, elements_ONT_in]; % each column: from node, to node
    
    % F_internal_ONT: each row: Fx, Fy
    syms F_internal_ONT [0,0]
    for j = 1:size(elements_ONT, 2)
        F_internal_ONT(j, :) = elements_ONT(3,j) * ... % mean out/in
                                (nodes(elements_ONT(2,j), :) - nodes(elements_ONT(1,j), :))... % direction
                                /norm((nodes(elements_ONT(2,j), :) - nodes(elements_ONT(1,j), :))) * ... % normalization
                                F_internal(elements_ONT(1,j), elements_ONT(2,j));
    end
        
    % Fx = 0, Fy = 0
    newEqn = sum(F_internal_ONT, 1)+ExF==[0,0];
    EqnSet = [EqnSet; newEqn(:)];
end

sol = struct2cell(solve(EqnSet));
sol = double(cat(1,sol{:}));

F_internal_result = zeros(size(F_internal));
c=1;
for i=1:size(F_internal_result,1)
    for j=1:size(F_internal_result,2)
        if F_internal(i,j)==0
            continue
        else
            TorC = '';
            if sol(c)>0
                TorC = 'Tension';
            elseif sol(c)<0
                TorC = 'Compression';
            else
                TorC = 'zero';
            end
            
            disp([char(F_internal(i,j)), string(sol(c)), TorC]);
            F_internal_result(i,j) = sol(c);
            c = c+1;
        end
    end
end











function gen_fig(nodes, elements, SupportTypesOnNodes)
    hold on
    for i=1:size(nodes, 1)
        if SupportTypesOnNodes(i)==2
            plot(nodes(i,1), nodes(i,2), 'Marker', '^', 'MarkerFaceColor', 'black', 'MarkerSize', 10)
        elseif SupportTypesOnNodes(i)==1
            plot(nodes(i,1), nodes(i,2), 'Marker', 'o', 'MarkerFaceColor', 'black', 'MarkerSize', 10)
        end
    end
    
    for i=1:size(elements, 2)
        plot( nodes(elements(:, i), 1), nodes(elements(:, i), 2), 'b' ) % plot each element
    end
    
    xlim([-2, 2 + max(nodes(:,1))])
    ylim([-2 + min(nodes(:,2)), 2 + max(nodes(:,2))])
    grid on
end

function draw_loads(total_F, ExM)
    M_mag = abs(ExM(:,1));
    if sum(M_mag) > 1e-5
        max_M = max(M_mag);
        for i=size(ExM,1)
            if ExM(i,1)<=0
                plot(ExM(i,2),ExM(i,3), 'Marker', 'x', 'Color', 'g', 'MarkerSize', 20*M_mag(i)/max_M, 'LineWidth', 4*M_mag(i)/max_M)
            elseif ExM(i,1)>0
                plot(ExM(i,2),ExM(i,3), 'Marker', 'o', 'Color', 'r', 'MarkerSize', 20*M_mag(i)/max_M, 'LineWidth', 4*M_mag(i)/max_M)
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
                plot(linspace(min(total_F(i,3), total_F(i,3) - total_F(i,1)/max_F),...
                    max(total_F(i,3), total_F(i,3) - total_F(i,1)/max_F),50),...
                    linspace(min(total_F(i,4), total_F(i,4)-total_F(i,2)/max_F),...
                    max(total_F(i,4), total_F(i,4)-total_F(i,2)/max_F),50), 'g.')
                plot(total_F(i,3), total_F(i,4), 'Marker', 'd', 'Color', 'g', 'MarkerSize', 7, 'MarkerFaceColor', 'green')
            elseif Mo(3)>=0
                plot(linspace(min(total_F(i,3), total_F(i,3) - total_F(i,1)/max_F),...
                    max(total_F(i,3), total_F(i,3) - total_F(i,1)/max_F),50),...
                    linspace(min(total_F(i,4), total_F(i,4)-total_F(i,2)/max_F),...
                    max(total_F(i,4), total_F(i,4)-total_F(i,2)/max_F),50), 'r.')
                plot(total_F(i,3), total_F(i,4), 'Marker', 'd', 'Color', 'r', 'MarkerSize', 7, 'MarkerFaceColor', 'red')
            end
        end
    end
    
end









