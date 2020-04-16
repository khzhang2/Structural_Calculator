clear
nodes = [
    0,0;
    10,5;
    10,0
    ]; % input nodes here, first column of matrix 'nodes' is x cord, second colomn is y

elements = [1 1 2;
            2 3 3]; % each column represents an element

SupportTypesOnNodes = [2; 0; 1]; % a vector that indicates how many unknown forces on each node respectively

ExF = [0, -10, 10, 5;
       0, 5, 5, 2.5;
       0, -10, 5, 0
       ]; % external forces, first column is Fx, second is Fy, third and forth are the position

ExM = [-10, 0, 0]; % [mag, x, y]

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

n_0_ind = find(F~=0);
for i=1:size(n_0_ind)
    disp([char(F(n_0_ind(i))), string(F_sol(i))])
end
    

F(n_0_ind)=F_sol;
F = double(F);

gen_fig(nodes, elements, SupportTypesOnNodes)
total_F = [ExF; [F(:,1:2) nodes]];
draw_loads(total_F, ExM)


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
                plot([min(total_F(i,3), total_F(i,3) - total_F(i,1)/max_F):.01:...
                    max(total_F(i,3), total_F(i,3) - total_F(i,1)/max_F)],...
                    [min(total_F(i,4), total_F(i,4)-total_F(i,2)/max_F): .01 :...
                    max(total_F(i,4), total_F(i,4)-total_F(i,2)/max_F)], 'g.')
                plot(total_F(i,3), total_F(i,4), 'Marker', 'd', 'Color', 'g', 'MarkerSize', 7, 'MarkerFaceColor', 'green')
            elseif Mo(3)>=0
                plot([min(total_F(i,3), total_F(i,3) - total_F(i,1)/max_F):.01:...
                    max(total_F(i,3), total_F(i,3) - total_F(i,1)/max_F)],...
                    [min(total_F(i,4), total_F(i,4)-total_F(i,2)/max_F): .01 :...
                    max(total_F(i,4), total_F(i,4)-total_F(i,2)/max_F)], 'r.')
                plot(total_F(i,3), total_F(i,4), 'Marker', 'd', 'Color', 'r', 'MarkerSize', 7, 'MarkerFaceColor', 'red')
            end
        end
    end
    
end







