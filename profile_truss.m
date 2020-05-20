clear; clc;
if_internal = 1;
nodes = [
    0 0;
    0 1.5;
    2 1.5;
    2 0
    ]; % input nodes here, first column of matrix 'nodes' is x cord, second colomn is y

elements = [
            1 2;
            2 3;
            3 4;
            4 1;
            1 3
            ]'; % each ![row]! represents an element, from node a to node b

SupportTypesOnNodes = [1 0 0 2]'; % a vector that indicates how many unknown forces on each node respectively

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% External loads
ExF = [
       1 1 0 1.5
       ]; % each row: [Fx,Fy,x,y]

ExM = [0, 0, 0]; % [mag, x, y]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('profile_truss.mat')


F_internal_res = Truss()








