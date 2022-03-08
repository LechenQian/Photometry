function [prt_colors, prt_conditions, num_CSTT, CSTT_mostR, R_colors] = get_protocol_info(protocol)

switch protocol
    case 'Revers'
        prt_conditions{1} = 'Rew \rightarrow Airpuff';
        prt_conditions{2} = 'Rew \rightarrow Nothing';
        prt_conditions{3} = 'Nothing \rightarrow Rew';
        prt_conditions{4} = 'Airpuff \rightarrow Rew';
        prt_conditions{5} = 'FreeRew';
        prt_colors = [[0.74,0.83,0.37]; [0.15,0.74,0]; 0.45*[0,1,1]+0.15; 0.45*[0,1,1]+0.4; [0 0 0]];
        num_CSTT = 4;         % number of trial types with CS
        CSTT_mostR = nan;     % trial type with CS with the most number of reward sizes being delivered
        R_colors = nan;
    case 'DRL6Rv'
        prt_conditions{1} = 'Nothing \rightarrow RewFx';
        prt_conditions{2} = 'Nothing \rightarrow RewD1';
        prt_conditions{3} = 'Nothing \rightarrow RewD2';
        prt_conditions{4} = 'RewFx \rightarrow Nothing';
        prt_conditions{5} = 'RewD1 \rightarrow Nothing';
        prt_conditions{6} = 'RewD2 \rightarrow Nothing';
        prt_conditions{7} = 'FreeRew';
        prt_colors = [0.45*[0,1,1]; 0.45*[0,1,1]+0.15; 0.45*[0,1,1]+0.4; [1,0.8,0.5]; [1,0.7,0.4]-0.1; [1,0.4,0.1]-0.1; [0.15,0.74,0]];
        num_CSTT = 6;         % number of trial types with CS
        CSTT_mostR = [5 6];     % trial type with CS with the most number of reward sizes being delivered      
        % Rewards are: 0; 0; 0; 4; 0.2 0.5 1 2 4 6 7 7.5 7.8; 0.2 0.5 1 2 4 6 7 7.5 7.8
        clrs_R = brewermap(9, 'Paired');  
        clrs_R = cat(1,[0.2 0 0],clrs_R); % these colors will correspond to Rs of sizes 0,0.2,0.5,1,2,4,6,7,7.5,7.8
        % Rearrange the colors to go from blue-green-pink
        aux = clrs_R(6,:);
        aux2 = clrs_R(7,:);
        clrs_R(6,:) = [1 0.9 0.1];      % mean reward (in the middle of the distribution) in orange
        clrs_R(7,:) = min(1,clrs_R(5,:)+0.2);      % higher side of the mean: green
        clrs_R(8,:) = clrs_R(5,:);
        clrs_R(9,:) = max(0,clrs_R(5,:)-0.2);
        clrs_R(10,:) = max(0,clrs_R(5,:)-0.4);
        clrs_R(5,:) = aux;              % lower side of the man: pink
        clrs_R(3,:) = aux2;
        clrs_R(4,:) = mean([clrs_R(3,:);clrs_R(5,:)]);        
        clrs_R(2,:) = max(0,clrs_R(3,:)-0.2);
       
%         clrs_R(2,:) = clrs_R(3,:)-0.1;
%         clrs_R(3,:) = aux-0.2;
%         clrs_R(4,:) = mean([clrs_R(3,:);clrs_R(5,:)]);
%         
%         clrs_R(end,:) = [0.83 0 0.33];
%         clrs_R(7,:) = min(1,clrs_R(8,:)+0.1);
        % Atribute the same colors to rewards of the same size delivered by
        % different TTs
        R_colors{1} = clrs_R(1,:);     
        R_colors{2} = clrs_R(1,:);  
        R_colors{3} = clrs_R(1,:);
        R_colors{4} = clrs_R(6,:);
        R_colors{5} = clrs_R(2:end,:);
        R_colors{6} = clrs_R(2:end,:); 
        R_colors{7} = clrs_R(2:end,:);
    case 'DRL6Od'
        prt_conditions = {'TT1','TT2','TT3','TT4','TT5','TT6','TT7(umpR)'};
        prt_colors = [0.45*[0,1,1]; [1,0.8,0.5]; [1,0.7,0.4]-0.1; [1,0.4,0.1]-0.1; 0.45*[0,1,1]+0.15; 0.45*[0,1,1]+0.4; [0.15,0.74,0]];
        num_CSTT = 6;       
        CSTT_mostR = 4;
        % Rewards are: 0; 2,4,6; 2,6; 1,2,3,4,5,6; 2; 6; 1,2,3,4,5,6
        clrs_R = brewermap(7, 'Paired');  % these colors will correspond to Rs of sizes 0,1,2,3,4,5,6,7
        clrs_R = cat(1,[0 0 0.5],clrs_R);
        % Rearrange the colors to go from blue-green-pink
        aux = clrs_R(2,:);
        clrs_R(2,:) = clrs_R(3,:)-0.1;
        clrs_R(3,:) = aux-0.2;
        clrs_R(4,:) = mean([clrs_R(3,:);clrs_R(5,:)]);
        clrs_R(6,:) = mean([1 0.6 0.7; clrs_R(5,:)]);
        clrs_R(end,:) = [0.83 0 0.33];
        clrs_R(7,:) = [1 0.6 0.7];
        % Atribute the same colors to rewards of the same size delivered by
        % different TTs
        R_colors{1} = clrs_R(1,:);     
        R_colors{2} = [clrs_R(3,:); clrs_R(5,:); clrs_R(7,:)];     
        R_colors{3} = [clrs_R(3,:); clrs_R(7,:)]; 
        R_colors{4} = clrs_R(2:end,:);
        R_colors{5} = clrs_R(3,:);
        R_colors{6} = clrs_R(7,:); 
        R_colors{7} = clrs_R(2:end,:);
    case 'DRL3Od'
        prt_conditions = {'TT1','TT3','TT4'};
        prt_colors = [0.45*[0,1,1];[1,0.7,0.4]-0.1; [1,0.4,0.1]-0.1];
        num_CSTT = 3;      
        CSTT_mostR = 3;   
        % colors for rewards (use the same as in DRL6Od
        clrs_R = brewermap(7, 'Spectral');  % these colors will correspond to Rs of sizes 0,1,2,3,4,5,6,7
        clrs_R = cat(1,[0 0 0],clrs_R);
        R_colors{1} = clrs_R(1,:);     
        R_colors{2} = clrs_R(5,:);     
        R_colors{3} = [clrs_R(2,:); clrs_R(8,:)]; 
    case 'VarMag'
        prt_conditions = {'TT1','TT2','TT3(umpR)'};
        prt_colors = [0.45*[0,1,1]; [0.74,0.83,0.37]; [0.15,0.74,0]];
        num_CSTT = 2;       
        CSTT_mostR = 2;    
        % Rewards are: 0; 0.3 0.5 1.2 2.5 5.0 8.0 11.0; 0.3 0.5 1.2 2.5 5.0 8.0 11.0; 
        clrs_R = brewermap(7, 'Paired');  
        clrs_R = cat(1,[0 0 0],clrs_R);
        % Rearrange the colors to go from blue-green-pink
        aux = clrs_R(2,:);
        clrs_R(2,:) = clrs_R(3,:)-0.1;
        clrs_R(3,:) = aux-0.2;
        clrs_R(4,:) = mean([clrs_R(3,:);clrs_R(5,:)]);
        clrs_R(6,:) = mean([1 0.6 0.7; 1 0.6 0.7; clrs_R(4,:)]);
        clrs_R(end,:) = [0.83 0 0.33];
        clrs_R(7,:) = [1 0.6 0.7];
        % 
        R_colors{1} = clrs_R(1,:);     
        R_colors{2} = clrs_R(2:end,:);
        R_colors{3} = clrs_R(2:end,:);
    case 'Eshl16'
        prt_conditions = {'TT1','TT2','TT3(umpR)'};
        prt_colors = [0.45*[0,1,1]; [0.74,0.83,0.37]; [0.15,0.74,0]];
        num_CSTT = 2;
        CSTT_mostR = 2;
        % Rewards are: 0; 0.1 0.3 1.2 2.5 5.0 10.0 20.0; 0.1 0.3 1.2 2.5 5.0 10.0 20.0
        clrs_R = brewermap(7, 'Paired');
        clrs_R = cat(1,[0 0 0],clrs_R);
        % Rearrange the colors to go from blue-green-pink
        aux = clrs_R(2,:);
        clrs_R(2,:) = clrs_R(3,:)-0.1;
        clrs_R(3,:) = aux-0.2;
        clrs_R(4,:) = mean([clrs_R(3,:);clrs_R(5,:)]);
        clrs_R(6,:) = mean([1 0.6 0.7; clrs_R(5,:)]);
        clrs_R(end,:) = [0.83 0 0.33];
        clrs_R(7,:) = [1 0.6 0.7];
        %
        R_colors{1} = clrs_R(1,:);
        R_colors{2} = clrs_R(2:end,:);
        R_colors{3} = clrs_R(2:end,:);
    case 'Es163O'
        prt_conditions = {'TT1','TT2','TT3','TT4(umpR)'};
        prt_colors = [0.45*[0,1,1]; [0.74,0.83,0.37]; [1,0.4,0.1]-0.1; [0.15,0.74,0]];
        num_CSTT = 3;
        CSTT_mostR = 2;
        % Rewards are: 0; 0.1 0.3 1.2 2.5 5.0 10.0 20.0; 0.1 0.3 1.2 2.5 5.0 10.0 20.0; 5.1
        clrs_R = brewermap(7, 'Paired');
        clrs_R = cat(1,[0 0 0],clrs_R);
        % Rearrange the colors to go from blue-green-pink
        aux = clrs_R(2,:);
        clrs_R(2,:) = clrs_R(3,:)-0.1;
        clrs_R(3,:) = aux-0.2;
        clrs_R(4,:) = mean([clrs_R(3,:);clrs_R(5,:)]);
        clrs_R(6,:) = mean([1 0.6 0.7; clrs_R(5,:)]);
        clrs_R(end,:) = [0.83 0 0.33];
        clrs_R(7,:) = [1 0.6 0.7];
        %
        R_colors{1} = clrs_R(1,:);
        R_colors{2} = clrs_R(2:end,:);
        R_colors{3} = clrs_R(2:end,:);
    case 'WDMSFP'
        prt_conditions = {'TT1'};
        prt_colors = [0.74,0.83,0.37];
        num_CSTT = 0;       
        CSTT_mostR = [];  
         % Rewards are: 0; 0.1 0.3 1.2 2.5 5.0 10.0 20.0; 0.1 0.3 1.2 2.5 5.0 10.0 20.0
        clrs_R = brewermap(7, 'Paired');  
        clrs_R = cat(1,[0 0 0],clrs_R);
        % Rearrange the colors to go from blue-green-pink
        aux = clrs_R(2,:);
        clrs_R(2,:) = clrs_R(3,:)-0.1;
        clrs_R(3,:) = aux-0.2;
        clrs_R(4,:) = mean([clrs_R(3,:);clrs_R(5,:)]);
        clrs_R(6,:) = mean([1 0.6 0.7; clrs_R(5,:)]);
        clrs_R(end,:) = [0.83 0 0.33];
        clrs_R(7,:) = [1 0.6 0.7];
        % 
        R_colors{1} = clrs_R(1,:);     
        R_colors{2} = clrs_R(2:end,:);
        R_colors{3} = clrs_R(2:end,:);
end
end