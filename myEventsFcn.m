function [v,isterminal,direction] = myEventsFcn(t,y)
v = [sqrt((y(1)-y(4))^2+(y(2)-y(5))^2+(y(3)-y(6))^2)-50,... 
    sqrt((y(1)-y(7))^2+(y(2)-y(8))^2+(y(3)-y(9))^2)-50,... 
    sqrt((y(4)-y(7))^2+(y(5)-y(8))^2+(y(6)-y(9))^2)-50];          % detect when position crosses zero
isterminal = [1 1 1];    % stop the integration
direction = [0 0 0];     % detect all crossings
end