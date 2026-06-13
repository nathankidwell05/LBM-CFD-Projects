function C_w = find_the_wall_point(x1,y1,x2,y2,R,center_x,center_y)

    % Check if Point A is inside 
    s1 = test_circle(x1,y1,R,center_x,center_y);
    % check if Point B is inside 
    s2 = test_circle(x2,y2,R,center_x,center_y);
    % error if both points are inside or outside
    if s1 == s2
        error('Point A and Point B must be on opposite sides of the circle.')
    end

    % direction from B to A
    dx = x1 - x2;
    dy = y1 - y2;

    % coefficients
    a = dx^2 + dy^2;
    b = 2*((x2 - center_x)*dx + (y2 - center_y)*dy);
    c = (x2 - center_x)^2 + (y2 - center_y)^2 - R^2;

    % solve quadratic
    t_values = roots([a b c]);

    % t value between 0 and 1
    for i = 1:length(t_values)
        if t_values(i) >= 0 && t_values(i) <= 1
            t = t_values(i);
        end
    end

    % Find wall point
    x = x2 + t*dx;
    y = y2 + t*dy;

    C_w=[x;y];

end