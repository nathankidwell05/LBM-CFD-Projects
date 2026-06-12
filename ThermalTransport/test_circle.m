function s = test_circle(x,y,R,center_x,center_y)

    distance = sqrt((x - center_x)^2 + (y - center_y)^2);

    if distance <= R
        s = 1;
    else
        s = 0;
    end

end