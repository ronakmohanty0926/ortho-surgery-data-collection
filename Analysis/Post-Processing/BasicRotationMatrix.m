function matrix = BasicRotationMatrix(axis, angle)

matrix = eye(3);

if axis == 'X'
    matrix(2,2) = cos(angle);
    matrix(2,3) = -sin(angle);
    matrix(3,2) = sin(angle);
    matrix(3,3) = cos(angle);
elseif axis == 'Y'
    matrix(1,1) = cos(angle);
    matrix(1,3) = sin(angle);
    matrix(3,1) = -sin(angle);
    matrix(3,3) = cos(angle);
elseif axis == 'Z'
    matrix(1,1) = cos(angle);
    matrix(1,2) = -sin(angle);
    matrix(2,1) = sin(angle);
    matrix(2,2) = cos(angle);
end