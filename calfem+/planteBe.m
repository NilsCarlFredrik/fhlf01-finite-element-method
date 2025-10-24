function [Be A] = planteBe(ex,ey) 
    % Calculates the element B matrix for a three node triangular element.
    % See equation ?? in the report.
    A = 0.5*det([ones(3,1) ex' ey']);
    Be = 1/(2*A)*[ey(2)-ey(3) 0 ey(3)-ey(1) 0 ey(1)-ey(2) 0;
          0 ex(3)-ex(2) 0 ex(1)-ex(3) 0 ex(2)-ex(1);
          ex(3)-ex(2) ey(2)-ey(3) ex(1)-ex(3) ey(3)-ey(1) ex(2)-ex(1) ey(1)-ey(2)];
end