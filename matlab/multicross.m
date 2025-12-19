function c = multicross(a,b)
  % a is [Nx3], b is [3x1].
  % returns a vector of cross products
  c = [a(:,2)*b(3)-a(:,2)*b(2), a(:,2)*b(1)-a(:,1)*b(3), a(:,1)*b(2)-a(:,2)*b(1)];
end

