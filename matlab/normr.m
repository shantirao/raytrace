
function [retval] = normr (input1)
 %input1 = input1';
%  retval = bsxfun(@times,input1,1./sqrt(sum(input1.^2,2)));
  retval = input1 ./ sqrt(sum(input1.^2,2));
end
