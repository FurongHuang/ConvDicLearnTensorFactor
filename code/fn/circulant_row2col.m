function c = circulant_row2col(r)
a = length(r);
r = r(:)';
c = [r(1),fliplr(r(2:end))];
