function val = interpVal(x1, y1, x2, y2)
    val = x1+(x2-x1)*(1-y1)/(y2-y1);
end