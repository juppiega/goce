function satCoverage(beginStr, endStr, ch, gr)

t1 = datenum(beginStr);
t2 = datenum(endStr);

chCover = sum(t1<ch & ch<t2) * mode(diff(ch(t1<ch & ch<t2))) / (t2-t1)
grCover = sum(t1<gr & gr<t2) * median(diff(gr(1:1E5))) / (t2-t1)

end