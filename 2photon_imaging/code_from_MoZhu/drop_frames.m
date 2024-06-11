% change parameters!!!
x_off = find(ops.xoff>20 | ops.xoff<-20);  %35/20
y_off = find(ops.yoff>20 | ops.yoff<-20);
xy_off = sort([x_off,y_off]);
xy_uniq = unique(xy_off);

% replace uncontinuous dropped frames with regression value
noncont_all = [];
for i = 1:length(xy_uniq)-2
    if xy_uniq(i+1)-xy_uniq(i)~=1 && xy_uniq(i+1)-xy_uniq(i+2)~=-1
        noncont = xy_uniq(i+1);
        noncont_all(end+1) = noncont;
    end
end
noncont_all = [xy_uniq(1),noncont_all];
for j = 2:length(noncont_all)-1
    for k = 1:size(F,1)
        F(k,noncont_all(j)) = mean([F(k,noncont_all(j)+1),F(k,noncont_all(j)-1)]);
        Fneu(k,noncont_all(j)) = mean([Fneu(k,noncont_all(j)+1),Fneu(k,noncont_all(j)-1)]);
        ops.xoff(1,noncont_all(j)) = mean([ops.xoff(1,noncont_all(j)+1),ops.xoff(1,noncont_all(j)-1)]);
        ops.yoff(1,noncont_all(j)) = mean([ops.yoff(1,noncont_all(j)+1),ops.yoff(1,noncont_all(j)-1)]);
    end
end

% replace 2 continuous dropped frames with regression value
cont2_all = [];
diff_uniq = diff(xy_uniq)==1;
for m = 1:length(diff_uniq)-2
    if diff_uniq(m+1)-diff_uniq(m)==1 && diff_uniq(m+1)-diff_uniq(m+2)==1
        cont2 = xy_uniq(m+1);
        cont2_all(end+1) = cont2;
    end
end
for n = 1:length(cont2_all)
    for k = 1:size(F,1)
        F(k,cont2_all(n)) = mean([F(k,cont2_all(n)-1),F(k,cont2_all(n)+2)]);
        F(k,cont2_all(n)+1) = mean([F(k,cont2_all(n)-1),F(k,cont2_all(n)+2)]);
        Fneu(k,cont2_all(n)) = mean([Fneu(k,cont2_all(n)-1),Fneu(k,cont2_all(n)+2)]);
        Fneu(k,cont2_all(n)+1) = mean([Fneu(k,cont2_all(n)-1),Fneu(k,cont2_all(n)+2)]);
    end
end
cont2_num = sort([cont2_all,cont2_all+1]);

cont3 = sort([setdiff(xy_uniq,[noncont_all,cont2_num])]);
zero_vec = zeros(1,size(F,2));
zero_vec(cont3)=1;
% zero_mat = reshape(zero_vec,[size(F,2)/32,32]);
zero_mat = reshape(zero_vec,[size(F,2)/1,1]);


col_all = {};
for p = 1:size(zero_mat,2)
    col = find(zero_mat(:,p) == 1);
    col_all{end+1} = col;
end


%%
FACC6 = F(:,30611:36732);
FSAT2 = F(:,48977:55098);
F(:,30611:36732) = FSAT2;
F(:,48977:55098) = FACC6;

FneuACC6 = Fneu(:,30611:36732);
FneuSAT2 = Fneu(:,48977:55098);
Fneu(:,30611:36732) = FneuSAT2;
Fneu(:,48977:55098) = FneuACC6;

XACC6 = ops.xoff(:,30611:36732);
XSAT2 = ops.xoff(:,48977:55098);
ops.xoff(:,30611:36732) = XSAT2;
ops.xoff(:,48977:55098) = XACC6;

YACC6 = ops.yoff(:,30611:36732);
YSAT2 = ops.yoff(:,48977:55098);
ops.yoff(:,30611:36732) = YSAT2;
ops.yoff(:,48977:55098) = YACC6;







