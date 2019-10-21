function Y=flatten_RDM(X,diag)
tria_ones=triu(ones(size(X)),diag);
tria_ones(tria_ones==0)=nan;
tria_X=X.*tria_ones;
Y_temp=tria_X(:);
Y_temp(isnan(Y_temp))=[];
Y=Y_temp;
end 