function entropy= entropy(p)
    if p*(1-p)==0
        entropy=0;
    else
        entropy= -p*log2(p)- (1-p)*log2(1-p);
    end
end