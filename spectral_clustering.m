function idx= spectral_clustering(W, k)
    D = diag(sum(W));
    L = D-W;
    opt = struct('issym', true, 'isreal', true);
    [V ,dummy] = eigs(L, D, k, 'SM', opt);
    idx = kmeans(V, k); 
  
    %idx = dbscan(V, 0.24, 21);
end

