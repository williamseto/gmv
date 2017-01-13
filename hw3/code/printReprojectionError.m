function [ error ] = printReprojectionError( graph )

K = f2K(graph.f);
if isfield(graph, 'px')
    K = [graph.f, 0, graph.px; 0, graph.f, graph.py; 0 0 1];
end
error = 0;
for i=1:length(graph.frames)
    
    % logical indicating which 3D points we can observe in this frame
    str_idx = (graph.ObsIdx(i,:) ~= 0);
    
    % indices into the observed feature points, indicating the image
    % features in this frame
    obs_idx = graph.ObsIdx(i, str_idx);
    
    % p is image coordinates (centered)
    p = graph.ObsVal(:, obs_idx);
  
    M = graph.Mot(:,:,i);
    P = graph.Str(:, str_idx);
    p_hat = K*M*[P; ones(1, size(P,2))];
    p_hat = [p_hat(1,:)./p_hat(3,:); p_hat(2,:)./p_hat(3,:)];
    
    error = error + sum(sum((p_hat - p).^2));
%     fprintf('error for frame %d: %f\n', i,  sum(sum((p_hat - p).^2)));
end
fprintf('reprojection error: %f\n', error);

end

