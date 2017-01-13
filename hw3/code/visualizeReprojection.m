function [ ] = visualizeReprojection( graph, frames )

K = f2K(graph.f);
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
    
    im = imread(frames.images{graph.frames(i)});
    % uncenter feature points
    p(1,:) = size(im,2)/2 - p(1,:);
    p(2,:) = size(im,1)/2 - p(2,:);
    
    p_hat(1,:) = size(im,2)/2 - p_hat(1,:);
    p_hat(2,:) = size(im,1)/2 - p_hat(2,:);
    
    figure;
    imshow(im); hold on;
    scatter(p(1,:), p(2,:), 'rx')
    scatter(p_hat(1,:), p_hat(2,:), 'g+')
    plot([p(1,:); p_hat(1,:)], [p(2,:); p_hat(2,:)], 'b');
    
    % also display unobserved points
    P_u = graph.Str(:, ~str_idx);
    p_u = K*M*[P_u; ones(1, size(P_u,2))];
    p_u = [p_u(1,:)./p_u(3,:); p_u(2,:)./p_u(3,:)];
    p_u(1,:) = size(im,2)/2 - p_u(1,:);
    p_u(2,:) = size(im,1)/2 - p_u(2,:);
    scatter(p_u(1,:), p_u(2,:), 'yo')
    
    
end

end

