function [rectI, H] = rectifyImage(filename, debug)

close all

im = imread(filename);

if debug == 0
  [d, fname] = fileparts(filename);
  load(fullfile(strrep(d, '/images', '/data'), [fname '.mat']));
   
else
  imshow(im);
  
  %get parallel lines
  [x, y] = ginput(8);
  
  
  p1 = [x(1) y(1) 1]';
  p2 = [x(2) y(2) 1]';
  p3 = [x(3) y(3) 1]';
  p4 = [x(4) y(4) 1]';
  p5 = [x(5) y(5) 1]';
  p6 = [x(6) y(6) 1]';
  p7 = [x(7) y(7) 1]';
  p8 = [x(8) y(8) 1]';
  
  l1 = cross(p1, p2);
  %p1'*l1
  %p2'*l1

  l2 = cross(p3, p4);
  %p3'*l2
  %p3'*l2

  l3 = cross(p5, p6);
  l4 = cross(p7, p8);
  
  vp1 = cross(l1, l2);
  vp2 = cross(l3, l4);

  l_inf = cross(vp1, vp2);
  l_inf = l_inf ./ l_inf(3);
  
  %show annotations
  hold on;
  plot(x(1:2), y(1:2), 'Marker', 'x')
  plot(x(3:4), y(3:4), 'Marker', 'x', 'Color', 'b')
  plot(x(5:6), y(5:6), 'Marker', 'x', 'Color', 'r')
  plot(x(7:8), y(7:8), 'Marker', 'x', 'Color', 'r')
  
  Ha = [1 0 0; 0 1 0; l_inf'];

  %inv(Ha)'*l_inf

  affine_im = applyH(im, Ha);
  figure;
  imshow(affine_im)
  
  %get orthogonal lines
  %last 4 just to check cosine
  [ax, ay] = ginput(12);

  pa1 = [ax(1) ay(1) 1]';
  pa2 = [ax(2) ay(2) 1]';
  pa3 = [ax(3) ay(3) 1]';
  pa4 = [ax(4) ay(4) 1]';
  pa5 = [ax(5) ay(5) 1]';
  pa6 = [ax(6) ay(6) 1]';
  pa7 = [ax(7) ay(7) 1]';
  pa8 = [ax(8) ay(8) 1]';
  pa9 = [ax(9) ay(9) 1]';
  pa10 = [ax(10) ay(10) 1]';
  pa11 = [ax(11) ay(11) 1]';
  pa12 = [ax(12) ay(12) 1]';

  la1 = cross(pa1, pa2);
  la2 = cross(pa3, pa4);
  la3 = cross(pa5, pa6);
  la4 = cross(pa7, pa8);
  la5 = cross(pa9, pa10);
  la6 = cross(pa11, pa12);
  
  %show annotations
  hold on;
  plot(ax(1:2), ay(1:2), 'Marker', 'x')
  plot(ax(3:4), ay(3:4), 'Marker', 'x', 'Color', 'b')
  plot(ax(5:6), ay(5:6), 'Marker', 'x', 'Color', 'r')
  plot(ax(7:8), ay(7:8), 'Marker', 'x', 'Color', 'r')
  
  [~, fname] = fileparts(filename);
  save(['data/', fname, '.mat'], 'x', 'y', 'ax', 'ay', 'Ha',...
        'l1', 'l2', 'l3', 'l4', 'la1', 'la2', 'la3', 'la4', 'la5', 'la6');
  
end

  A = [la1(1)*la2(1), la1(1)*la2(2)+la1(2)*la2(1), la1(2)*la2(2);
       la3(1)*la4(1), la3(1)*la4(2)+la3(2)*la4(1), la3(2)*la4(2)];
       
  s = null(A);

  S = [s(1) s(2); s(2) s(3)];

  C_inf = [S(1,:) 0; S(2,:) 0; 0 0 0];
%  la1'*C_inf*la2
%  la3'*C_inf*la4

  try
    % use lower since we assume S = KK',
    % and thats how matlab defines lower
    K = chol(S, 'lower');
  catch
    fprintf('cholesky failed!, trying SVD\n');
    % alternate method
    [U, D, V] = svd(S);
    K = U*sqrt(D)*V';
  end
  
  %[la3(1:2)]'*(K*K')*[la4(1:2)]

  Hs = [K(1,:) 0; K(2,:) 0; 0 0 1];

  H = inv(Hs)*Ha;
  rectI = applyH(im, H);
  
  c_inf_o = [1 0 0; 0 1 0; 0 0 0];
 before = [cosine(Ha'*la1, Ha'*la2, c_inf_o);
           cosine(Ha'*la3, Ha'*la4, c_inf_o);
           cosine(Ha'*la5, Ha'*la6, c_inf_o)]
 after = [cosine(Hs'*la1, Hs'*la2, c_inf_o);
          cosine(Hs'*la3, Hs'*la4, c_inf_o);
          cosine(Hs'*la5, Hs'*la6, c_inf_o)]
           
  if debug == 0
    %show annotated images
    figure; subplot(1,3,1); imshow(im);
    hold on;
    plot(x(1:2), y(1:2), 'Marker', 'x')
    plot(x(3:4), y(3:4), 'Marker', 'x', 'Color', 'b')
    plot(x(5:6), y(5:6), 'Marker', 'x', 'Color', 'r')
    plot(x(7:8), y(7:8), 'Marker', 'x', 'Color', 'r')
    
    affine_im = applyH(im, Ha);
    subplot(1,3,2); imshow(affine_im);
    hold on;
    plot(ax(1:2), ay(1:2), 'Marker', 'x')
    plot(ax(3:4), ay(3:4), 'Marker', 'x', 'Color', 'b')
    plot(ax(5:6), ay(5:6), 'Marker', 'x', 'Color', 'r')
    plot(ax(7:8), ay(7:8), 'Marker', 'x', 'Color', 'r')
    
    subplot(1,3,3); imshow(rectI);
  end
  
  %figure; imshow(rectI);
 

end