function [res, eps] = reinit(phi, dt, h)
  % Solve reinitialization equation for the Signed Diatance Function
  % phi_t + sign(phi)(1 - |grad(phi)|)=0, which gives SDF with property |grad(phi)|=1
  % Fast marching method from http://persson.berkeley.edu/pub/persson05levelset.pdf, formula (16)
  % forward/backward differences
  phi_x_b = (phi - backward_x(phi)) / h(1);
  phi_x_f = (forward_x(phi) - phi) / h(1);
  phi_y_b = (phi - backward_y(phi)) / h(2);
  phi_y_f = (forward_y(phi) - phi) / h(2);
  phi_z_b = (phi - backward_z(phi)) / h(3);
  phi_z_f = (forward_z(phi) - phi) / h(3);
  
  x_b_p = phi_x_b;  x_b_n = phi_x_b;
  x_f_p = phi_x_f;  x_f_n = phi_x_f;
  y_b_p = phi_y_b;  y_b_n = phi_y_b; 
  y_f_p = phi_y_f;  y_f_n = phi_y_f;
  z_b_p = phi_z_b;  z_b_n = phi_z_b;
  z_f_p = phi_z_f;  z_f_n = phi_z_f;
  
  x_b_p(phi_x_b < 0) = 0;
  x_b_n(phi_x_b > 0) = 0;
  x_f_p(phi_x_f < 0) = 0;
  x_f_n(phi_x_f > 0) = 0;
  
  y_b_p(phi_y_b < 0) = 0;
  y_b_n(phi_y_b > 0) = 0;
  y_f_p(phi_y_f < 0) = 0;
  y_f_n(phi_y_f > 0) = 0;
  
  z_b_p(phi_z_f < 0) = 0;
  z_b_n(phi_z_f > 0) = 0;
  z_f_p(phi_z_f < 0) = 0;
  z_f_n(phi_z_f > 0) = 0;
  
  df = zeros(size(phi));
  phi_neg_ind = find(phi < 0);
  phi_pos_ind = find(phi > 0);
  df(phi_pos_ind) = 1 - sqrt(max(y_b_p(phi_pos_ind).^2, y_f_n(phi_pos_ind).^2) ...
                       + max(x_b_p(phi_pos_ind).^2, x_f_n(phi_pos_ind).^2) ...
                       + max(z_b_p(phi_pos_ind).^2, z_f_n(phi_pos_ind).^2));
 
  df(phi_neg_ind) = 1 - sqrt(max(y_b_n(phi_neg_ind).^2, y_f_p(phi_neg_ind).^2) ...
                       + max(x_b_n(phi_neg_ind).^2, x_f_p(phi_neg_ind).^2) ...
                       + max(z_b_n(phi_neg_ind).^2, z_f_p(phi_neg_ind).^2));

  delta = dt .* sign(phi) .* df;

  res = phi + delta;
  eps = min(min(min(abs(df)))); % should be |grad(phi)| = 1
end

function [shift] = forward_x(M)
  shift = cat(1, M(2:end, :, :), M(end, :, :));
end

function [shift] = backward_x(M)
  shift = cat(1, M(1, :, :), M(1:end-1, :, :));
end

function [shift] = forward_y(M)
  shift = cat(2, M(:, 2:end, :), M(:, end, :));
end

function [shift] = backward_y(M)
  shift = cat(2, M(:, 1, :), M(:, 1:end-1, :));
end

function [shift] = forward_z(M)
  shift = cat(3, M(:, :, 2:end), M(:, :, end));
end

function [shift] = backward_z(M)
  shift = cat(3, M(:, :, 1), M(:, :, 1:end-1));
end

function [sgn] = sign(phi)
  sgn = phi ./ sqrt(phi.^2 + 1);
end
