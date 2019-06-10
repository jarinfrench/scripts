function csl_boundary (_size, misorientation, _axis=100)

  if (mod(_size,2)==1)
    _size -= 1;
    printf('Changed size to %d',_size);
  endif

  c1 = cos(misorientation * pi/180);
  s1 = sin(misorientation * pi/180);

  % Generate the original lattice
  x1 = ones(_size+1).*[0:_size] - _size/2;
  y1 = transpose(x1);

  % Turn them into column vectors
  x1 = x1(:);
  y1 = y1(:);

  if (mod(_size*1.5,2)==1)
    modifier = 1;
  else
    modifier = 0;
  endif
  x2 = ones(_size*1.5 + 1).*[0:_size*1.5] - (_size*1.5 + modifier)/2;
  y2 = transpose(x2);
  x2 = x2(:);
  y2 = y2(:);

  if (_axis == 100)
    r = [cos(misorientation * pi/180.0) sin(misorientation * pi / 180.0);...
        -sin(misorientation * pi / 180.0) cos(misorientation * pi / 180.0)];

      % rotate the second matrix
    for i = 1:length(x2)
      point = r * [x2(i); y2(i)];
      if (point(1) < min(x1) || point(1) > max(x1) || ...
          point(2) < min(y1) || point(2) > max(y1))
        x2(i) = NaN;
        y2(i) = NaN;
      else
        x2(i) = point(1);
        y2(i) = point(2);
      end
    end
  elseif (_axis == 110)
    a_vec = [1 1 0]/norm([1 1 0]);

    % Create the cube for lattice 1
    x1_2 = repmat(x1,length(min(x1):max(x1)),1);
    y1_2 = repmat(y1,length(min(x1):max(x1)),1);
    z1_2 = min(x1_2):max(x1_2);
    z1_2 = repelems(z1_2, [1:length(z1_2); (length(z1_2)^2)*ones(size(z1_2))])';

    % Create the cube for lattice 2
    x2_2 = repmat(x2, length(min(x2):max(x2)), 1);
    y2_2 = repmat(y2, length(min(x2):max(x2)), 1);
    z2_2 = min(x2_2):max(x2_2);
    z2_2 = repelems(z2_2, [1:length(z2_2); (length(z2_2)^2)*ones(size(z2_2))])';

    % This rotates about the [110] axis by misorientation
    r = [c1+a_vec(1)^2*(1-c1)                  a_vec(1)*a_vec(2)*(1-c1)-a_vec(3)*s1   a_vec(1)*a_vec(3)*(1-c1)+a_vec(2)*s1
         a_vec(1)*a_vec(2)*(1-c1)+a_vec(3)*s1  c1+a_vec(2)^2*(1-c1)                   a_vec(2)*a_vec(3)*(1-c1)-a_vec(1)*s1
         a_vec(1)*a_vec(3)*(1-c1)-a_vec(2)*s1  a_vec(2)*a_vec(3)*(1-c1)+a_vec(1)*s1   c1+a_vec(3)^2*(1-c1)];

    tmp2 = transpose([x2_2 y2_2 z2_2]);

    rot_lat2 = r*tmp2;

    tmp_x = rot_lat2(1,:);
    tmp_y = rot_lat2(2,:);
    tmp_z = rot_lat2(3,:);

    x1_updated = [];
    y1_updated = [];
    z1_updated = [];
    x2_updated = [];
    y2_updated = [];
    z2_updated = [];

    % The (110) plane is defined by the equation -x -y = 0, so:
    for i=1:length(tmp_z)
      if (abs(-tmp_x(i) - tmp_y(i)) < 1e-8)
        x2_updated = [x2_updated; tmp_x(i)];
        y2_updated = [y2_updated; tmp_y(i)];
        z2_updated = [z2_updated; tmp_z(i)];
      end
      if (i > length(x1_2))
        continue
      else
        if (abs(-x1_2(i) - y1_2(i)) < 1e-8)
          x1_updated = [x1_updated; x1_2(i)];
          y1_updated = [y1_updated; y1_2(i)];
          z1_updated = [z1_updated; z1_2(i)];
        end
      end
    end

    x1 = x1_updated;
    y1 = y1_updated;
    z1 = z1_updated;

    x2 = x2_updated;
    y2 = y2_updated;
    z2 = z2_updated;

    t = [0 0 1];
    v = cross(a_vec,t);
    u = v / norm(v);
    c = dot(a_vec, t);
    h = (1 - c) / dot(v, v);

    rot2z = [c + h*v(1)^2      h*v(1)*v(2)-v(3)  h*v(1)*v(3)+v(2);
             h*v(1)*v(2)+v(3)  c + h*v(2)^2      h*v(2)*v(3)-v(1);
             h*v(1)*v(3)-v(2)  h*v(2)*v(3)+v(1)  c+h*v(3)^2];

    single_plane_1 = rot2z * transpose([x1,y1,z1]);
    single_plane_2 = rot2z * transpose([x2,y2,z2]);

    x1 = single_plane_1(1,:)';
    y1 = single_plane_1(2,:)';

    x2 = single_plane_2(1,:)';
    y2 = single_plane_2(2,:)';

    for i = 1:length(x2)
      if (x2(i) < min(x1) || x2(i) > max(x1) || y2(i) < min(y1) || y2(i) > max(y1))
        x2(i) = NaN;
        y2(i) = NaN;
      end
    end
  elseif (_axis == 111)
    a_vec = [1 1 1]/norm([1 1 1]);

    % Create the cube for lattice 1
    x1_2 = repmat(x1,length(min(x1):max(x1)),1);
    y1_2 = repmat(y1,length(min(x1):max(x1)),1);
    z1_2 = min(x1_2):max(x1_2);
    z1_2 = repelems(z1_2, [1:length(z1_2); (length(z1_2)^2)*ones(size(z1_2))])';

    % Create the cube for lattice 2
    x2_2 = repmat(x2, length(min(x2):max(x2)), 1);
    y2_2 = repmat(y2, length(min(x2):max(x2)), 1);
    z2_2 = min(x2_2):max(x2_2);
    z2_2 = repelems(z2_2, [1:length(z2_2); (length(z2_2)^2)*ones(size(z2_2))])';

    % This rotates about the [111] axis by misorientation
    r = [c1+a_vec(1)^2*(1-c1)                  a_vec(1)*a_vec(2)*(1-c1)-a_vec(3)*s1   a_vec(1)*a_vec(3)*(1-c1)+a_vec(2)*s1
         a_vec(1)*a_vec(2)*(1-c1)+a_vec(3)*s1  c1+a_vec(2)^2*(1-c1)                   a_vec(2)*a_vec(3)*(1-c1)-a_vec(1)*s1
         a_vec(1)*a_vec(3)*(1-c1)-a_vec(2)*s1  a_vec(2)*a_vec(3)*(1-c1)+a_vec(1)*s1   c1+a_vec(3)^2*(1-c1)];

    tmp2 = transpose([x2_2 y2_2 z2_2]);

    rot_lat2 = r*tmp2;

    tmp_x = rot_lat2(1,:);
    tmp_y = rot_lat2(2,:);
    tmp_z = rot_lat2(3,:);

    x1_updated = [];
    y1_updated = [];
    z1_updated = [];
    x2_updated = [];
    y2_updated = [];
    z2_updated = [];

    % The (111) plane is defined by the equation x + y + z = 0, so:
    for i=1:length(tmp_z)
      if (abs(tmp_x(i) + tmp_y(i) + tmp_z(i)) < 1e-8)
        x2_updated = [x2_updated; tmp_x(i)];
        y2_updated = [y2_updated; tmp_y(i)];
        z2_updated = [z2_updated; tmp_z(i)];
      end
      if (i > length(x1_2))
        continue
      else
        if (abs(x1_2(i) + y1_2(i) + z1_2(i)) < 1e-8)
          x1_updated = [x1_updated; x1_2(i)];
          y1_updated = [y1_updated; y1_2(i)];
          z1_updated = [z1_updated; z1_2(i)];
        end
      end
    end

    x1 = x1_updated;
    y1 = y1_updated;
    z1 = z1_updated;

    x2 = x2_updated;
    y2 = y2_updated;
    z2 = z2_updated;

    % Now lets rotate everything so [111] is aligned with [001]
    t = [0 0 1];
    v = cross(a_vec,t);
    u = v / norm(v);
    c = dot(a_vec, t);
    h = (1 - c) / dot(v, v);

    rot2z = [c + h*v(1)^2      h*v(1)*v(2)-v(3)  h*v(1)*v(3)+v(2);
             h*v(1)*v(2)+v(3)  c + h*v(2)^2      h*v(2)*v(3)-v(1);
             h*v(1)*v(3)-v(2)  h*v(2)*v(3)+v(1)  c+h*v(3)^2];

    single_plane_1 = rot2z * transpose([x1,y1,z1]);
    single_plane_2 = rot2z * transpose([x2,y2,z2]);

    x1 = single_plane_1(1,:)';
    y1 = single_plane_1(2,:)';

    x2 = single_plane_2(1,:)';
    y2 = single_plane_2(2,:)';

    for i = 1:length(x2)
      if (x2(i) < min(x1) || x2(i) > max(x1) || y2(i) < min(y1) || y2(i) > max(y1))
        x2(i) = NaN;
        y2(i) = NaN;
      end
    end
  else
    disp('This can only do the high symmetry rotation axes (100, 110, or 111)');
    return;
  end

  % create the original point matrix
  p1 = [x1 y1];

  x2(isnan(x2)) = [];
  y2(isnan(y2)) = [];

  p2 = [x2 y2];


  all_lattice_points = [p1; p2];
  p3 = NaN(size(p1));

  % generate a cell-linked list between the two point lists
  lx = max(x1) - min(x1);
  ly = max(y1) - min(y1);
  rcut = 1.5;

  ncellx = floor(lx/rcut) + 1; % +1 so we don't divide by 0
  ncelly = floor(ly/rcut) + 1;

  lcellx = lx/ncellx;
  lcelly = ly/ncelly;

  n_per_cell = 100;

  icell = zeros(ncellx,ncelly);
  pcell = zeros(ncellx,ncelly,n_per_cell);
  iatom = zeros(n_per_cell,length(all_lattice_points));

  for i = 1:length(all_lattice_points)
    id = floor((all_lattice_points(i,:) - [min(x1) min(y1)]./[lcellx lcelly]))+1;
    if (id(1) > ncellx) id(1) = ncellx; end
    if (id(2) > ncelly) id(2) = ncelly; end
    if (id(1) < 1) id(1) = 1; end
    if (id(2) < 1) id(2) = 1; end
    ++icell(id(1),id(2));
    pcell(id(1),id(2),icell(id(1),id(2))) = i;
  end

  for i = 1:ncellx
    for j = 1:ncelly
      for l = 1:icell(i,j)
        id = pcell(i,j,l);

        % check surrounding cells
        for ii = -1:1
          for jj = -1:1
            ia = i + ii; % min value: 0.  Max value: ncellx+1
            ja = j + jj;

            if (ia > ncellx) ia = 1; end
            if (ja > ncelly) ja = 1; end
            if (ia < 1) ia = ncellx; end
            if (ja < 1) ja = ncelly; end

            for m = 1:icell(ia,ja)
              jd = pcell(ia,ja,m);
              if (jd <= id) continue; end

              rxij = all_lattice_points(id,1) - all_lattice_points(jd,1);
              ryij = all_lattice_points(id,2) - all_lattice_points(jd,2);

              rxij = rxij - round(rxij / lx) * lx;
              ryij = ryij - round(ryij / ly) * ly;

              drij_sq = rxij^2 + ryij^2;
              if (drij_sq > rcut^2) continue; end

              ++iatom(1,id); % number of neighbors increases
              iatom(iatom(1,id),id) = jd; % give the neighbor id
              ++iatom(1,jd); % similarly for this atom
              iatom(iatom(1,jd),jd) = id;
            end
          end
        end
      end
    end
  end

  clear pcell icell;

  for i = 1:length(all_lattice_points)
    x = all_lattice_points(i,1);
    y = all_lattice_points(i,2);
    for l = 2:iatom(1,i)
      id = iatom(l,i);
      if (id <= length(p1)) % we don't want to include the original lattice
        continue;
      end
      rxij = all_lattice_points(id,1) - x;
      ryij = all_lattice_points(id,2) - y;

      rxij = rxij - round(rxij / lx) * lx;
      ryij = ryij - round(ryij / ly) * ly;

      drij_sq = rxij^2 + ryij^2;
      if (drij_sq < 1.0e-5)

        p3(i,1) = all_lattice_points(i,1);
        p3(i,2) = all_lattice_points(i,2);
      end
    end
  end

  p3(any(isnan(p3),2),:) = [];
  p3 = unique(p3, 'rows');

  figure();
  xvals = {x1, x2, p3(:,1), 0};
  yvals = {y1, y2, p3(:,2), 0};
  sizes = {1, 1 , 2, 3};
  colors = {[1 0 0], [0 0 1], [0 1 0], [0 0 0]};
  styles = {'o','o','o','o'};
  scatter_series_set(xvals, yvals, sizes, colors, styles);
  axis([min(x1)-5 max(x2)+5 min(y1)-5 max(y1)+5], 'square');


  % scatter(x1,y1,25,"r","filled")
  % hold on
  % scatter(x2,y2,25,"b")
  % scatter(p3(:,1),p3(:,2),36,"g","filled")
  % scatter(0, 0, 49, "k", "filled");

  legend('Original lattice', 'Rotated lattice', 'CSL Points', 'Origin', 'location', 'northeast')
end
