function csl_boundary (_size, misorientation, _axis=100)
  if (_axis ~= 100)
    disp('Different axes not yet implemented');
    return
  end
  x1 = ones(_size+1).*[0:_size] - _size/2;
  y1 = transpose(x1);
  rotation_matrix = [cos(misorientation * pi/180.0) sin(misorientation * pi / 180.0); -sin(misorientation * pi / 180.0) cos(misorientation * pi / 180.0)];

  % Turn them into column vectors
  x1 = x1(:);
  y1 = y1(:);

  % create the point matrix
  p1 = [x1 y1];

  % Likewise for the second point matrix
  x2 = ones(_size*2 + 1).*[0:_size*2] - _size;
  y2 = transpose(x2);
  x2 = x2(:);
  y2 = y2(:);

  % Specify the end points of the lattice lines for the first point matrix
  horz_x1 = [min(x1)*ones(size(unique(y1))) max(x1) * ones(size(unique(y1)))];
  horz_y1 = [unique(y1) unique(y1)];
  vert_x1 = [unique(x1) unique(x1)];
  vert_y1 = [min(y1) * ones(size(unique(x1))) max(y1) * ones(size(unique(x1)))];

  lines_x1 = [horz_x1; vert_x1];
  lines_y1 = [horz_y1; vert_y1];

  % Likewise for the second point matrix
  horz_x2 = [min(x2)*ones(size(unique(y2))) max(x2) * ones(size(unique(y2)))];
  horz_y2 = [unique(y2) unique(y2)];
  vert_x2 = [unique(x2) unique(x2)];
  vert_y2 = [min(y2) * ones(size(unique(x2))) max(y2) * ones(size(unique(x2)))];

  lines_x2 = [horz_x2; vert_x2];
  lines_y2 = [horz_y2; vert_y2];

  tmp = (rotation_matrix * [lines_x2(:) lines_y2(:)]')'; % rotate the second lattice
  lines_x2 = reshape(tmp(:,1), length(tmp)/2, 2);
  lines_y2 = reshape(tmp(:,2), length(tmp)/2, 2);


  % rotate the second matrix
  for i = 1:length(x2)
    point = rotation_matrix * [x2(i); y2(i)];
    if (point(1) < min(x1) || point(1) > max(x1) || ...
        point(2) < min(y1) || point(2) > max(y1))
      x2(i) = NaN;
      y2(i) = NaN;
    else
      x2(i) = point(1);
      y2(i) = point(2);
    end
  end

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
      if (drij_sq < 1.0e-8)

        p3(i,1) = all_lattice_points(i,1);
        p3(i,2) = all_lattice_points(i,2);
      end
    end
  end

  p3(any(isnan(p3),2),:) = [];
  p3 = unique(p3, 'rows');

  figure();
  scatter(x1,y1,25,"r","filled")
  axis([min(x1) max(x2) min(y1) max(y1)], "square");
  hold on
  scatter(x2,y2,25,"b")
  scatter(p3(:,1),p3(:,2),36,"g","filled")
  scatter(0, 0, 49, "k", "filled");

  line(lines_x1', lines_y1', "linestyle", "-", "color", "r");
  line(lines_x2', lines_y2', "linestyle", "-", "color", "b");




end
