
% ray tracing

%%%%%%%%%%%%%%%%%%%%%%%%%%
% initiate all variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% the number of rays
number_of_rays = 200;

% set geometry of the mirror

% normal vector of the plane of the mirror (x,y,z)
mirror_normal = [0.2, 0, 0.7071];

% point on the plane of the mirror(x,y,z)
mirror_point = [0, 0, 1];

% extent of the mirror
mirror_min_x = 0.25;
mirror_max_x = 0.75;
mirror_min_y = 0.25;
mirror_max_y = 0.75;

% define an aperture through which all relevant rays will pass
aperture_min_x = 0;
aperture_max_x = 1;
aperture_min_y = 0;
aperture_max_y = 1;
aperture_z = 2;

% initiate vectors and matrices for the rays

% point on each incident ray located on the aperture
% there are many rays, so this is a 2 dimensional matrix
% the first index indentifies the ray
% the second index provides the spatial dimension (x,y,z)
% using the monte-carlo method, so pick random points
ray_incident_point = zeros(number_of_rays,3);
ray_incident_point(:,1) = aperture_min_x + rand(number_of_rays,1)*(aperture_max_x - aperture_min_x); % x
ray_incident_point(:,2) = aperture_min_y + rand(number_of_rays,1)*(aperture_max_y - aperture_min_y); % y
ray_incident_point(:,3) = aperture_z; % z

% direction of the incident ray (x,y,z)
% azimuth (compass direction)of each ray will be random, and evenly
% distributed, in rad
% the zenith is the angle from the vertical in mrad
ray_incident_direction = zeros(number_of_rays,3);
ray_intensity = zeros(number_of_rays,1);
% intensity data with CSR = 10
% taken from Neumann et al, Transactions of the ASME, p198, v124, 2002.
intensity_data = [1,0.99,0.96,0.90,0.79,0.049,0.025,0.016,0.009,0.007,0.005];
for ray = 1:number_of_rays
    ray_incident_azimuth = 2*pi*rand; % in rad
    ray_incident_zenith = 10*rand; % in mrad
    
    ray_incident_unit_radius = tan(ray_incident_zenith/1000);
    [ray_incident_direction(ray,1),ray_incident_direction(ray,2)] = pol2cart(ray_incident_azimuth,ray_incident_unit_radius);
    ray_incident_direction(ray,3) = -1; % z
    
    % associate a relative intensity with each ray
    ray_intensity(ray,1) = intensity_data(ceil(ray_incident_zenith));

end

% calculate the power scaling factor


% does the ray intersect the mirror? true(1) or false(0)
% initiate to false(0)
mirror_intersect = zeros(number_of_rays,1);

% point of intersection between ray and mirror plane
% intiate but calculated later
ray_intersect_point = zeros(number_of_rays,3);

% point on reflected ray
ray_reflected_point = zeros(number_of_rays,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform the ray tracing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ray = 1:number_of_rays
    
    % find intersection between ray and mirror plane
    m = dot((mirror_point - ray_incident_point(ray,:)),mirror_normal) / dot(ray_incident_direction(ray,:),mirror_normal);
    ray_intersect_point(ray,:) = ray_incident_point(ray,:) + m*ray_incident_direction(ray,:);

    % find the reflected point (on the other side of the mirror)
    q = dot((mirror_point - ray_incident_point(ray,:)),mirror_normal) / dot(mirror_normal,mirror_normal);
    ray_reflected_point_b = ray_incident_point(ray,:) + 2*q*mirror_normal;

    % find the reflected point (on the original side of the mirror)
    ray_reflected_point(ray,:) = 2*ray_intersect_point(ray,:) - ray_reflected_point_b;

    % calculate if the intersection is on the mirror
    % and set mirror_intersect to true or false
    if (ray_intersect_point(ray,1) > mirror_min_x) && (ray_intersect_point(ray,1) < mirror_max_x) && (ray_intersect_point(ray,2) > mirror_min_y) && (ray_intersect_point(ray,2) < mirror_max_y)
        mirror_intersect(ray) = true;
    else
        mirror_intersect(ray) = false;
    end
    
end
    
% the Matlab line function requires separate vectors for x, y, and z
% set these up
x = zeros(1,2);
y = zeros(1,2);
z = zeros(1,2);

figure(1);
axis equal;

% plot the incident rays
for ray = 1:number_of_rays
    x(1) = ray_incident_point(ray,1);
    y(1) = ray_incident_point(ray,2);
    z(1) = ray_incident_point(ray,3);
    x(2) = ray_intersect_point(ray,1);
    y(2) = ray_intersect_point(ray,2);
    z(2) = ray_intersect_point(ray,3);
    % plot a green line
    line(x,y,z,'color','g');
end

% plot the reflected rays
for ray = 1:number_of_rays
    if(mirror_intersect(ray) == true)
        x(1) = ray_intersect_point(ray,1);
        y(1) = ray_intersect_point(ray,2);
        z(1) = ray_intersect_point(ray,3);
        x(2) = ray_reflected_point(ray,1);
        y(2) = ray_reflected_point(ray,2);
        z(2) = ray_reflected_point(ray,3);
        % plot a blue line
        line(x,y,z,'color','b');
    end
  %  if(mirror_intersect(ray) == false)
  %      x(1) = ray_intersect_point(ray,1);
  %      y(1) = ray_intersect_point(ray,2);
  %      z(1) = ray_intersect_point(ray,3);
  %      x(2) = ray_intersect_point(ray,1);
  %      y(2) = ray_intersect_point(ray,2);
   %     z(2) = 0;
        % plot a green line
    %    line(x,y,z,'color','g');
  %  end
end

% plot normal to the mirror
x(1) = mirror_point(1);
y(1) = mirror_point(2);
z(1) = mirror_point(3);
x(2) = mirror_point(1) + mirror_normal(1);
y(2) = mirror_point(2) + mirror_normal(2);
z(2) = mirror_point(3) + mirror_normal(3);
% plot a red line
line(x,y,z,'color','r');


