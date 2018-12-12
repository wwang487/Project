% LUV2RGB performs Luv to RGB color conversion.
% This code is downloaded and modified from CS766 website.
%   RGB = LUV2RGB( LUV ) converts the LUV color vectors into RGB.
%
% See also RGB2LUV.
%

% 6.891 Problem Set 4: Problem 1
% C. Mario Christoudias

function rgb_img = luvtorgb(luv_img)
    height = size(luv_img,1);
    width = size(luv_img,2);
    temp_img = zeros(height,width,3);
    temp_vec = zeros(3,1);
    for i = 1:height
        for j = 1:width
            temp_vec(1,1) = luv_img(i,j,1);
            temp_vec(2,1) = luv_img(i,j,2);
            temp_vec(3,1) = luv_img(i,j,3);
            temp_rgb = vec_luv2rgb(temp_vec);
            temp_img(i,j,1) = temp_rgb(1,1);
            temp_img(i,j,2) = temp_rgb(2,1);
            temp_img(i,j,3) = temp_rgb(3,1);
        end
    end 
    rgb_img = temp_img;
end

function rgb = vec_luv2rgb( luv )

% convert to xyz

xyz = luv;

Yn = 1;
Un_prime = 0.19784977571475;
Vn_prime = 0.46834507665248;

u_prime	= luv(2,:) ./ (13 * luv(1,:)) + Un_prime;
v_prime	= luv(3,:) ./ (13 * luv(1,:)) + Vn_prime;

xyz(2,:) = (luv(1,:) < 8.0) .* (Yn * luv(1,:) / 903.3) + (luv(1,:) >= 8.0) .* (Yn * (((luv(1,:) + 16.0) / 116.0) .^ 3));
xyz(1,:) = 9 * u_prime .* xyz(2,:) ./ (4 * v_prime);
xyz(3,:) = (12 - 3 * u_prime - 20 * v_prime) .* xyz(2,:) ./ (4 * v_prime);

% convert to RGB

RGB = [3.2405, -1.5371, -0.4985; ...
       -0.9693,  1.8760,  0.0416; ...
       0.0556, -0.2040,  1.0573];

rgb = RGB * xyz;

% bounds check

rgb(1,:) = (luv(1,:) >= 0.1) .* rgb(1,:);
rgb(2,:) = (luv(1,:) >= 0.1) .* rgb(2,:);
rgb(3,:) = (luv(1,:) >= 0.1) .* rgb(3,:);

rgb(1,:) = (rgb(1,:) >= 0 & rgb(1,:) <= 1) .* rgb(1,:) + (rgb(1,:) > 1);
rgb(2,:) = (rgb(2,:) >= 0 & rgb(2,:) <= 1) .* rgb(2,:) + (rgb(2,:) > 1);
rgb(3,:) = (rgb(3,:) >= 0 & rgb(3,:) <= 1) .* rgb(3,:) + (rgb(3,:) > 1);
end