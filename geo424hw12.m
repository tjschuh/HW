function M = geo424hw12(n)
% M = GEO424HW12(n)
% This is a function to produce the moment
% tensors requested in HW 12 Problem 4
%
% INPUT:
%
% n      1   1999, Hector Mine Earthquake
%        2   1999, Chi-Chi, Taiwan Earthquake
%        3   2001, Bhuj, India Earthquake
%
% OUTPUT:
%
% M      3x3 moment tensor [N*m]
%
% TESTED ON: 9.8.0.1417392 (R2020a) Update 4
%
% Originally written by tschuh@princeton.edu, 5/4/2021

% depending on function input n, use certain dataset
switch n
    case 1
        disp('1999, Hector Mine Earthquake')
        % dip [degrees]
        delta = 85;
        % rake [degrees]
        lambda = 179;
        % strike [degrees]
        phi = 336;
        % seismic moment [N*m]
        M0 = 5.9e19;
    case 2
        disp('1999, Chi-Chi, Taiwan Earthquake')
        % dip [degrees]
        delta = 27;
        % rake [degrees]
        lambda = 82;
        % strike [degrees]
        phi = 26;
        % seismic moment [N*m]
        M0 = 4.1e20;
    case 3
        disp('2001, Bhuj, India Earthquake')
        % dip [degrees]
        delta = 50;
        % rake [degrees]
        lambda = 50;
        % strike [degrees]
        phi = 65;
        % seismic moment [N*m]
        M0 = 3.6e20;
    otherwise
        disp('Not an expected input. No calculation performed.')
        delta = NaN;
        lambda = NaN;
        phi = NaN;
        M0 = NaN;
end

% calculate the components of the moment tensor M
% M12=M21, M23=M32, M13=M31
M11 = 2*sind(delta)*sind(phi)*[cosd(lambda)*cosd(phi)-sind(lambda)*cosd(delta)*sind(phi)];

M12 = sind(delta)*sind(phi)*[cosd(lambda)*sind(phi)+sind(lambda)*cosd(delta)*cosd(phi)] ...
      - sind(delta)*cosd(phi)*[cosd(lambda)*cosd(phi)-sind(lambda)*cosd(delta)*sind(phi)];

M13 = sind(lambda)*sind(phi)*sind(delta)*sind(delta) ...
      + cosd(delta)*[cosd(lambda)*cosd(phi)-sind(lambda)*cosd(delta)*sind(phi)];

M22 = -2*sind(delta)*cosd(phi)*[cosd(lambda)*sind(phi)+sind(lambda)*cosd(delta)*cosd(phi)];

M23 = -sind(lambda)*cosd(phi)*sind(delta)*sind(delta) ...
      + cosd(delta)*[cosd(lambda)*sind(phi)+sind(lambda)*cosd(delta)*cosd(phi)];

M33 = 2*sind(lambda)*sind(delta)*cosd(delta);

% put all the elements together to make M the moment tensor
M = M0.*[M11 M12 M13; M12 M22 M23; M13 M23 M33];