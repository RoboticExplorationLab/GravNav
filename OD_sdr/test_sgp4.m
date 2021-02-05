clear


%% read TLE's in 
clear

tlecand =read_tle_data();

name_idx = 1:3:(length(tlecand)-2);
Nsats = length(name_idx);

TLEs.names = cell(length(name_idx),1);
TLEs.sats = cell(length(name_idx),1);
for i = 1:length(TLEs.names)
    TLEs.names{i} = strtrim(tlecand(name_idx(i)));
    TLEs.sats{i} = [tlecand(name_idx(i)+1),tlecand(name_idx(i)+2)];
end

%%

jd_init = 2.4592512800462963e6;
jd_final = jd_init + 20;
N = 10000;
traj.ECI = cell(Nsats,1);
traj.ECEF = cell(Nsats,1);
traj.JD = cell(Nsats,1);
for i = 1:1
    [traj.ECI{i},traj.JD{i}] = prop_orbit(TLEs.sats{i}(1),TLEs.sats{i}(2),jd_init,jd_final,N);
end


%% plotting
R_EARTH = 6378.1; %km
[x,y,z] = sphere(40);
imgRGB = imread('earth.jpg');
figure
hold on
% title('Molniya Orbit in ECEF')
warp(R_EARTH*x,R_EARTH*y,R_EARTH*z,circshift(rot90(imgRGB,2),569,2))
%plot3($gs1(1),$gs1(2),$gs1(3),'g.','MarkerSize',20)
%plot3($gsr(1),$gsr(2),$gsr(3),'b.','MarkerSize',20)
% plot3(eci_hist(1,:),eci_hist(2,:),eci_hist(3,:),'r','linewidth',2)
for i = 1:1
    eci_hist = traj.ECI{i};
    plot3(eci_hist(1,:),eci_hist(2,:),eci_hist(3,:),'r','linewidth',2)
end
% plot3($ecef_hist(1,:),eci_hist(2,:),eci_hist(3,:),'r','linewidth',2)
view(150,34)
xlabel('ECI X (m)')
ylabel('ECI Y (m)')
zlabel('ECI Z (m)')
axis equal
hold off


function ecef_hist = getecefhist(eci_hist,jd_hist)
    ecef_hist = zeros(3,size(eci_hist,2));
    for i = 1:length(jd_hist)
        ecef_hist(:,i) = getecef(eci_hist(1:3,i),jd_hist(i));
    end
end

function r_ecef =  getecef(r_eci,jd)
    r_ecef = dcmeci2ecef('IAU-2000/2006',datevec(datetime(jd,'convertfrom','juliandate')))*r_eci;
end

function [eci_hist, jd_vec] = prop_orbit(tle_line1,tle_line2,jd_init,jd_final,N)
    tle_line1 = convertStringsToChars(tle_line1);
    tle_line2 = convertStringsToChars(tle_line2);
    [~, ~, ~, satrec] = twoline2rv(tle_line1,tle_line2,'c','m','a',84);
    jd1 = satrec.jdsatepoch + satrec.jdsatepochf;
    jd_vec = linspace(jd_init,jd_final,N);
    t_vec_mins = (jd_vec - jd1) * 1440;
    eci_hist = zeros(6,N);
    for i = 1:N
    [satrec, eci_hist(1:3,i), eci_hist(4:6,i)] = sgp4(satrec, t_vec_mins(i));
    end
end