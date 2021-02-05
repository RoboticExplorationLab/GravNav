clear


%% read TLE's in 
clear

tlecand =read_tle_data();

name_idx = 1:3:(length(tlecand)-2);


TLEs.names = cell(length(name_idx),1);
TLEs.sats = cell(length(name_idx),1);
for i = 1:length(TLEs.names)
    TLEs.names{i} = strtrim(tlecand(name_idx(i)));
    TLEs.sats{i} = [tlecand(name_idx(i)+1),tlecand(name_idx(i)+2)];
end


%%

R_EARTH = 6378.1; %km


tle.l1 = '1 47463U 21006BC  21035.45745765  .00001106  00000-0  69163-4 0  9994';
tle.l2 = '2 47463  97.5023  98.9325 0008900 214.3497 145.7159 15.11436064  1942';

test_jd = 2459251.28005;
eci_hist = prop_orbit(tle,test_jd,test_jd + 95/1440,100);


[x,y,z] = sphere(40);
imgRGB = imread('earth.jpg');
figure
hold on
% title('Molniya Orbit in ECEF')
warp(R_EARTH*x,R_EARTH*y,R_EARTH*z,circshift(rot90(imgRGB,2),569,2))
%plot3($gs1(1),$gs1(2),$gs1(3),'g.','MarkerSize',20)
%plot3($gsr(1),$gsr(2),$gsr(3),'b.','MarkerSize',20)
plot3(eci_hist(1,:),eci_hist(2,:),eci_hist(3,:),'r','linewidth',2)
% plot3($ecef_hist(1,:),eci_hist(2,:),eci_hist(3,:),'r','linewidth',2)
view(150,34)
xlabel('ECI X (m)')
ylabel('ECI Y (m)')
zlabel('ECI Z (m)')
axis equal
hold off


function r_ecef =  getecef(r_eci,jd)
r_ecef = dcmeci2ecef('IAU-2000/2006',datevec(datetime(jd,'convertfrom','juliandate')))*r_eci;
end

function eci_hist = prop_orbit(tle,jd_init,jd_final,N)
[~, ~, ~, satrec] = twoline2rv(tle.l1,tle.l2,'c','m','a',84);
jd1 = satrec.jdsatepoch + satrec.jdsatepochf;
jd_vec = linspace(jd_init,jd_final,N);
t_vec_mins = (jd_vec - jd1) * 1440;
eci_hist = zeros(6,N);
for i = 1:N
[satrec, eci_hist(1:3,i), eci_hist(4:6,i)] = sgp4(satrec, t_vec_mins(i));
end
end