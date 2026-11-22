// Phase portraits of relative motion using the Hill-Clohessy-Wiltshire model.
// Independent circular reference orbit, several separation scales.

clear;
clc;

function c = cross3(a, b)
    c = [a(2)*b(3)-a(3)*b(2);
        a(3)*b(1)-a(1)*b(3);
        a(1)*b(2)-a(2)*b(1)];
endfunction

function [r, v] = keplerState(a, e, inc, raan, argp, nu, mu)
    // Position and velocity in the inertial frame from Kepler elements.
    p = a*(1-e^2);
    radius = p/(1+e*cos(nu));
    u = argp+nu;

    radialSpeed = sqrt(mu/p)*e*sin(nu);
    transverseSpeed = sqrt(mu/p)*(1+e*cos(nu));

    cO = cos(raan);
    sO = sin(raan);
    cu = cos(u);
    su = sin(u);
    ci = cos(inc);
    si = sin(inc);

    r = radius*[cO*cu-sO*su*ci;
        sO*cu+cO*su*ci;
        su*si];

    v = [radialSpeed*(cO*cu-sO*su*ci)-transverseSpeed*(cO*su+sO*cu*ci);
        radialSpeed*(sO*cu+cO*su*ci)-transverseSpeed*(sO*su-cO*cu*ci);
        radialSpeed*su*si+transverseSpeed*cu*si];
endfunction

// Сценарии из задания 1: №1, №3 и №65.
function ds = gravity(t, state)
    mu = 398600.4418;
    ds = [state(4:6); -mu*state(1:3)/norm(state(1:3))^3];
endfunction

mu = 398600.4418;
a = 7000;
n = sqrt(mu/a^3);
period = 2*%pi/n;
time = linspace(0, 2*period, 721);
theta = n*time;
iRef = 51.6*%pi/180;
raan = 30*%pi/180;
nuGridDeg = [0 90 180 270];
deltaIDeg = [0.01 0.1 0.5];
caseNumbers = [1 3 65];
root = get_absolute_file_path("hcw_phase_portraits.sce");
out = fullfile(root, "results");
if ~isdir(out) then
    mkdir(out);
end
fileId = mopen(fullfile(out, "grid_summary.csv"), "wt");
mfprintf(fileId, "case1,delta_i_deg,nu_deg,x0_km,y0_km,z0_km,vx0_km_s,vy0_km_s,vz0_km_s,max_error_km,max_separation_km,ratio,drift_km_s\n");

for k = 1:length(caseNumbers)
    for i = 1:length(deltaIDeg)
        for j = 1:length(nuGridDeg)
            nu = nuGridDeg(j)*%pi/180;
            deltaI = deltaIDeg(i)*%pi/180;
            [rRef, vRef] = keplerState(a, 0, iRef, raan, 0, nu, mu);
            [rTilt, vTilt] = keplerState(a, 0, iRef+deltaI, raan, 0, nu, mu);
            R = rRef/norm(rRef);
            N = cross3(rRef, vRef); N = N/norm(N);
            T = cross3(N, R);
            Q = [R T N];
            rho0 = Q'*(rTilt-rRef);
            vel0 = Q'*(vTilt-vRef)-cross3([0;0;n], rho0);
            // Сохраняем поперечные компоненты от разности наклонений,
            // а в плоскости задаём начальные условия из задания 1.
            rho0(2) = 0;
            if k == 1 then
                rho0(1) = 0.1; vel0(1) = 0; vel0(2) = -2*n*rho0(1);
            elseif k == 2 then
                rho0(1) = 0.1; vel0(1) = 0; vel0(2) = 0.1;
            else
                rho0(1) = 100; vel0(1) = 0.1; vel0(2) = 0.1;
            end
            x0 = rho0(1); y0 = rho0(2); z0 = rho0(3);
            vx0 = vel0(1); vy0 = vel0(2); vz0 = vel0(3);
            x = (4-3*cos(theta))*x0 + sin(theta)/n*vx0 + 2*(1-cos(theta))/n*vy0;
            y = y0 + 6*(sin(theta)-theta)*x0 - 2*(1-cos(theta))/n*vx0 + (4*sin(theta)-3*theta)/n*vy0;
            z = z0*cos(theta) + vz0/n*sin(theta);
            vy = 6*n*(cos(theta)-1)*x0 - 2*sin(theta)*vx0 + (4*cos(theta)-3)*vy0;
            vz = -n*sin(theta)*z0 + cos(theta)*vz0;
            ref0 = [rRef; vRef];
            dep0 = [rRef+Q*rho0; vRef+Q*(vel0+cross3([0;0;n],rho0))];
            ref = ode(ref0, 0, time, 1d-10, 1d-10, gravity);
            dep = ode(dep0, 0, time, 1d-10, 1d-10, gravity);
            actual = zeros(3,length(time));
            for m = 1:length(time)
                er = ref(1:3,m)/norm(ref(1:3,m));
                en = cross3(ref(1:3,m),ref(4:6,m)); en = en/norm(en);
                et = cross3(en,er);
                actual(:,m) = [er et en]'*(dep(1:3,m)-ref(1:3,m));
            end
            err = sqrt(sum(([x;y;z]-actual).^2,"r"));
            sep = sqrt(sum(actual.^2,"r"));
            C = 2*n*x0+vy0;
            mfprintf(fileId, "%d,%.2f,%d,%.6f,%.6f,%.6f,%.9f,%.9f,%.9f,%.6f,%.6f,%.8f,%.9f\n", caseNumbers(k),deltaIDeg(i),nuGridDeg(j),x0,y0,z0,vx0,vy0,vz0,max(err),max(sep),max(sep)/a,C);
            figure(10); clf(); gcf().figure_size = [1100 850]; gcf().background = color("white");
            subplot(2,2,1); plot(actual(1,:),actual(2,:),"b-"); plot(x,y,"r--"); xgrid();
            xtitle(msprintf("Случай %d, di=%.2f град, nu=%d град",caseNumbers(k),deltaIDeg(i),nuGridDeg(j)),"x, км","y, км"); legend(["Кеплер" "ХКУ"]);
            subplot(2,2,2); plot(y,vy); xgrid(); xtitle("Продольный фазовый портрет ХКУ","y, км","vy, км/с");
            subplot(2,2,3); plot(z,vz); xgrid(); xtitle("Поперечный фазовый портрет ХКУ","z, км","vz, км/с");
            subplot(2,2,4); plot(time/period,err); xgrid(); xtitle("Расхождение положений","Время, витки","Ошибка, км");
            xs2pdf(gcf(), fullfile(out, msprintf("case_%02d_i%02d_nu%03d.pdf",caseNumbers(k),i,nuGridDeg(j))));
        end
    end
end
mclose(fileId);
mprintf("Готово: %d сценария, %d наклонения, %d аномалии.\n",length(caseNumbers),length(deltaIDeg),length(nuGridDeg));
